"""Turn a PDB entry into the inputs `rbfe hybrid` needs.

Everything downstream of `rbfe hybrid` assumes the inputs already exist:
`systems/4ylj/inputs/` and `systems/6i5i/inputs/` are frozen trees that someone
prepared by hand.  This module is that missing first step, and exists because
preparing a new system by hand is where the quiet mistakes live -- an element
column read from the wrong field, a hydrogen count that never got checked, a
ligand whose bond orders were perceived as something plausible but wrong.

    rbfe-prep report   3HTB.pdb                  what is in here?
    rbfe-prep fetch    3HTB                      download it
    rbfe-prep strip    3HTB.pdb -o inputs/protein.pdb --chain A
    rbfe-prep ligand   3HTB.pdb JZ4 -o inputs/mut.sdf
    rbfe-prep ligand   3HTB.pdb JZ4 -o inputs/ref.sdf --drop C4
    rbfe-prep params   inputs/ref.sdf --name ref -o inputs/
    rbfe-prep scaffold systems/3htb --chain A --resname UNL

What it deliberately does NOT do
--------------------------------
It does not choose the mutation.  Which ligand is the reference, which way the
transformation runs, and whether the pair is even suitable for a relative free
energy are chemistry judgements; `--drop` covers the commonest edit (removing a
substituent) and anything beyond that is RDKit's job, not this tool's.  What
this removes is the file shuffling around that judgement, and the chances to get
it silently wrong.

Chemistry is perceived by **obabel**, not by this module.  `acpype` already
requires obabel, so it is not a new dependency, and it was checked against a
hand-built reference: for 3HTB's ligand both routes give the same 22 atoms, the
same 12 hydrogens and the same SMILES.  Perception is still perception, so every
step prints what it produced -- formula, SMILES, atom count -- for you to
disagree with, and `ligand --smiles` overrides it when it is wrong.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

from rbfe.errors import RbfeError

RCSB = "https://files.rcsb.org/download/{id}.pdb"
REPO = Path(__file__).resolve().parent.parent

# Solvent and common crystallisation additives.  Not a judgement call: these are
# stripped by default and `--keep` puts any of them back.
SOLVENT = {
    "HOH", "DOD", "WAT", "H2O",          # water
    "PO4", "SO4", "NO3", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CU",
    "NI", "CO", "CD", "BR", "IOD", "F",  # ions
    "GOL", "EDO", "PEG", "MPD", "DMS", "TRS", "EPE", "MES", "BME", "DTT",
    "ACT", "ACY", "IMD", "CIT", "TLA",  # cryo/buffer
}

# Bond orders obabel cannot always get right, where a wrong answer is silent.
# Reported rather than policed -- they are a reason to look, not to refuse.
WATCH = {"MSE": "selenomethionine: psfgen needs it built as MET + SE",
         "SEP": "phosphoserine: not a CHARMM36 residue",
         "TPO": "phosphothreonine: not a CHARMM36 residue",
         "CSO": "modified cysteine",
         "KCX": "modified lysine"}


class PrepError(RbfeError):
    """Preparation cannot proceed, or produced something that cannot be right."""


# --------------------------------------------------------------------------
# external tools
# --------------------------------------------------------------------------

def _tool(name: str, env: str) -> str:
    """Locate `name`, honouring `$env`, PATH, then the AI Studio layout."""
    override = os.environ.get(env)
    if override:
        if not Path(override).exists():
            raise PrepError(f"${env} points at {override}, which does not exist")
        return override
    found = shutil.which(name)
    if found:
        return found
    fallback = Path.home() / "external-libraries" / "bin" / name
    if fallback.exists():
        return str(fallback)
    raise PrepError(
        f"{name} not found. Install it, put it on PATH, or set ${env}.\n"
        f"  On this box it lives in ~/external-libraries/bin, which is not on "
        f"PATH by default.")


def _tool_env(tool: str) -> dict:
    """Environment for running a tool, with its sibling binaries reachable.

    acpype shells out to obabel, and on this box both sit in the same
    non-PATH directory, so the directory is prepended rather than requiring the
    caller to have exported it.
    """
    env = dict(os.environ)
    d = str(Path(tool).resolve().parent)
    env["PATH"] = d + os.pathsep + env.get("PATH", "")
    return env


# --------------------------------------------------------------------------
# PDB reading
# --------------------------------------------------------------------------

def _read(path: Path) -> list[str]:
    if not Path(path).exists():
        raise PrepError(f"{path} does not exist")
    return Path(path).read_text().splitlines()


def _altloc(line: str) -> str:
    return line[16] if len(line) > 16 else " "


def _element(line: str) -> str:
    """The element for a PDB line, from columns 77-78 with a name fallback.

    The column is authoritative when it is populated, because an atom name
    cannot always be read as one (`CA` is calcium in a HETATM record and an
    alpha carbon in an ATOM record).  When the column is blank -- which older
    files and some writers leave -- fall back to the name.
    """
    col = line[76:78].strip() if len(line) >= 78 else ""
    if col:
        return col.upper()
    name = line[12:16].strip()
    m = re.match(r"([A-Za-z]{1,2})", name)
    return m.group(1).upper() if m else ""


def inventory(path: Path) -> dict:
    """What a PDB actually contains, before anything is stripped or kept."""
    lines = _read(path)
    chains: dict[str, int] = {}
    het: dict[str, dict] = {}
    altlocs: dict[str, int] = {}
    residues: dict[tuple[str, str], set] = {}
    for l in lines:
        if l.startswith(("ATOM  ", "HETATM")):
            # Only coordinate records have an altLoc column.  Counting column 16
            # on every line would read the header records too and report a
            # scattering of letters as "alternate conformations".
            a = _altloc(l).strip()
            if a:
                altlocs[a] = altlocs.get(a, 0) + 1
        if l.startswith("ATOM  "):
            ch = l[21]
            chains[ch] = chains.get(ch, 0) + 1
            residues.setdefault((ch, l[22:27].strip()), set()).add(l[17:20].strip())
        elif l.startswith("HETATM"):
            res = l[17:20].strip()
            e = het.setdefault(res, {"atoms": 0, "chains": set(), "elements": set()})
            e["atoms"] += 1
            e["chains"].add(l[21])
            e["elements"].add(_element(l))
    gaps = {}
    for ch in chains:
        nums = sorted(int(n) for (c, n), r in residues.items()
                      if c == ch and n.lstrip("-").isdigit())
        if not nums:
            continue
        gaps[ch] = [n for n in range(nums[0], nums[-1] + 1) if n not in nums]
    return {"chains": chains, "het": het, "altlocs": altlocs, "gaps": gaps,
            "residues": {c: sorted({r for (cc, _), rs in residues.items()
                                    if cc == c for r in rs})
                         for c in chains}}


def cmd_report(a) -> int:
    inv = inventory(Path(a.pdb))
    print(f"{a.pdb}")
    print(f"  protein   : {len(inv['chains'])} chain(s)")
    for ch, n in sorted(inv["chains"].items()):
        gap = inv["gaps"].get(ch, [])
        probs = [r for r in inv["residues"][ch] if r in WATCH]
        print(f"      chain {ch!r}: {n} atoms, " +
              (f"{len(gap)} missing residue number(s) {gap[:6]}"
               if gap else "no numbering gaps"))
        for r in probs:
            print(f"        WARNING {r}: {WATCH[r]}")
    if inv["het"]:
        print("  HETATM    :")
        for res, e in sorted(inv["het"].items(), key=lambda kv: -kv[1]["atoms"]):
            mark = "  (stripped as solvent/ion)" if res in SOLVENT else ""
            print(f"      {res:<4} {e['atoms']:>5} atoms  chain(s) "
                  f"{''.join(sorted(e['chains']))}  "
                  f"{','.join(sorted(e['elements']))}{mark}")
    else:
        print("  HETATM    : none")
    if inv["altlocs"]:
        print(f"  altlocs   : {inv['altlocs']}  (strip keeps blank + one, "
              f"default A)")
    return 0


def cmd_fetch(a) -> int:
    dest = Path(a.out or f"{a.pdb_id.upper()}.pdb")
    if dest.exists() and not a.force:
        raise PrepError(f"{dest} already exists (use --force to overwrite)")
    url = RCSB.format(id=a.pdb_id.upper())
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read()
    except (urllib.error.URLError, OSError) as exc:
        raise PrepError(
            f"could not fetch {url}: {exc}\n"
            f"  Download it by hand and pass the file to the other subcommands; "
            f"everything except `fetch` works on a local PDB.") from None
    if not data.startswith((b"HEADER", b"ATOM", b"HETATM", b"EXPDTA", b"TITLE")):
        raise PrepError(f"{url} did not return a PDB file "
                        f"(got {data[:40]!r}) -- is the ID right?")
    dest.write_bytes(data)
    print(f"wrote {dest} ({len(data.splitlines())} lines) from {url}")
    return 0


def cmd_strip(a) -> int:
    """Protein only: drop every HETATM, resolve altlocs, optionally one chain."""
    lines = _read(Path(a.pdb))
    out, dropped_alt, dropped_res, residues = [], 0, 0, set()
    for l in lines:
        if not l.startswith("ATOM  "):
            if l.startswith("HETATM"):
                dropped_res += 1
            continue
        if a.chain and l[21] not in a.chain:
            continue
        alt = _altloc(l)
        if alt.strip() and alt != a.altloc:
            dropped_alt += 1
            continue
        # normalise the altloc to blank so psfgen sees one conformer
        out.append(l[:16] + " " + l[17:])
        residues.add((l[21], l[22:27].strip()))
    if not out:
        raise PrepError(
            f"no protein atoms selected from {a.pdb}"
            + (f" for chain(s) {a.chain}" if a.chain else "")
            + ".\n  Run `rbfe-prep report` to see what the file contains.")
    dest = Path(a.out or "protein.pdb")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(out) + "\nEND\n")
    print(f"wrote {dest}")
    print(f"  {len(out)} atoms, {len(residues)} residues, "
          f"chain(s) {sorted({c for c, _ in residues})}")
    print(f"  dropped {dropped_res} HETATM atom(s) and {dropped_alt} "
          f"altloc-{a.altloc} atom(s)")
    return 0


def cmd_ligand(a) -> int:
    """One HETATM residue -> a hydrogen-complete SDF in the crystal frame."""
    if a.smiles and a.drop:
        raise PrepError("--smiles and --drop are different ways to say what the "
                        "molecule is; use one or the other")
    lines = _read(Path(a.pdb))
    sel = [l for l in lines
           if l.startswith("HETATM") and l[17:20].strip() == a.resname]
    if not sel:
        raise PrepError(
            f"no HETATM residue {a.resname!r} in {a.pdb}.\n"
            f"  `rbfe-prep report` lists what is there.")
    if a.chain:
        sel = [l for l in sel if l[21] in a.chain]
        if not sel:
            raise PrepError(f"residue {a.resname} is not in chain(s) {a.chain}")

    dropped = []
    if a.drop:
        names = {l[12:16].strip() for l in sel}
        unknown = [n for n in a.drop if n not in names]
        if unknown:
            raise PrepError(
                f"--drop names atom(s) not in {a.resname}: {', '.join(unknown)}\n"
                f"  available: {', '.join(sorted(names))}")
        keep = [l for l in sel if l[12:16].strip() not in set(a.drop)]
        dropped = list(a.drop)
        sel = keep
        if not sel:
            raise PrepError("--drop removed every atom")

    with tempfile.TemporaryDirectory(prefix="rbfe-prep-") as tmp:
        src = Path(tmp) / f"{a.resname}.pdb"
        src.write_text("\n".join(sel) + "\nEND\n")
        dest = Path(a.out or f"{a.resname}.sdf")
        dest.parent.mkdir(parents=True, exist_ok=True)
        if a.smiles:
            _from_smiles(a.smiles, src, dest)
        else:
            obabel = _tool("obabel", "OBABEL")
            # -h completes the valences with hydrogens.  The crystal has none,
            # and acpype's AM1-BCC charges are meaningless without them.
            r = subprocess.run([obabel, str(src), "-O", str(dest), "-h"],
                               capture_output=True, text=True,
                               env=_tool_env(obabel))
            if r.returncode != 0 or not dest.exists():
                raise PrepError(
                    f"obabel could not convert {src.name}:\n"
                    f"{(r.stderr or r.stdout)[-800:]}\n"
                    f"  If the bond orders it perceived are wrong, pass "
                    f"--smiles.")

    print(f"wrote {dest}")
    _describe_sdf(dest, dropped)
    return 0


def _from_smiles(smiles: str, src: Path, dest: Path) -> None:
    """Assign bond orders from a supplied SMILES, keeping the crystal frame."""
    try:
        from rdkit import Chem
    except ImportError:
        raise PrepError(
            "--smiles needs RDKit, which is not installed. Install it, or omit "
            "--smiles and let obabel perceive the bonds.") from None
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise PrepError(f"RDKit could not parse --smiles {smiles!r}")
    probe = Chem.MolFromPDBFile(str(src), removeHs=False, sanitize=False)
    if probe is None:
        raise PrepError(f"RDKit could not read {src}")
    # skeletonise both so the match is by element and connectivity only
    def skel(m):
        m = Chem.RWMol(m)
        for b in m.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for at in m.GetAtoms():
            at.SetIsAromatic(False)
            at.SetNoImplicit(True)
        out = m.GetMol()
        out.UpdatePropertyCache(strict=False)
        return out
    match = skel(probe).GetSubstructMatch(skel(template))
    if len(match) != template.GetNumAtoms():
        raise PrepError(
            f"--smiles {smiles!r} does not match the atoms in {src.name} "
            f"(matched {len(match)} of {template.GetNumAtoms()}). "
            f"Check the residue name and --drop.")
    conf = probe.GetConformer()
    crystal = {}
    for atom, idx in zip(template.GetAtoms(), match):
        p = conf.GetAtomPosition(idx)
        crystal[atom.GetIdx()] = (p.x, p.y, p.z)
    mol = Chem.AddHs(template)
    out = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        if i in crystal:
            out.SetAtomPosition(i, crystal[i])
    mol.AddConformer(out, assignId=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    Chem.MolToMolFile(mol, str(dest))


def _describe_sdf(path: Path, dropped: list[str]) -> None:
    """Say what the file actually is, so a wrong perception is visible."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors as D
    except ImportError:
        print("  (install RDKit to have the formula and SMILES checked here)")
        return
    mol = Chem.MolFromMolFile(str(path), removeHs=False)
    if mol is None:
        raise PrepError(f"RDKit could not read back {path} -- it is likely "
                        f"malformed; do not use it")
    n_h = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    print(f"  {mol.GetNumAtoms()} atoms ({n_h} explicit H), "
          f"{D.CalcMolFormula(mol)}")
    print(f"  SMILES {Chem.MolToSmiles(Chem.RemoveHs(mol))}")
    if dropped:
        print(f"  derived by deleting {', '.join(dropped)}")
    print("  ^ CHECK THIS. These are perceived, and ACYPYPE's charges depend on "
          "them.")


# --------------------------------------------------------------------------
# parameterisation
# --------------------------------------------------------------------------

def cmd_params(a) -> int:
    sdf = Path(a.sdf)
    if not sdf.exists():
        raise PrepError(f"{sdf} does not exist")
    out = Path(a.out or sdf.parent)
    out.mkdir(parents=True, exist_ok=True)
    acpype = _tool("acpype", "ACPYPE")
    env = _tool_env(acpype)
    _tool("obabel", "OBABEL")          # acpype needs it; fail early and clearly

    with tempfile.TemporaryDirectory(prefix="rbfe-prep-") as tmp:
        work = Path(tmp) / sdf.name
        work.write_bytes(sdf.read_bytes())
        cmd = [acpype, "-i", work.name, "-c", a.charge_method, "-n", str(a.net),
               "-a", a.forcefield, "-b", a.name, "-o", "charmm"]
        r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True,
                           env=env, timeout=a.timeout)
        log = (r.stdout or "") + (r.stderr or "")
        sub = Path(tmp) / f"{a.name}.acpype"
        if r.returncode != 0 or not sub.exists():
            tail = "\n".join(log.strip().splitlines()[-15:])
            raise PrepError(f"acpype failed (exit {r.returncode}):\n{tail}")

        produced = {}
        for src, dst in ((sub / f"{a.name}_CHARMM.rtf", out / f"{a.name}.rtf"),
                         (sub / f"{a.name}_CHARMM.prm", out / f"{a.name}.prm"),
                         (sub / f"{a.name}_bcc_gaff2.mol2", out / f"{a.name}.mol2")):
            if not src.exists():
                raise PrepError(
                    f"acpype did not write {src.name}. Its output was:\n"
                    f"{', '.join(sorted(p.name for p in sub.iterdir()))}")
            shutil.copy(src, dst)
            produced[dst] = src
        n = _mol2_to_pdb(produced[out / f"{a.name}.mol2"], out / f"{a.name}.pdb")
        _check_names(out / f"{a.name}.rtf", out / f"{a.name}.pdb")

    print(f"wrote {out}/{a.name}.{{rtf,prm,pdb,mol2}}  ({n} atoms)")
    if a.keep_log:
        Path(a.keep_log).write_text(log)
        print(f"wrote {a.keep_log}")
    if re.search(r"charge to be balanced: total ([\-0-9.]+)", log):
        drift = float(re.search(r"charge to be balanced: total ([\-0-9.]+)", log).group(1))
        if abs(drift) > 0.01:
            print(f"  note: acpype balanced a charge drift of {drift:+.4f} e")
    total = _rtf_charge(out / f"{a.name}.rtf")
    print(f"  RESI total charge {total:+.4f} e  (rbfe checks this is near an integer)")
    return 0


def _mol2_to_pdb(mol2: Path, dest: Path) -> int:
    """A ligand PDB whose atom names are the .rtf's -- what `load_ligand` reads.

    The coordinates are the ones acpype carried through from the input SDF, so
    the ligand keeps the pose it was extracted in.
    """
    lines = mol2.read_text().splitlines()
    try:
        i = lines.index("@<TRIPOS>ATOM") + 1
    except ValueError:
        raise PrepError(f"{mol2} has no @<TRIPOS>ATOM block") from None
    rows, n = [], 0
    while i < len(lines) and not lines[i].startswith("@<TRIPOS>"):
        f = lines[i].split()
        if len(f) >= 6:
            n += 1
            name, x, y, z = f[1], float(f[2]), float(f[3]), float(f[4])
            elem = re.match(r"[A-Za-z]+", name).group(0)[:2]
            rows.append(f"HETATM{n:5d} {name:<4s} UNL Z   1    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}"
                        f"          {elem:>2s}")
        i += 1
    if not rows:
        raise PrepError(f"{mol2} has an empty atom block")
    dest.write_text("\n".join(rows) + "\nEND\n")
    return n


def _rtf_charge(rtf: Path) -> float:
    tot = 0.0
    for m in re.finditer(r"^ATOM\s+\S+\s+\S+\s+([-\d.]+)", rtf.read_text(), re.M):
        tot += float(m.group(1))
    return tot


def _check_names(rtf: Path, pdb: Path) -> None:
    """Every atom the .rtf declares must have a coordinate -- or the build fails."""
    rtf_names = set(re.findall(r"^ATOM\s+(\S+)", rtf.read_text(), re.M))
    pdb_names = {l[12:16].strip() for l in pdb.read_text().splitlines()
                 if l.startswith("HETATM")}
    missing = sorted(rtf_names - pdb_names)
    if missing:
        raise PrepError(
            f"{pdb.name} is missing {len(missing)} atom(s) the .rtf declares: "
            f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}\n"
            f"  `rbfe hybrid` would fail on this; the .pdb is not usable.")


# --------------------------------------------------------------------------
# scaffolding
# --------------------------------------------------------------------------

def cmd_scaffold(a) -> int:
    d = Path(a.dir).resolve()
    ini = d / "system.ini"
    if ini.exists() and not a.force:
        raise PrepError(f"{ini} already exists (use --force to overwrite)")
    (d / "inputs").mkdir(parents=True, exist_ok=True)
    toppar = os.path.relpath(REPO / "toppar", d)
    ini.write_text(_config(d.name, a, toppar))
    print(f"wrote {ini}")
    print(f"""
The config is deliberately incomplete: [mutation] strategy is blank, and `rbfe`
refuses to run until it is set -- it will not guess at chemistry.

Next, and this part is yours, not the tool's:
  1. decide which ligand is the reference and which the mutant
  2. EDIT {ini} and set [mutation] strategy (element_swap / atom_addition /
     mcs) and, for the first two, vanish/appear
  3. cd {d.name} && rbfe similarity   # is this pair suitable for RBFE at all?
     cd {d.name} && rbfe hybrid""")
    return 0


def _config(name: str, a, toppar: str) -> str:
    return f"""\
# One system, one config. `rbfe-prep scaffold` wrote this; the [mutation]
# block is deliberately blank because choosing it is a chemistry decision.

[system]
name        = {name}
description = {a.description or ''}

[ligand]
dir     = inputs
ref     = {a.ref}
mut     = {a.mut}
resname = {a.resname}

[protein]
pdb            = inputs/protein.pdb
chain          = {a.chain or 'A'}
segment_prefix = P
first          = NTER
last           = CTER

[mutation]
# What changes between the two ligands. An unrecognised strategy is a hard
# error -- the engine does not guess at chemistry.
#
#   element_swap  -- one atom becomes a different element (I -> Br)
#   atom_addition -- the mutant has an atom the reference lacks (N-CH3 -> N-H)
#   mcs           -- anything else; the mapping is DERIVED from the structures,
#                    so omit vanish/appear entirely
strategy      =
vanish        =
appear        =
placement     = native
missing_terms = synth_from_geometry

[hybrid]
title       = {a.description or name}
type_prefix = L
rtf_order   = sorted
prm_order   = source_then_new

[build]
padding    = 15.0
neutralize = yes
salt       = 0.0
pdbalias   = atom ILE CD1 CD
             atom SER HG HG1
             residue HIS HSD

[topology]
prot_rtf = {{root}}/{toppar}/top_all36_prot.rtf
params   = {{root}}/{toppar}/par_all36m_prot.prm
           {{root}}/{toppar}/par_all36_na.prm
           {{root}}/{toppar}/par_all36_carb.prm
           {{root}}/{toppar}/par_all36_lipid.prm
           {{root}}/{toppar}/par_all36_cgenff.prm
           {{root}}/{toppar}/par_water_ions_clean.prm
           {{hybrid}}

[run]
legs              = complex solvent
lambdas           = 0.0 0.045 0.09 0.14546 0.22425 0.30303 0.38182 0.46061
                    0.5394 0.61819 0.697 0.77576 0.85455 0.91 0.955 1.0
steps_per_window  = 500000
equil_steps       = 50000
alch_equil_steps  = 50000
outfreq           = 500
min_steps         = 5000
timestep          = 2.0
temperature       = 300.0
cutoff            = 12.0
switchdist        = 10.0
pairlistdist      = 14.0
elec_lambda_start = 0.5
vdw_lambda_end    = 1.0
vdw_shift_coeff   = 5.0
decouple          = off
max_hours         = 40.0

[binaries]
namd         = ${{NAMD:-/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3}}
vmd          = ${{VMD:-/home/aistudio/vmd-env/bin/vmd}}
flags        = +p1 +devices 0
gpu_resident = on
"""


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="rbfe-prep", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add(name, fn, help_):
        p = sub.add_parser(name, help=help_)
        p.set_defaults(fn=fn)
        return p

    p = add("report", cmd_report, "what is in a PDB, before changing anything")
    p.add_argument("pdb")

    p = add("fetch", cmd_fetch, "download a PDB entry from RCSB")
    p.add_argument("pdb_id")
    p.add_argument("-o", "--out")
    p.add_argument("--force", action="store_true")

    p = add("strip", cmd_strip, "protein only: drop HETATM, resolve altlocs")
    p.add_argument("pdb")
    p.add_argument("-o", "--out", help="default protein.pdb")
    p.add_argument("--chain", help="keep only this chain (default: all)")
    p.add_argument("--altloc", default="A",
                   help="alternate conformation to keep (default A)")

    p = add("ligand", cmd_ligand, "one HETATM residue -> a hydrogen-complete SDF")
    p.add_argument("pdb")
    p.add_argument("resname")
    p.add_argument("-o", "--out")
    p.add_argument("--chain")
    p.add_argument("--drop", nargs="+", metavar="ATOM",
                   help="delete these atoms first, to derive a variant")
    p.add_argument("--smiles", help="assign bond orders from this SMILES "
                                    "instead of letting obabel perceive them")

    p = add("params", cmd_params, "acpype -> ref.{rtf,prm,pdb,mol2}")
    p.add_argument("sdf")
    p.add_argument("--name", default="ref", help="output stem (default ref)")
    p.add_argument("-o", "--out", help="default: alongside the SDF")
    p.add_argument("--net", type=int, default=0, help="net charge (default 0)")
    p.add_argument("--forcefield", default="gaff2")
    p.add_argument("--charge-method", default="bcc")
    p.add_argument("--timeout", type=int, default=3600)
    p.add_argument("--keep-log", help="also write acpype's log here")

    p = add("scaffold", cmd_scaffold, "write a skeleton system.ini")
    p.add_argument("dir")
    p.add_argument("--chain")
    p.add_argument("--resname", default="UNL")
    p.add_argument("--ref", default="ref")
    p.add_argument("--mut", default="mut")
    p.add_argument("--description")
    p.add_argument("--force", action="store_true")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.fn(args)
    except RbfeError as e:
        print(str(e), file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
