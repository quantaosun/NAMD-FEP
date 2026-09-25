"""Build the complex and solvent systems: psfgen -> solvate -> autoionize.

Pure VMD + NAMD workflow, no FEPrepare and no Maestro.

  complex  : protein + hybrid ligand -> solvate -> ionize
  solvent  : hybrid ligand alone     -> solvate -> ionize

Ported from `4YLJ/build_system.py`, which had already generalised the two things
that made the 6I5I copy system-specific:

* **The segment split is derived, not declared.** psfgen cannot carry a
  numbering gap inside one segment, which is why 6I5I's copy needed a
  hand-written `SEGMENTS = [("P1", 148, 412), ("P2", 416, 483)]`. Computing it
  from the numbering breaks means any structure works, and for 6I5I it
  reproduces that list exactly.
* **The box is written out, not read back from a constant.** `<leg>/box.json`
  is what `rbfe inputs` reads, so a system cannot end up simulated in another
  system's cell.

Everything else -- padding, salt, the pdbalias trio, the terminus types, the
binary paths -- comes from `system.ini`.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

from rbfe.config import Config
from rbfe.errors import BuildError

TIMEOUT_S = 900


def protein_segments(path: Path, chain: str | None) -> list[tuple[str, str, int, int]]:
    """One psfgen segment per run of continuous residue numbering.

    Returns (segment_name, chain, first_resseq, last_resseq).
    """
    res, seen = [], set()
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain and line[21] != chain:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        key = (line[21], resseq)
        if key not in seen:
            seen.add(key)
            res.append(key)
    if not res:
        raise BuildError(
            f"no protein atoms for chain {chain!r} in {path}. "
            f"Set [protein] chain in system.ini to the chain you want, or '' "
            f"to take every chain.")
    segs: list[list] = []
    for i, (ch, rs) in enumerate(res):
        if i == 0 or ch != res[i - 1][0] or rs != res[i - 1][1] + 1:
            segs.append([len(segs) + 1, ch, rs, rs])
        else:
            segs[-1][3] = rs
    return [(f"P{n}", ch, lo, hi) for n, ch, lo, hi in segs]


def run_vmd(vmd: str, script: Path, cwd: Path) -> str:
    try:
        r = subprocess.run([vmd, "-dispdev", "text", "-e", str(script)],
                           capture_output=True, text=True, cwd=str(cwd),
                           timeout=TIMEOUT_S)
    except FileNotFoundError:
        raise BuildError(
            f"VMD not found at {vmd!r}. Set [binaries] vmd in system.ini, or "
            f"export VMD=/path/to/vmd.") from None
    except subprocess.TimeoutExpired:
        raise BuildError(f"VMD timed out after {TIMEOUT_S}s running {script}") from None
    return r.stdout + "\n" + r.stderr


def split_protein(src: Path, outdir: Path, segs, chain: str | None) -> list[Path]:
    """Write one PDB per segment (residues in the given numbering range)."""
    paths = []
    lines = src.read_text().splitlines()
    for seg, ch, lo, hi in segs:
        out = [l for l in lines
               if l.startswith("ATOM") and l[21] == ch
               and lo <= int(l[22:26].strip()) <= hi]
        out.append("END")
        p = outdir / f"protein_{seg}.pdb"
        p.write_text("\n".join(out) + "\n")
        paths.append(p)
    return paths


def psfgen_complex(outdir: Path, segs, seg_paths, *, prot_rtf: Path, hybrid_rtf: Path,
                   hybrid_pdb: Path, pdbalias: list[str], first: str, last: str) -> Path:
    lines = ["package require psfgen",
             f"topology {prot_rtf}",
             f"topology {hybrid_rtf}"]
    lines += ["pdbalias " + a for a in pdbalias]
    for (seg, ch, lo, hi), p in zip(segs, seg_paths):
        lines += [f"segment {seg} {{ first {first}; last {last}; pdb {p} }}",
                  f"coordpdb {p} {seg}"]
    lines += [f"segment LIG {{ first none; last none; pdb {hybrid_pdb} }}",
              f"coordpdb {hybrid_pdb} LIG",
              "guesscoord",
              "writepsf complex.psf",
              "writepdb complex.pdb",
              "exit"]
    p = outdir / "psfgen.tcl"
    p.write_text("\n".join(lines) + "\n")
    return p


def psfgen_solvent(outdir: Path, *, hybrid_rtf: Path, hybrid_pdb: Path) -> Path:
    lines = ["package require psfgen",
             f"topology {hybrid_rtf}",
             f"segment LIG {{ first none; last none; pdb {hybrid_pdb} }}",
             f"coordpdb {hybrid_pdb} LIG",
             "guesscoord",
             "writepsf ligand.psf",
             "writepdb ligand.pdb",
             "exit"]
    p = outdir / "psfgen.tcl"
    p.write_text("\n".join(lines) + "\n")
    return p


def solvate_tcl(outdir: Path, psf: str, pdb: str, padding: float,
                neutralize: bool, salt: float) -> Path:
    ion_flags = "-neutralize" if neutralize else ""
    if salt > 0:
        ion_flags += f" -c {salt} -s KCL"
    lines = ["package require solvate",
             "package require autoionize",
             "mol delete all",
             f"mol load psf {psf} pdb {pdb}",
             f"solvate {psf} {pdb} -t {padding:g} -o solvated",
             f"autoionize -psf solvated.psf -pdb solvated.pdb {ion_flags} -o ionized".strip(),
             "mol delete all",
             "mol load psf ionized.psf pdb ionized.pdb",
             "set all [atomselect top all]",
             'puts "CELL [molinfo top get {a b c}]"',
             'puts "ORIGIN [molinfo top get {center}]"',
             "exit"]
    p = outdir / "solvate.tcl"
    p.write_text("\n".join(lines) + "\n")
    return p


def parse_box(out: str) -> dict:
    """Pull the box out of the tcl's CELL/ORIGIN lines (Angstrom).

    Tcl prints a list of numbers bare for `molinfo get {a b c}` but braced for
    the `center` vector, so accept either.
    """
    def grab(tag: str) -> list[float]:
        m = re.search(rf"^{tag}\s+(.*)$", out, re.M)
        if not m:
            raise BuildError(f"no {tag} in VMD output:\n" + out[-2000:])
        return [float(x) for x in m.group(1).strip().strip("{}").split()]
    return {"cell": grab("CELL"), "origin": grab("ORIGIN")}


def build(cfg: Config, out: Path, legs: list[str] | None = None) -> None:
    """Build every leg.  Writes <out>/<leg>/{psf,pdb,box.json}."""
    out = Path(out)
    legs = legs or cfg.legs
    vmd = str(cfg.resolve(cfg.get("binaries", "vmd")))
    prot = cfg.resolve(cfg.get("protein", "pdb"))
    chain = cfg.get("protein", "chain") or None
    first, last = cfg.get("protein", "first"), cfg.get("protein", "last")
    pdbalias = cfg.get("build", "pdbalias")
    padding = cfg.get("build", "padding")
    neutralize = cfg.get("build", "neutralize")
    salt = cfg.get("build", "salt")

    hybrid_rtf = out / "hybrid" / "hybrid.rtf"
    hybrid_pdb = out / "hybrid" / "hybrid.pdb"
    if not hybrid_rtf.exists() or not hybrid_pdb.exists():
        raise BuildError(
            f"{hybrid_rtf} / {hybrid_pdb} missing -- run `rbfe hybrid` first.")

    prot_rtf = cfg.resolve(cfg.get("topology", "prot_rtf"))
    segs = protein_segments(prot, chain)
    nres = sum(hi - lo + 1 for _, _, lo, hi in segs)
    print(f"protein  : {prot}  chain(s) {sorted({c for _, c, _, _ in segs})}  "
          f"{nres} residues in {len(segs)} segment(s)")
    print(f"segments : {[(s, lo, hi) for s, _, lo, hi in segs]}"
          f"   (derived from numbering breaks, not declared)")

    for leg in legs:
        d = out / leg
        d.mkdir(parents=True, exist_ok=True)
        print(f"\n===== {leg.upper()} leg =====")

        if leg == "complex":
            seg_paths = split_protein(prot, d, segs, chain)
            script = psfgen_complex(d, segs, seg_paths, prot_rtf=prot_rtf,
                                    hybrid_rtf=hybrid_rtf, hybrid_pdb=hybrid_pdb,
                                    pdbalias=pdbalias, first=first, last=last)
            psf, pdb = "complex.psf", "complex.pdb"
        else:
            script = psfgen_solvent(d, hybrid_rtf=hybrid_rtf, hybrid_pdb=hybrid_pdb)
            psf, pdb = "ligand.psf", "ligand.pdb"

        vout = run_vmd(vmd, script, d)
        if not (d / psf).exists():
            print(vout[-3000:])
            raise BuildError(f"psfgen failed for {leg} (no {psf})")
        print(f"  psfgen OK -> {psf}")

        script = solvate_tcl(d, psf, pdb, padding, neutralize, salt)
        vout = run_vmd(vmd, script, d)
        if not (d / "ionized.psf").exists():
            print(vout[-3000:])
            raise BuildError(f"solvate/autoionize failed for {leg}")
        box = parse_box(vout)
        (d / "box.json").write_text(json.dumps(box, indent=2) + "\n")
        print("  solvate/ionize OK -> ionized.psf/pdb")
        print(f"    cell   {[round(x, 3) for x in box['cell']]}")
        print(f"    origin {[round(x, 3) for x in box['origin']]}  -> box.json "
              f"(rbfe inputs reads this; the cell is never hardcoded)")
