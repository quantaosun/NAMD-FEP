"""The upload contract — what a user must provide, and how we check it.

The backend is the existing, proven pipeline (`acpype` -> `prepare_hybrid.py` ->
`build_system.py` -> `write_fep_inputs.py`). This module does not generalise it;
it defines the *file layout* that pipeline needs and refuses anything else with a
precise reason.

Why a contract instead of guessing: every hard failure this pipeline has
produced so far has been silent. `prepare_hybrid.py` picks the perturbation by
looking for "the `c3` carbon with exactly three hydrogens" -- on a methoxy or a
tert-butyl that selects the wrong atom and then runs to completion. `build_system.py`
hardcodes one protein's residue range. A validator that says "your ligand has two
candidate methyl groups" is worth more than any amount of clever inference.

Layout
------
    <name>/
      system.json          manifest (see DEFAULTS)
      protein.pdb          receptor
      ref.sdf | .mol2 | .pdb    reference ligand, WITH hydrogens
      mut.sdf | .mol2 | .pdb    mutant ligand,    WITH hydrogens

Manifest keys (all optional except the ligands are discovered by filename):
    segments    : "auto" | [[segid, first_res, last_res], ...]
                  "auto" splits on chain ID and CA-CA gap -- the thing that is
                  hardcoded to 148-412/416-483 today.
    temperature, padding, salt, n_windows   panel settings, defaulted.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

MANIFEST = "system.json"
PROTEIN = "protein.pdb"
LIGAND_STEMS = ("ref", "mut")
LIGAND_EXTS = (".sdf", ".mol2", ".pdb", ".mol")

DEFAULTS = {
    "segments": "auto",
    "temperature": 300.0,
    "padding": 15.0,          # solvate -t, matches build_system.py
    "salt": 0.0,              # autoionize neutralise-only
    "n_windows": 15,
}

# A ligand with fewer atoms than this is almost certainly a failed export.
MIN_LIGAND_ATOMS = 6
MIN_PROTEIN_ATOMS = 100
MIN_HYDROGENS = 1


@dataclass
class Problem:
    level: str          # "error" blocks the build; "warning" does not
    where: str
    message: str

    def __str__(self) -> str:
        return f"[{self.level}] {self.where}: {self.message}"


@dataclass
class Validation:
    ok: bool
    root: Path | None = None
    manifest: dict = field(default_factory=dict)
    ligands: dict[str, Path] = field(default_factory=dict)
    problems: list[Problem] = field(default_factory=list)

    @property
    def errors(self) -> list[Problem]:
        return [p for p in self.problems if p.level == "error"]

    @property
    def warnings(self) -> list[Problem]:
        return [p for p in self.problems if p.level == "warning"]

    def summary(self) -> str:
        if self.ok and not self.warnings:
            return "input looks valid"
        head = "valid, with warnings" if self.ok else "input rejected"
        return f"{head}: " + "; ".join(str(p) for p in self.problems)


# --------------------------------------------------------------------------
# light-weight structure readers (no rdkit dependency: this must run anywhere)
# --------------------------------------------------------------------------
def _count_pdb_atoms(path: Path) -> tuple[int, int, set[str]]:
    """(heavy_atom_count, hydrogen_count, chain_ids)."""
    heavy = hydro = 0
    chains: set[str] = set()
    with open(path, errors="replace") as fh:
        for ln in fh:
            if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                continue
            name = ln[12:16].strip() if len(ln) >= 16 else ""
            element = ln[76:78].strip() if len(ln) >= 78 else ""
            chains.add(ln[21:22].strip() if len(ln) > 21 else "")
            if element.upper() == "H" or (not element and re.match(r"^\d*H", name)):
                hydro += 1
            else:
                heavy += 1
    return heavy, hydro, chains


def _count_sdf_atoms(path: Path) -> tuple[int, int]:
    """(heavy, hydrogen) from the SDF counts line, with a title-block check."""
    text = path.read_text(errors="replace")
    if "$$$$" not in text and "V2000" not in text and "V3000" not in text:
        raise ValueError("not an SDF (no $$$$ or V2000/V3000 marker)")
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError("truncated SDF")
    counts = lines[3]
    try:
        natoms = int(counts[0:3])
    except ValueError:
        raise ValueError(f"unreadable SDF counts line: {counts!r}")
    heavy = hydro = 0
    for ln in lines[4:4 + natoms]:
        sym = ln[31:34].strip() if len(ln) >= 34 else ln.split()[-1]
        if sym.upper() == "H":
            hydro += 1
        else:
            heavy += 1
    return heavy, hydro


def _count_mol2_atoms(path: Path) -> tuple[int, int]:
    text = path.read_text(errors="replace")
    if "@<TRIPOS>ATOM" not in text:
        raise ValueError("not a MOL2 (no @<TRIPOS>ATOM)")
    heavy = hydro = 0
    started = False
    for ln in text.splitlines():
        if ln.startswith("@<TRIPOS>"):
            started = ln.strip() == "@<TRIPOS>ATOM"
            continue
        if not started or not ln.strip():
            continue
        parts = ln.split()
        sym = parts[5].split(".")[0] if len(parts) > 5 else ""
        if sym.upper() == "H":
            hydro += 1
        elif sym:
            heavy += 1
    return heavy, hydro


def read_ligand(path: Path) -> tuple[int, int]:
    """(heavy, hydrogen) atom counts for whichever format this is."""
    ext = path.suffix.lower()
    if ext == ".sdf" or ext == ".mol":
        return _count_sdf_atoms(path)
    if ext == ".mol2":
        return _count_mol2_atoms(path)
    if ext == ".pdb":
        heavy, hydro, _ = _count_pdb_atoms(path)
        return heavy, hydro
    raise ValueError(f"unsupported ligand format {ext!r}")


# --------------------------------------------------------------------------
def find_ligand(root: Path, stem: str) -> Path | None:
    for ext in LIGAND_EXTS:
        p = root / f"{stem}{ext}"
        if p.exists():
            return p
    return None


def validate(root: Path) -> Validation:
    """Check a candidate upload directory. Never raises."""
    root = Path(root)
    v = Validation(ok=False, root=root)

    if not root.is_dir():
        v.problems.append(Problem("error", str(root), "not a directory"))
        return v

    # -- manifest ------------------------------------------------------
    mpath = root / MANIFEST
    if not mpath.exists():
        v.problems.append(Problem(
            "error", MANIFEST,
            f"missing. Create it ({{}} is enough for defaults: "
            f"{json.dumps(DEFAULTS)})"))
        manifest = dict(DEFAULTS)
    else:
        try:
            manifest = json.loads(mpath.read_text())
        except json.JSONDecodeError as exc:
            v.problems.append(Problem("error", MANIFEST, f"invalid JSON: {exc}"))
            return v
        if not isinstance(manifest, dict):
            v.problems.append(Problem("error", MANIFEST, "must be a JSON object"))
            return v
        for key, default in DEFAULTS.items():
            manifest.setdefault(key, default)
        unknown = set(manifest) - set(DEFAULTS) - {"name", "ligand_a", "ligand_b"}
        if unknown:
            v.problems.append(Problem(
                "warning", MANIFEST,
                f"ignoring unrecognised key(s): {', '.join(sorted(unknown))}"))
    v.manifest = manifest

    # -- protein -------------------------------------------------------
    prot = root / PROTEIN
    if not prot.exists():
        v.problems.append(Problem("error", PROTEIN, "missing"))
    else:
        heavy, hydro, chains = _count_pdb_atoms(prot)
        if heavy < MIN_PROTEIN_ATOMS:
            v.problems.append(Problem(
                "error", PROTEIN,
                f"only {heavy} heavy atoms — is this really the receptor?"))
        if hydro == 0:
            # NOT an error. `build_system.py` feeds this to psfgen as
            # `segment { pdb ... }` plus `guesscoord`, so psfgen builds each
            # residue from the topology and adds the hydrogens itself -- the
            # 6I5I protein.pdb that produced the published result has zero
            # hydrogens. Still worth surfacing, because protonation *states*
            # are then decided entirely by build_system.py's pdbalias rules
            # (it assigns HSD to every histidine).
            v.problems.append(Problem(
                "warning", PROTEIN,
                "no hydrogens — psfgen will add them from the topology. "
                "Histidine protonation is then set by build_system.py, which "
                "assigns HSD to every HIS."))
        if len(chains) == 0:
            v.problems.append(Problem("error", PROTEIN, "no chain identifiers"))
        segs = manifest.get("segments")
        if segs == "auto" and len(chains) > 1:
            v.problems.append(Problem(
                "warning", PROTEIN,
                f"{len(chains)} chains detected — 'segments: auto' will treat "
                f"each as its own psfgen segment"))

    # -- ligands -------------------------------------------------------
    sizes: dict[str, tuple[int, int]] = {}   # only for ligands that parsed
    for stem in LIGAND_STEMS:
        p = find_ligand(root, stem)
        if p is None:
            v.problems.append(Problem(
                "error", f"{stem}.*",
                f"missing — need one of {', '.join(stem + e for e in LIGAND_EXTS)}"))
            continue
        v.ligands[stem] = p
        try:
            heavy, hydro = read_ligand(p)
        except (ValueError, OSError) as exc:
            # Record it and move on. Re-reading it later (as an earlier version
            # did) raised out of validate() and turned a clear message into a
            # traceback.
            v.problems.append(Problem("error", p.name, f"cannot parse: {exc}"))
            continue
        sizes[stem] = (heavy, hydro)
        if heavy < MIN_LIGAND_ATOMS:
            v.problems.append(Problem(
                "error", p.name,
                f"only {heavy} heavy atoms — the ligand looks truncated"))
        if hydro < MIN_HYDROGENS:
            v.problems.append(Problem(
                "error", p.name,
                "no hydrogens. acpype derives charges and atom types from an "
                "all-atom structure; a heavy-atom-only export gives wrong "
                "parameters, not an error."))

    # -- cross-checks (only over ligands that actually parsed) ---------
    if len(sizes) == 2 and sizes["ref"] == sizes["mut"]:
        v.problems.append(Problem(
            "warning", "ref/mut",
            "the two ligands have identical atom counts — if this is meant "
            "to be a real perturbation, check you uploaded the right files"))

    v.ok = not v.errors
    return v
