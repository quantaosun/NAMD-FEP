"""CHARMM/GAFF file parsing: .rtf, .prm, .pdb, and atom-type helpers.

Ported from `4YLJ/prepare_hybrid.py` (itself already generic) with two fixes,
both of which are about not being quietly wrong:

1. `element()` used to return the first letter of the GAFF type uppercased.
   That is correct for `c3` -> C, `hn` -> H, `i` -> I, and wrong for every
   two-letter element: `br` -> B (boron), `cl` -> C (carbon).  It is derived
   from the MASS instead, which the .rtf always carries, and the fallback
   parses the atom NAME.  Ambiguity is a hard error rather than a coin flip.

2. `_nt()` lower-cases types and maps X to itself, but the wildcard check was
   `t.upper() == "X"` on a single letter -- a type like `x1` would slip through.
   Wildcards are matched on the whole token.
"""

from __future__ import annotations

import math
from pathlib import Path

from rbfe.errors import ChemistryError
from rbfe.model import Atom, Ligand, XYZ

# Standard atomic weights.  Only used to turn a MASS line into an element
# symbol, so precision beyond the first decimal is irrelevant -- what matters is
# that no two entries in this table are close enough to confuse.
ELEMENT_MASSES: dict[str, float] = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38, "Se": 78.971, "Br": 79.904, "Mo": 95.95, "Ag": 107.868,
    "I": 126.904, "Ba": 137.327, "Pt": 195.084, "Au": 196.967, "Hg": 200.592,
    "Pb": 207.2,
}

# Two candidate elements count as ambiguous when their masses are this close to
# equidistant from the queried mass.  The tightest real pair in the table is
# Co/Ni (58.933/58.693, 0.24 apart) -- neither is a plausible ligand type, and
# 0.2 Da keeps them apart while still catching a genuinely unknown mass.
AMBIGUITY_DA = 0.2


def element_from_mass(mass: float, typ: str) -> str | None:
    """The element whose standard atomic weight is nearest `mass`.

    Returns None when the nearest two candidates are within `AMBIGUITY_DA` of
    each other -- the caller then falls back to the atom name, or fails.
    """
    if not mass or mass <= 0:
        return None
    ranked = sorted(ELEMENT_MASSES.items(), key=lambda kv: abs(kv[1] - mass))
    (e1, m1), (e2, m2) = ranked[0], ranked[1]
    if abs(abs(m2 - mass) - abs(m1 - mass)) < AMBIGUITY_DA:
        return None
    return e1


def element_from_name(name: str) -> str | None:
    """Element symbol from a PDB atom name: `Br1` -> Br, `C9` -> C, `I1` -> I.

    Longest prefix wins, and a two-letter symbol is only accepted when what
    follows is a digit or nothing -- so `Br1` is bromine while `C9` is carbon,
    not a nonexistent "C9" element.

    This is the *fallback*.  The mass is the authority: a ligand `.rtf` from
    acpype always carries MASS lines, and without them an atom named `CA` is
    genuinely ambiguous between calcium and an alpha carbon.  When the mass is
    available this function is never reached.
    """
    s = name.strip()
    if not s:
        return None
    if len(s) >= 2:
        cand = s[0].upper() + s[1].lower()
        if cand in ELEMENT_MASSES and (len(s) == 2 or s[2].isdigit()):
            return cand
    head = s[0].upper()
    if head in ELEMENT_MASSES and (len(s) == 1 or not s[1].isalpha()):
        return head
    return None


def element(atom: Atom) -> str:
    """The element of one atom, from its mass, falling back to its name.

    Fails loudly rather than guessing: a wrong element writes a wrong PDB
    element column, and psfgen then builds the wrong atom.
    """
    by_mass = element_from_mass(atom.mass, atom.type)
    if by_mass:
        return by_mass
    by_name = element_from_name(atom.name)
    if by_name:
        return by_name
    raise ChemistryError(
        f"cannot determine the element of atom '{atom.name}' "
        f"(type '{atom.type}', mass {atom.mass}).\n"
        f"  Neither the mass nor the name identifies an element unambiguously. "
        f"Fix the input .rtf MASS line or rename the atom.")


# ---------------------------------------------------------------------------
# .rtf / .prm / .pdb
# ---------------------------------------------------------------------------

def _nt(t: str) -> str:
    """Normalise a parameter-table key: wildcards stay X, everything else lower."""
    return "X" if t == "X" else t.lower()


def is_wildcard(t: str) -> bool:
    return t == "X"


def parse_rtf(path: Path) -> dict:
    mass: dict[str, float] = {}
    atoms: list[tuple[str, str, float]] = []
    bonds, angles, dihe, impr = [], [], [], []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        p = s.split()
        kw = p[0].upper()
        if kw == "MASS" and len(p) >= 4:
            mass[p[2].lower()] = float(p[3])
        elif kw == "ATOM" and len(p) >= 4:
            atoms.append((p[1], p[2].lower(), float(p[3])))
        elif kw == "BOND" and len(p) >= 3:
            bonds.append((p[1], p[2]))
        elif kw == "ANGL" and len(p) >= 4:
            angles.append((p[1], p[2], p[3]))
        elif kw == "DIHE" and len(p) >= 5:
            dihe.append((p[1], p[2], p[3], p[4]))
        elif kw in ("IMPR", "IMPH") and len(p) >= 5:
            impr.append((p[1], p[2], p[3], p[4]))
    return {"mass": mass, "atoms": atoms, "bonds": bonds,
            "angles": angles, "dihedrals": dihe, "impropers": impr}


def parse_prm(path: Path) -> dict:
    bonds, angles, dihe, impr, nonb = {}, {}, {}, {}, {}
    mode = None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith(("*", "!")):
            continue
        u = s.upper()
        if u.startswith("BOND"):
            mode = "bond"; continue
        if u.startswith("ANGL"):
            mode = "angl"; continue
        if u.startswith(("DIHE", "DIHEDRAL")):
            mode = "dihe"; continue
        if u.startswith(("IMPR", "IMPHI")):
            mode = "impr"; continue
        if u.startswith("NONBONDED"):
            mode = "nonb"; continue
        if u.startswith(("CUTNB", "HBOND", "END", "NBFIX")):
            continue
        p = s.split()
        if mode == "bond" and len(p) >= 4:
            bonds[(_nt(p[0]), _nt(p[1]))] = (float(p[2]), float(p[3]))
        elif mode == "angl" and len(p) >= 5:
            angles[(_nt(p[0]), _nt(p[1]), _nt(p[2]))] = (float(p[3]), float(p[4]))
        elif mode == "dihe" and len(p) >= 7:
            dihe[tuple(_nt(x) for x in p[:4])] = (float(p[4]), int(p[5]), float(p[6]))
        elif mode == "impr" and len(p) >= 7:
            impr[tuple(_nt(x) for x in p[:4])] = (float(p[4]), int(p[5]), float(p[6]))
        elif mode == "nonb" and len(p) >= 7:
            nonb[_nt(p[0])] = s
    return {"bonds": bonds, "angles": angles, "dihedrals": dihe,
            "impropers": impr, "nonbonded": nonb}


def parse_pdb(path: Path) -> dict[str, XYZ]:
    xyz: dict[str, XYZ] = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            xyz[line[12:16].strip()] = (float(line[30:38]), float(line[38:46]),
                                        float(line[46:54]))
    return xyz


def load_ligand(stem: Path) -> Ligand:
    """Read `<stem>.{rtf,prm,pdb}` into one Ligand."""
    stem = Path(stem)
    rtf = parse_rtf(stem.with_suffix(".rtf"))
    prm = parse_prm(stem.with_suffix(".prm"))
    xyz = parse_pdb(stem.with_suffix(".pdb"))

    atoms: dict[str, Atom] = {}
    for name, typ, charge in rtf["atoms"]:
        atoms[name] = Atom(name=name, type=typ, charge=charge,
                           mass=rtf["mass"].get(typ, 0.0))

    missing = [n for n in atoms if n not in xyz]
    if missing:
        raise ChemistryError(
            f"{stem.with_suffix('.pdb')} has no coordinates for "
            f"{', '.join(missing[:8])}"
            f"{' ...' if len(missing) > 8 else ''} "
            f"({len(missing)} atom(s) declared in the .rtf but absent here).")

    return Ligand(stem=str(stem), atoms=atoms, xyz={n: xyz[n] for n in atoms},
                  bonds=tuple(rtf["bonds"]), angles=tuple(rtf["angles"]),
                  dihedrals=tuple(rtf["dihedrals"]),
                  impropers=tuple(rtf["impropers"]))


def dist(a: XYZ, b: XYZ) -> float:
    return math.dist(a, b)


def angle_deg(a: XYZ, b: XYZ, c: XYZ) -> float:
    """Angle a-b-c in degrees."""
    v = [a[i] - b[i] for i in range(3)]
    w = [c[i] - b[i] for i in range(3)]
    nv = math.sqrt(sum(x * x for x in v))
    nw = math.sqrt(sum(x * x for x in w))
    if nv == 0.0 or nw == 0.0:
        raise ChemistryError("angle_deg: two of the three points coincide")
    cos = sum(v[i] * w[i] for i in range(3)) / (nv * nw)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))
