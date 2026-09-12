"""
Ligand alignment and hybrid PDB/FEP file generation.

Adapted from the independent window approach: uses atom name matching
for alignment (Kabsch via Horn's quaternion) rather than MCS, making it
dependency-free beyond the standard library.

Generates:
  - Aligned mutant ligand PDB
  - Hybrid PDB (both ligands in one file, B-factor marks λ-dependence)
  - FEP file (PDB with B-factor column: 0=common, -1=ref-only, 1=mut-only)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Atom:
    """A single ATOM/HETATM record from a PDB file."""

    line: str
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    xyz: tuple[float, float, float]
    element: str
    bfactor: float = 0.0

    @property
    def key(self) -> str:
        """Canonical key for atom matching: atom_name:element."""
        return f"{self.name.upper()}:{self.element.upper()}"


# ---------------------------------------------------------------------------
# PDB I/O
# ---------------------------------------------------------------------------

def _guess_element(name: str) -> str:
    """Guess element from PDB atom name."""
    name = name.strip()
    if not name:
        return "C"
    # Two-letter elements (first char is digit in PDB naming: 1H, 2C, etc.)
    if len(name) >= 2 and name[0].isdigit():
        elem = name[1]
    else:
        elem = name[0]
    return elem.upper()


def read_atoms(path: Path) -> list[Atom]:
    """Read all ATOM/HETATM records from a PDB file."""
    atoms: list[Atom] = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            name = line[12:16].strip()
            element = line[76:78].strip().upper() or _guess_element(name)
            atoms.append(
                Atom(
                    line=line,
                    serial=int(line[6:11]),
                    name=name,
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip(),
                    resseq=int(line[22:26]),
                    xyz=(
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ),
                    element=element,
                    bfactor=float(line[60:66]) if len(line) > 66 else 0.0,
                )
            )
        except (ValueError, IndexError):
            continue
    if not atoms:
        raise ValueError(f"{path} contains no parsable ATOM/HETATM records")
    return atoms


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    n = len(points)
    return tuple(sum(p[i] for p in points) / n for i in range(3))


def _horn_rotation(
    source: list[tuple[float, float, float]],
    target: list[tuple[float, float, float]],
    source_center: tuple[float, float, float],
    target_center: tuple[float, float, float],
) -> tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]:
    """Compute optimal rotation matrix using Horn's quaternion method.

    This is a closed-form solution to the orthogonal Procrustes problem
    (Kabsch algorithm without reflection).
    """
    # Build cross-covariance matrix
    c = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    for (sx, sy, sz), (tx, ty, tz) in zip(source, target):
        dx, dy, dz = sx - source_center[0], sy - source_center[1], sz - source_center[2]
        ex, ey, ez = tx - target_center[0], ty - target_center[1], tz - target_center[2]
        c[0][0] += dx * ex; c[0][1] += dx * ey; c[0][2] += dx * ez
        c[1][0] += dy * ex; c[1][1] += dy * ey; c[1][2] += dy * ez
        c[2][0] += dz * ex; c[2][1] += dz * ey; c[2][2] += dz * ez

    # Build 4x4 symmetric matrix for quaternion eigenvalue problem
    xx, xy, xz = c[0]
    yx, yy, yz = c[1]
    zx, zy, zz = c[2]
    m = [
        [xx + yy + zz,  yz - zy,      zx - xz,      xy - yx],
        [yz - zy,       xx - yy - zz,  xy + yx,      zx + xz],
        [zx - xz,       xy + yx,      -xx + yy - zz, yz + zy],
        [xy - yx,       zx + xz,       yz + zy,     -xx - yy + zz],
    ]

    # Power iteration for dominant eigenvector
    q = [1.0, 0.0, 0.0, 0.0]
    for _ in range(100):
        nq = [sum(m[i][j] * q[j] for j in range(4)) for i in range(4)]
        norm = math.sqrt(sum(v * v for v in nq))
        if norm < 1e-15:
            raise ValueError("degenerate point set — atoms may be collinear")
        nq = [v / norm for v in nq]
        if sum(abs(nq[i] - q[i]) for i in range(4)) < 1e-14:
            break
        q = nq

    w, x, y, z = q
    return (
        (w*w + x*x - y*y - z*z, 2*(x*y - w*z),       2*(x*z + w*y)),
        (2*(x*y + w*z),        w*w - x*x + y*y - z*z, 2*(y*z - w*x)),
        (2*(x*z - w*y),        2*(y*z + w*x),        w*w - x*x - y*y + z*z),
    )


def apply_transform(
    point: tuple[float, float, float],
    source_center: tuple[float, float, float],
    target_center: tuple[float, float, float],
    rotation: tuple,
) -> tuple[float, float, float]:
    """Apply rotation + translation to a point."""
    shifted = tuple(point[i] - source_center[i] for i in range(3))
    return tuple(
        target_center[i] + sum(rotation[i][j] * shifted[j] for j in range(3))
        for i in range(3)
    )


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def align_ligands(
    reference: list[Atom],
    mobile: list[Atom],
) -> dict[str, tuple[float, float, float]]:
    """Compute optimal rigid-body transform from mobile to reference.

    Uses common atoms (matched by atom_key = name:element) as anchors.
    Requires at least 3 common atoms.
    """
    ref_keys = {a.key for a in reference}
    mob_keys = {a.key for a in mobile}
    common = sorted(ref_keys & mob_keys)

    if len(common) < 3:
        raise ValueError(
            f"Need at least 3 common atoms for alignment, found {len(common)}. "
            f"Common: {common}"
        )

    ref_map = {a.key: a.xyz for a in reference}
    mob_map = {a.key: a.xyz for a in mobile}

    ref_pts = [ref_map[k] for k in common]
    mob_pts = [mob_map[k] for k in common]

    ref_ctr = centroid(ref_pts)
    mob_ctr = centroid(mob_pts)

    rotation = _horn_rotation(mob_pts, ref_pts, mob_ctr, ref_ctr)

    return {
        a.key: apply_transform(a.xyz, mob_ctr, ref_ctr, rotation)
        for a in mobile
    }


# ---------------------------------------------------------------------------
# PDB coordinate manipulation
# ---------------------------------------------------------------------------

def _set_xyz(line: str, xyz: tuple[float, float, float]) -> str:
    """Replace coordinates in a PDB ATOM line."""
    return f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def _set_bfactor(line: str, value: float) -> str:
    """Replace B-factor in a PDB ATOM line."""
    return f"{line[:60]}{value:6.2f}{line[66:]}"


def write_aligned_pdb(
    input_path: Path,
    output_path: Path,
    reference_path: Path,
) -> None:
    """Align mobile ligand (input) to reference ligand, write result."""
    mobile = read_atoms(input_path)
    reference = read_atoms(reference_path)
    transforms = align_ligands(reference, mobile)

    lines_out: list[str] = []
    for line in input_path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            element = _guess_element(name)
            key = f"{name.upper()}:{element.upper()}"
            if key in transforms:
                lines_out.append(_set_xyz(line, transforms[key]))
            else:
                lines_out.append(line)
        else:
            lines_out.append(line)

    output_path.write_text("\n".join(lines_out) + "\n")


# ---------------------------------------------------------------------------
# Hybrid PDB & FEP file generation (dual topology at coordinate level)
# ---------------------------------------------------------------------------

def build_hybrid(
    ref_atoms: list[Atom],
    mut_atoms: list[Atom],
    transforms: dict[str, tuple[float, float, float]],
) -> tuple[list[str], list[tuple[str, float]]]:
    """Build a hybrid PDB with both ligands' atoms.

    Common atoms (matched by key) appear once with B-factor 0.
    Reference-unique atoms appear with B-factor -1.
    Mutant-unique atoms appear with B-factor +1.

    Returns (pdb_lines, [(atom_key, bfactor), ...]) for FEP file generation.
    """
    ref_keys = {a.key for a in ref_atoms}
    mut_keys = {a.key for a in mut_atoms}
    common = ref_keys & mut_keys

    hybrid_lines: list[str] = []
    fep_marks: list[tuple[str, float]] = []

    # Reference atoms first
    for a in ref_atoms:
        bfac = 0.0 if a.key in common else -1.0
        hybrid_lines.append(_set_bfactor(a.line, bfac))
        fep_marks.append((a.key, bfac))

    # Mutant-unique atoms (already transformed to reference frame)
    for a in mut_atoms:
        if a.key not in ref_keys:
            xformed_line = _set_xyz(a.line, transforms[a.key])
            hybrid_lines.append(_set_bfactor(xformed_line, 1.0))
            fep_marks.append((a.key, 1.0))

    return hybrid_lines, fep_marks


def write_hybrid_and_fep(
    ref_path: Path,
    mut_path: Path,
    hybrid_out: Path,
    fep_out: Path,
) -> tuple[int, int, int]:
    """Generate hybrid PDB and FEP file from two aligned ligand PDBs.

    Returns (n_common, n_ref_only, n_mut_only) atom counts.
    """
    ref_atoms = read_atoms(ref_path)
    mut_atoms = read_atoms(mut_path)
    transforms = align_ligands(ref_atoms, mut_atoms)

    hybrid_lines, fep_marks = build_hybrid(ref_atoms, mut_atoms, transforms)

    # Write hybrid PDB
    hybrid_out.write_text("\n".join(hybrid_lines) + "\n")

    # Write FEP file (same as hybrid but all atoms marked)
    # For the complex/solvent .fep files, we need to mark ALL atoms in the
    # full system. This is the ligand-only FEP reference.
    fep_lines: list[str] = []
    for (key, bfac), line in zip(fep_marks, hybrid_lines):
        fep_lines.append(_set_bfactor(line, bfac))
    fep_out.write_text("\n".join(fep_lines) + "\n")

    # Counts
    ref_keys = {a.key for a in ref_atoms}
    mut_keys = {a.key for a in mut_atoms}
    common = ref_keys & mut_keys
    n_common = len(common)
    n_ref_only = len(ref_keys - mut_keys)
    n_mut_only = len(mut_keys - ref_keys)

    return n_common, n_ref_only, n_mut_only
