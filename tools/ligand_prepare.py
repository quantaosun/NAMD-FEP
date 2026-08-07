#!/usr/bin/env python3
"""Local ligand alignment and NAMD FEP PDB preparation utilities.

This tool operates on already parameterized ligand PDB files.  It does not
assign atom types or generate force-field parameters.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Atom:
    line: str
    key: str
    xyz: tuple[float, float, float]


def read_atoms(path: Path) -> list[Atom]:
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            element = line[76:78].strip().upper() or name[0].upper()
            key = f"{name.upper()}:{element}"
            atoms.append(Atom(line, key, (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    if not atoms:
        raise ValueError(f"{path} contains no ATOM/HETATM records")
    if len({atom.key for atom in atoms}) != len(atoms):
        raise ValueError(f"{path} contains duplicate atom names/elements")
    return atoms


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    return tuple(sum(point[i] for point in points) / len(points) for i in range(3))


def covariance(a, b):
    return [[sum((x[i] - a[i]) * (y[j] - b[j]) for x, y in zip(a, b)) for j in range(3)] for i in range(3)]


def rotation_from_covariance(c):
    # Horn's quaternion method, solved by power iteration on the 4x4 matrix.
    xx, xy, xz = c[0]
    yx, yy, yz = c[1]
    zx, zy, zz = c[2]
    matrix = [
        [xx + yy + zz, yz - zy, zx - xz, xy - yx],
        [yz - zy, xx - yy - zz, xy + yx, zx + xz],
        [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
        [xy - yx, zx + xz, yz + zy, -xx - yy + zz],
    ]
    q = [1.0, 0.0, 0.0, 0.0]
    for _ in range(80):
        next_q = [sum(matrix[i][j] * q[j] for j in range(4)) for i in range(4)]
        norm = math.sqrt(sum(value * value for value in next_q))
        if norm == 0:
            raise ValueError("matching atoms are collinear or coincident")
        q = [value / norm for value in next_q]
    w, x, y, z = q
    return (
        (1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)),
        (2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)),
        (2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)),
    )


def transform(point, source_center, target_center, rotation):
    shifted = tuple(point[i] - source_center[i] for i in range(3))
    return tuple(target_center[i] + sum(rotation[i][j] * shifted[j] for j in range(3)) for i in range(3))


def align_atoms(reference: list[Atom], mobile: list[Atom]) -> dict[str, tuple[float, float, float]]:
    common = sorted(set(atom.key for atom in reference) & {atom.key for atom in mobile})
    if len(common) < 3:
        raise ValueError("at least three common atoms are required for alignment")
    ref = {atom.key: atom.xyz for atom in reference}
    mob = {atom.key: atom.xyz for atom in mobile}
    ref_center = centroid([ref[key] for key in common])
    mob_center = centroid([mob[key] for key in common])
    rotation = rotation_from_covariance(covariance([mob[key] for key in common], [ref[key] for key in common]))
    return {atom.key: transform(atom.xyz, mob_center, ref_center, rotation) for atom in mobile}


def set_coordinates(line: str, xyz: tuple[float, float, float]) -> str:
    return f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def set_bfactor(line: str, value: float) -> str:
    return f"{line[:60]}{value:6.2f}{line[66:]}"


def write_aligned(input_path: Path, output_path: Path, reference_path: Path) -> None:
    mobile = read_atoms(input_path)
    reference = read_atoms(reference_path)
    transformed = align_atoms(reference, mobile)
    records = iter(transformed.values())
    output = []
    for line in input_path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = next(atom for atom in mobile if atom.line == line)
            output.append(set_coordinates(line, transformed[atom.key]))
        else:
            output.append(line)
    output_path.write_text("\n".join(output) + "\n")


def write_hybrid(a_path: Path, b_path: Path, output_path: Path, fep_path: Path) -> None:
    atoms_a, atoms_b = read_atoms(a_path), read_atoms(b_path)
    transformed_b = align_atoms(atoms_a, atoms_b)
    by_key_a = {atom.key: atom for atom in atoms_a}
    by_key_b = {atom.key: atom for atom in atoms_b}
    common = set(by_key_a) & set(by_key_b)
    lines = []
    flags = []
    for atom in atoms_a:
        lines.append(set_bfactor(atom.line, 0.0 if atom.key in common else -1.0))
        flags.append((atom.key, 0.0 if atom.key in common else -1.0))
    for atom in atoms_b:
        if atom.key not in by_key_a:
            lines.append(set_bfactor(set_coordinates(atom.line, transformed_b[atom.key]), 1.0))
            flags.append((atom.key, 1.0))
    output_path.write_text("\n".join(lines) + "\n")
    fep_path.write_text("\n".join(
        set_bfactor(line, flag) for line, (key, flag) in zip(lines, flags)
    ) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    align = sub.add_parser("align")
    align.add_argument("--reference", type=Path, required=True)
    align.add_argument("--mobile", type=Path, required=True)
    align.add_argument("--output", type=Path, required=True)
    hybrid = sub.add_parser("hybrid")
    hybrid.add_argument("--ligand-a", type=Path, required=True)
    hybrid.add_argument("--ligand-b", type=Path, required=True)
    hybrid.add_argument("--output", type=Path, required=True)
    hybrid.add_argument("--fep-output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "align":
        write_aligned(args.mobile, args.output, args.reference)
    else:
        write_hybrid(args.ligand_a, args.ligand_b, args.output, args.fep_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
