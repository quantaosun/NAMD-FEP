"""
FEP file generation for NAMD alchemical simulations.

The .fep file is a PDB file where the B-factor (occupancy) column marks
which atoms are λ-dependent:

  B-factor =  0  →  common atoms (always present, always interacting)
  B-factor = -1  →  unique to reference ligand (disappears as λ→1)
  B-factor = +1  →  unique to mutant ligand (appears as λ→1)

This module takes the hybrid ligand atom map and applies the correct
B-factor markings to atoms in the full solvated system PDB.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .alignment import read_atoms, _set_bfactor


def _atom_identity_key(line: str) -> tuple[str, str, str]:
    """Extract (name, resname, chain) from a PDB ATOM/HETATM line.

    These three fields together uniquely identify an atom in the context
    of a protein-ligand complex. For matching against the hybrid ligand
    reference, we use atom name + element.
    """
    name = line[12:16].strip()
    resname = line[17:20].strip()
    chain = line[21:22].strip()
    return (name, resname, chain)


def build_system_fep(
    system_pdb: Path,
    hybrid_fep_ref: Path,
    output_path: Path,
    *,
    ligand_resname: str = "LIG",
    ligand_chain: str = "X",
) -> Path:
    """Generate the .fep file for a full solvated system.

    Reads the full system PDB (protein + hybrid ligand + water + ions)
    and marks atoms based on the hybrid ligand FEP reference.

    Ligand atoms get their B-factor from the reference (-1, 0, or +1).
    All other atoms (protein, water, ions) get B-factor = 0 (always present).

    Args:
        system_pdb: Full solvated/ionized system PDB
        hybrid_fep_ref: Ligand-only FEP reference PDB (from build_hybrid)
        output_path: Where to write the .fep file
        ligand_resname: Residue name of the hybrid ligand in the system
        ligand_chain: Chain ID of the hybrid ligand

    Returns:
        Path to the generated .fep file
    """
    # Read the ligand FEP reference to build a mapping of atom_name → bfactor
    lig_bfactor_map: dict[str, float] = {}
    for line in hybrid_fep_ref.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            element = line[76:78].strip().upper()
            if not element:
                # Guess element
                n = name.strip()
                element = (n[1] if len(n) >= 2 and n[0].isdigit() else n[0]).upper()
            key = f"{name.upper()}:{element.upper()}"
            bfac = float(line[60:66]) if len(line) > 66 else 0.0
            lig_bfactor_map[key] = bfac

    # Process the full system PDB
    out_lines: list[str] = []
    for line in system_pdb.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            resname = line[17:20].strip()

            if resname == ligand_resname:
                # This is a ligand atom — look up its B-factor
                name = line[12:16].strip()
                element = line[76:78].strip().upper()
                if not element:
                    n = name.strip()
                    element = (n[1] if len(n) >= 2 and n[0].isdigit() else n[0]).upper()
                key = f"{name.upper()}:{element.upper()}"
                bfac = lig_bfactor_map.get(key, 0.0)
                out_lines.append(_set_bfactor(line, bfac))
            else:
                # Protein, water, ion — always present
                out_lines.append(_set_bfactor(line, 0.0))
        else:
            out_lines.append(line)

    output_path.write_text("\n".join(out_lines) + "\n")
    return output_path


def validate_fep_file(fep_path: Path) -> dict:
    """Check that a .fep file has the expected structure.

    Returns summary dict with counts of -1, 0, +1 B-factor atoms.
    """
    counts = {"neg_one": 0, "zero": 0, "pos_one": 0, "other": 0}
    for line in fep_path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            bfac = float(line[60:66]) if len(line) > 66 else 0.0
            if bfac < -0.5:
                counts["neg_one"] += 1
            elif bfac > 0.5:
                counts["pos_one"] += 1
            elif abs(bfac) < 0.01:
                counts["zero"] += 1
            else:
                counts["other"] += 1
    return counts
