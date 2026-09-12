"""
Ligand parameterization using OpenMM ForceField (GAFF/OpenFF).

Replaces LigParGen web server for generating ligand force field
parameters. Two approaches are supported:

1. OpenMM + GAFF2 (via openmmforcefields) — recommended
   Uses the OpenFF Toolkit + AMBER GAFF2 for small molecule parameters.
   Compatible with CHARMM36 protein force field when used with the
   appropriate mixing rules.

2. Antechamber/ACPYPE — if ambertools is installed
   Direct GAFF parameterization via antechamber + acpype for
   conversion to CHARMM-compatible formats.

Output format: CHARMM-style RTF (topology) + PRM (parameters) files,
compatible with VMD psfgen and NAMD.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# OpenMM / OpenFF route
# ---------------------------------------------------------------------------

def parameterize_with_openmm(
    input_pdb: Path,
    output_dir: Path,
    *,
    charge: int = 0,
    forcefield: str = "gaff-2.2",
    water_model: str = "tip3p",
) -> dict:
    """Parameterize a ligand using OpenMM + GAFF.

    Uses openmmforcefields which wraps the OpenFF Toolkit to assign
    GAFF/GAFF2 atom types and generate AMBER-format parameters.

    Args:
        input_pdb: Ligand PDB file (must have CONECT records or be an SDF)
        output_dir: Where to write RTF/PRM files
        charge: Net charge of the ligand
        forcefield: GAFF version ("gaff-2.2" or "gaff-1.8")
        water_model: Water model for nonbonded parameters

    Returns:
        dict with paths to generated files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        from openmmforcefields.generators import SystemGenerator
        from openmm import unit, app
        import parmed as pmd
    except ImportError as e:
        raise RuntimeError(
            "openmmforcefields and parmed are required for ligand "
            "parameterization. Install with: pip install openmmforcefields parmed"
        ) from e

    # Load molecule
    pdb = app.PDBFile(str(input_pdb))

    # Generate system with GAFF
    barostat = None  # No barostat for parameterization
    system_generator = SystemGenerator(
        forcefields=[forcefield],
        small_molecule_forcefield=forcefield,
        molecules=[pdb.topology],
        cache=None,
    )

    # Create OpenMM System
    system = system_generator.create_system(
        pdb.topology,
        molecules=[pdb.topology],
    )

    # Convert to ParmEd for RTF/PRM output
    structure = pmd.openmm.load_topology(
        pdb.topology, system, xyz=pdb.positions
    )

    # Write CHARMM-style files
    rtf_path = output_dir / f"{input_pdb.stem}.rtf"
    prm_path = output_dir / f"{input_pdb.stem}.prm"

    structure.save(str(rtf_path), format="charmm_rtf")
    structure.save(str(prm_path), format="charmm_prm")

    return {
        "rtf": str(rtf_path),
        "prm": str(prm_path),
        "pdb": str(input_pdb),
        "charge": charge,
        "forcefield": forcefield,
    }


# ---------------------------------------------------------------------------
# Antechamber + ACPYPE route
# ---------------------------------------------------------------------------

def parameterize_with_antechamber(
    input_pdb: Path,
    output_dir: Path,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    gaff_version: str = "gaff2",
) -> dict:
    """Parameterize a ligand using antechamber + acpype.

    Requires ambertools to be installed (antechamber, sqm).
    acpype converts AMBER format to CHARMM-compatible RTF/PRM.

    Args:
        input_pdb: Ligand PDB file
        output_dir: Where to write output files
        charge: Net charge of the ligand
        multiplicity: Spin multiplicity (1=singlet, 2=doublet, etc.)
        gaff_version: "gaff" or "gaff2"

    Returns:
        dict with paths to generated files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check for antechamber
    try:
        subprocess.run(["antechamber", "-h"], capture_output=True, timeout=10)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        raise RuntimeError(
            "antechamber not found. Install ambertools: "
            "conda install -c conda-forge ambertools"
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Step 1: antechamber — assign GAFF atom types and AM1-BCC charges
        mol2_path = tmp / f"{input_pdb.stem}.mol2"
        result = subprocess.run([
            "antechamber",
            "-i", str(input_pdb),
            "-fi", "pdb",
            "-o", str(mol2_path),
            "-fo", "mol2",
            "-c", "bcc",
            "-nc", str(charge),
            "-m", str(multiplicity),
            "-at", gaff_version,
            "-pf", "y",
        ], capture_output=True, text=True, timeout=300, cwd=tmp)

        if result.returncode != 0 or not mol2_path.exists():
            raise RuntimeError(
                f"antechamber failed:\n{result.stderr[-1000:]}"
            )

        # Step 2: parmchk2 — check for missing parameters
        frcmod_path = tmp / f"{input_pdb.stem}.frcmod"
        subprocess.run([
            "parmchk2",
            "-i", str(mol2_path),
            "-f", "mol2",
            "-o", str(frcmod_path),
            "-s", gaff_version,
        ], capture_output=True, timeout=120, cwd=tmp)

        # Step 3: tleap — generate AMBER prmtop/inpcrd
        # (Skip if acpype can handle mol2 directly)
        try:
            subprocess.run(["acpype", "-h"], capture_output=True, timeout=10)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            raise RuntimeError(
                "acpype not found. Install with: pip install acpype"
            )

        result = subprocess.run([
            "acpype",
            "-i", str(mol2_path),
            "-b", input_pdb.stem,
            "-o", "charmm",
        ], capture_output=True, text=True, timeout=120, cwd=tmp)

        if result.returncode != 0:
            raise RuntimeError(
                f"acpype failed:\n{result.stderr[-1000:]}"
            )

        # Collect output files
        base = input_pdb.stem
        rtf_src = tmp / f"{base}.acpype" / f"{base}_charmm.rtf"
        prm_src = tmp / f"{base}.acpype" / f"{base}_charmm.prm"

        rtf_path = output_dir / f"{base}.rtf"
        prm_path = output_dir / f"{base}.prm"

        if rtf_src.exists():
            rtf_path.write_text(rtf_src.read_text())
        if prm_src.exists():
            prm_path.write_text(prm_src.read_text())

    return {
        "rtf": str(rtf_path),
        "prm": str(prm_path),
        "pdb": str(input_pdb),
        "charge": charge,
        "forcefield": gaff_version,
        "method": "antechamber+acpype",
    }


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def parameterize_ligand(
    input_pdb: Path,
    output_dir: Path,
    *,
    method: str = "openmm",
    charge: int = 0,
    **kwargs,
) -> dict:
    """Parameterize a ligand using the specified method.

    Args:
        input_pdb: Ligand PDB file
        output_dir: Output directory for RTF/PRM files
        method: "openmm" (default) or "antechamber"
        charge: Net charge of the ligand
        **kwargs: Forwarded to the specific parameterization function

    Returns:
        dict with 'rtf', 'prm', 'pdb' paths and metadata
    """
    if method == "openmm":
        return parameterize_with_openmm(input_pdb, output_dir, charge=charge, **kwargs)
    elif method == "antechamber":
        return parameterize_with_antechamber(input_pdb, output_dir, charge=charge, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'openmm' or 'antechamber'.")
