"""
VMD integration for system building: psfgen, solvation, ionization.

Generates VMD Tcl scripts that can be run with:
    vmd -dispdev text -e <script>.tcl

Handles:
  1. Complex PSF generation (protein + hybrid ligand)
  2. Solvent PSF generation (hybrid ligand only)
  3. Solvation (water box)
  4. Ionization (neutralize or add salt)
  5. Box dimension extraction

VMD must be installed and available on PATH. The scripts are written
to disk and can be inspected before execution.
"""

from __future__ import annotations

import subprocess
import re
import tempfile
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# VMD Tcl script generators
# ---------------------------------------------------------------------------

def generate_complex_psfgen(
    *,
    protein_pdb: Path,
    hybrid_pdb: Path,
    topology_files: list[Path],
    output_dir: Path,
    output_prefix: str = "complex",
    segname: str = "LIG",
    protein_chains: Optional[list[str]] = None,
) -> Path:
    """Generate a VMD script that builds the PSF for the protein-ligand complex.

    The script handles:
    - Multiple protein chains (detected automatically if not specified)
    - Hybrid ligand as a separate segment
    - Common PDB aliases for CHARMM (HIS→HSD, ILE CD1→CD, etc.)

    Returns path to the generated .tcl script.
    """
    # Auto-detect protein chains from PDB
    if protein_chains is None:
        protein_chains = _detect_chains(protein_pdb)

    script_lines = [
        "mol delete all",
        "package require psfgen",
        "",
        "# Topology files",
    ]
    for topo in topology_files:
        script_lines.append(f"topology {topo}")
    script_lines.append(f"topology {hybrid_pdb}")

    script_lines.extend([
        "",
        "# PDB aliases for CHARMM",
        "pdbalias HIS HSD",
        "pdbalias atom SER HG HG1",
        "pdbalias residue HIS HSE",
        "pdbalias atom ILE CD1 CD",
        "",
        "# Build protein segments",
    ])

    for chain in protein_chains:
        script_lines.extend([
            f"segment {chain} {{",
            f"  pdb {protein_pdb}",
            f"  first NONE",
            f"  last NONE",
            "}",
            f"coordpdb {protein_pdb} {chain}",
        ])

    script_lines.extend([
        "",
        "# Build ligand segment",
        f"segment {segname} {{",
        f"  pdb {hybrid_pdb}",
        "  first NONE",
        "  last NONE",
        "}",
        f"coordpdb {hybrid_pdb} {segname}",
        "",
        "guesscoord",
        "regenerate angles dihedrals",
        "",
        f"writepsf {output_dir / f'{output_prefix}.psf'}",
        f"writepdb {output_dir / f'{output_prefix}.pdb'}",
        "exit",
    ])

    script_path = output_dir / f"psfgen_{output_prefix}.tcl"
    script_path.write_text("\n".join(script_lines) + "\n")
    return script_path


def generate_ligand_psfgen(
    *,
    hybrid_pdb: Path,
    topology_files: list[Path],
    output_dir: Path,
    output_prefix: str = "ligand",
) -> Path:
    """Generate a VMD script that builds the PSF for the hybrid ligand alone."""
    script_lines = [
        "mol delete all",
        "package require psfgen",
        "",
        "# Topology files",
    ]
    for topo in topology_files:
        script_lines.append(f"topology {topo}")
    script_lines.append(f"topology {hybrid_pdb}")

    script_lines.extend([
        "",
        "segment LIG {",
        f"  pdb {hybrid_pdb}",
        "  first NONE",
        "  last NONE",
        "}",
        f"coordpdb {hybrid_pdb} LIG",
        "",
        "guesscoord",
        "",
        f"writepsf {output_dir / f'{output_prefix}.psf'}",
        f"writepdb {output_dir / f'{output_prefix}.pdb'}",
        "exit",
    ])

    script_path = output_dir / f"psfgen_{output_prefix}.tcl"
    script_path.write_text("\n".join(script_lines) + "\n")
    return script_path


def generate_solvate_ionize(
    *,
    psf_file: Path,
    pdb_file: Path,
    output_dir: Path,
    output_prefix: str = "ionized",
    padding: float = 15.0,
    neutralize: bool = True,
    salt_concentration: Optional[float] = None,
    cation: str = "SOD",
    anion: str = "CLA",
) -> Path:
    """Generate VMD script for solvation and ionization.

    Args:
        psf_file: Input PSF
        pdb_file: Input PDB
        output_dir: Where to write the script and outputs
        output_prefix: Prefix for output files
        padding: Water box padding in Å (default 15)
        neutralize: If True, add counter-ions
        salt_concentration: Salt concentration in mol/L (None = no salt)
        cation: Cation type for CHARMM
        anion: Anion type for CHARMM
    """
    script_lines = [
        "mol delete all",
        f"mol load psf {psf_file} pdb {pdb_file}",
        "",
        "package require solvate",
        f"solvate {psf_file} {pdb_file} -t {padding} -o solvated",
        "",
        "package require autoionize",
    ]

    out_psf = output_dir / f"{output_prefix}.psf"
    out_pdb = output_dir / f"{output_prefix}.pdb"

    if salt_concentration is not None and salt_concentration > 0:
        script_lines.append(
            f"autoionize -psf solvated.psf -pdb solvated.pdb "
            f"-sc {salt_concentration} -cation {cation} -anion {anion} "
            f"-o {output_prefix}"
        )
    elif neutralize:
        script_lines.append(
            f"autoionize -psf solvated.psf -pdb solvated.pdb "
            f"-neutralize -o {output_prefix}"
        )
    else:
        # Just rename solvated to output
        script_lines.extend([
            "mol delete all",
            "mol load psf solvated.psf pdb solvated.pdb",
            f"set all [atomselect top all]",
            f"$all writepsf {out_psf}",
            f"$all writepdb {out_pdb}",
        ])

    script_lines.extend([
        "",
        "# Extract box dimensions and center",
        "set all [atomselect top all]",
        "measure minmax $all",
        "measure center $all",
        "",
        "exit",
    ])

    script_path = output_dir / f"solvate_{output_prefix}.tcl"
    script_path.write_text("\n".join(script_lines) + "\n")
    return script_path


# ---------------------------------------------------------------------------
# Utility: detect protein chains from PDB
# ---------------------------------------------------------------------------

def _detect_chains(pdb_path: Path) -> list[str]:
    """Extract unique chain IDs from a PDB file."""
    chains = set()
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21:22].strip()
            if chain:
                chains.add(chain)
    return sorted(chains) if chains else ["A"]


# ---------------------------------------------------------------------------
# VMD execution helpers
# ---------------------------------------------------------------------------

def run_vmd_script(script_path: Path, vmd_binary: str = "vmd") -> str:
    """Execute a VMD script and return its combined stdout+stderr.

    Raises subprocess.CalledProcessError if VMD exits non-zero.
    """
    result = subprocess.run(
        [vmd_binary, "-dispdev", "text", "-e", str(script_path)],
        capture_output=True,
        text=True,
        cwd=script_path.parent,
        timeout=600,  # 10 min timeout for solvation
    )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            [vmd_binary, "-dispdev", "text", "-e", str(script_path)],
            output=result.stdout,
            stderr=result.stderr,
        )
    return result.stdout + "\n" + result.stderr


def extract_box_info(vmd_output: str) -> dict:
    """Parse VMD 'measure minmax' and 'measure center' output.

    Returns dict with keys: box_x, box_y, box_z, box_ox, box_oy, box_oz
    """
    info = {
        "box_x": 0.0, "box_y": 0.0, "box_z": 0.0,
        "box_ox": 0.0, "box_oy": 0.0, "box_oz": 0.0,
    }

    # Find the last two measure outputs
    lines = vmd_output.splitlines()
    minmax_lines: list[str] = []
    center_lines: list[str] = []

    in_minmax = False
    in_center = False
    for line in lines:
        if "minmax" in line.lower() or (in_minmax and "{" in line):
            in_minmax = True
            continue
        if in_minmax and "}" in line:
            in_minmax = False
            continue
        if "center" in line.lower() or (in_center and "measure center" in line.lower()):
            in_center = True
            continue

    # Simpler approach: find all float triplets at end of output
    all_floats = []
    for line in lines:
        # Look for lines that are just numbers
        parts = line.strip().split()
        if len(parts) >= 3:
            try:
                vals = [float(p) for p in parts[-3:]]
                all_floats.append(vals)
            except ValueError:
                continue

    # VMD outputs min, max, then center. Last 3 lines should have these.
    if len(all_floats) >= 3:
        # min = all_floats[-3], max = all_floats[-2], center = all_floats[-1]
        box_min = all_floats[-3]
        box_max = all_floats[-2]
        center = all_floats[-1]
        info["box_x"] = abs(box_max[0] - box_min[0])
        info["box_y"] = abs(box_max[1] - box_min[1])
        info["box_z"] = abs(box_max[2] - box_min[2])
        info["box_ox"] = center[0]
        info["box_oy"] = center[1]
        info["box_oz"] = center[2]

    return info
