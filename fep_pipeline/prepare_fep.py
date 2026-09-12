#!/usr/bin/env python3
"""
NAMD-FEP Preparation Pipeline — Main Orchestrator.

Replaces the Feprepare web server + Maestro + LigParGen workflow with a
single, self-contained Python 3 tool. Uses VMD for psfgen/solvate/ionize
and produces ready-to-run NAMD configuration files.

Usage:
    prepare_fep.py run \\
        --protein protein.pdb \\
        --ligand-a ref.pdb --ligand-b mut.pdb \\
        --workdir ./fep_run

The workflow:
  1. Align ligands via atom-name matching (Kabsch algorithm)
  2. Build hybrid PDB (dual topology at coordinate level)
  3. Generate VMD scripts for psfgen, solvation, ionization
  4. Run VMD to build systems (complex + solvent)
  5. Generate FEP files (B-factor column for λ-dependence)
  6. Write NAMD config files for equilibration + production
  7. (Optional) Run the simulation

Requirements:
  - Python 3.8+ with numpy
  - VMD (for system building)
  - NAMD 2.13+ or NAMD3 (pre-compiled, for simulation)
  - Force field files: CHARMM36 toppar directory
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

from . import __version__
from .alignment import (
    read_atoms,
    align_ligands,
    write_aligned_pdb,
    write_hybrid_and_fep,
)
from .namd_config import (
    NAMDConfigGenerator,
    write_fep_tcl,
    namd_gpu_command,
    namd_cpu_command,
    lambda_windows,
    FEP_TCL_CONTENT,
)
from .system_builder import (
    generate_complex_psfgen,
    generate_ligand_psfgen,
    generate_solvate_ionize,
    run_vmd_script,
    extract_box_info,
)


# ---------------------------------------------------------------------------
# Default force field parameters
# ---------------------------------------------------------------------------

# Standard CHARMM36 parameter files used by CHARMM-GUI for protein+ligand FEP
DEFAULT_PARAMETER_FILES = [
    "par_all36m_prot.prm",
    "par_all36_na.prm",
    "par_all36_carb.prm",
    "par_all36_lipid.prm",
    "par_all36_cgenff.prm",
    "toppar_water_ions.str",
]


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

class FEPPipeline:
    """Orchestrate the full FEP preparation workflow."""

    def __init__(
        self,
        workdir: Path,
        *,
        toppar_dir: Optional[Path] = None,
        vmd_binary: str = "vmd",
        temperature: float = 300.0,
        n_windows: int = 16,
        padding: float = 15.0,
        salt_concentration: Optional[float] = None,
    ):
        self.workdir = Path(workdir)
        self.toppar_dir = Path(toppar_dir) if toppar_dir else self.workdir / "toppar"
        self.vmd_binary = vmd_binary
        self.temperature = temperature
        self.n_windows = n_windows
        self.padding = padding
        self.salt_concentration = salt_concentration

        # Directory structure
        self.complex_dir = self.workdir / "complex"
        self.solvent_dir = self.workdir / "solvent"
        self.input_dir = self.workdir / "input_files"

    # ------------------------------------------------------------------
    # Step 1: Ligand alignment and hybrid PDB
    # ------------------------------------------------------------------

    def prepare_ligands(
        self,
        ref_pdb: Path,
        mut_pdb: Path,
    ) -> dict:
        """Align ligands and build the hybrid dual-topology PDB."""
        self.input_dir.mkdir(parents=True, exist_ok=True)

        ref_atoms = read_atoms(ref_pdb)
        mut_atoms = read_atoms(mut_pdb)
        transforms = align_ligands(ref_atoms, mut_atoms)

        # Write aligned mutant
        aligned_mut = self.input_dir / "mutant_aligned.pdb"
        write_aligned_pdb(mut_pdb, aligned_mut, ref_pdb)

        # Build hybrid PDB
        hybrid_path = self.input_dir / "hybrid_ligand.pdb"
        fep_ref_path = self.input_dir / "ligand_fep_ref.pdb"

        n_common, n_ref_only, n_mut_only = write_hybrid_and_fep(
            ref_pdb, mut_pdb, hybrid_path, fep_ref_path
        )

        info = {
            "ref_pdb": str(ref_pdb),
            "mut_pdb": str(mut_pdb),
            "aligned_mut": str(aligned_mut),
            "hybrid_pdb": str(hybrid_path),
            "fep_ref": str(fep_ref_path),
            "n_atoms_ref": len(ref_atoms),
            "n_atoms_mut": len(mut_atoms),
            "n_common": n_common,
            "n_ref_only": n_ref_only,
            "n_mut_only": n_mut_only,
            "n_hybrid_total": len(ref_atoms) + n_mut_only,
        }
        print(f"  Ligand alignment: {n_common} common, {n_ref_only} ref-only, "
              f"{n_mut_only} mut-only atoms")
        return info

    # ------------------------------------------------------------------
    # Step 2: System building (VMD psfgen + solvate + ionize)
    # ------------------------------------------------------------------

    def build_systems(
        self,
        protein_pdb: Path,
        hybrid_pdb: Path,
        *,
        parameter_files: Optional[list[str]] = None,
        topology_files: Optional[list[Path]] = None,
    ) -> dict:
        """Build both complex and solvent systems using VMD.

        Returns dict with psf/pdb paths and box info for each leg.
        """
        if parameter_files is None:
            parameter_files = DEFAULT_PARAMETER_FILES

        self.complex_dir.mkdir(parents=True, exist_ok=True)
        self.solvent_dir.mkdir(parents=True, exist_ok=True)

        # Resolve topology files
        if topology_files is None:
            topology_files = []
            for pf in parameter_files:
                candidate = self.toppar_dir / pf
                if candidate.exists():
                    topology_files.append(candidate)

        # Write fep.tcl into both directories
        write_fep_tcl(self.complex_dir)
        write_fep_tcl(self.solvent_dir)

        box_info = {}

        # --- Complex leg ---
        print("\n  Building COMPLEX system...")
        box_info["complex"] = self._build_one_system(
            output_dir=self.complex_dir,
            protein_pdb=protein_pdb,
            hybrid_pdb=hybrid_pdb,
            topology_files=topology_files,
            parameter_files=parameter_files,
            output_prefix="ionized",
            is_complex=True,
        )

        # --- Solvent leg ---
        print("\n  Building SOLVENT system...")
        box_info["solvent"] = self._build_one_system(
            output_dir=self.solvent_dir,
            protein_pdb=None,
            hybrid_pdb=hybrid_pdb,
            topology_files=topology_files,
            parameter_files=parameter_files,
            output_prefix="ionized",
            is_complex=False,
        )

        return box_info

    def _build_one_system(
        self,
        output_dir: Path,
        protein_pdb: Optional[Path],
        hybrid_pdb: Path,
        topology_files: list[Path],
        parameter_files: list[str],
        output_prefix: str,
        is_complex: bool,
    ) -> dict:
        """Build one system (complex or solvent)."""
        leg = "complex" if is_complex else "solvent"

        # --- PSF generation ---
        if is_complex:
            psfgen_script = generate_complex_psfgen(
                protein_pdb=protein_pdb,
                hybrid_pdb=hybrid_pdb,
                topology_files=topology_files,
                output_dir=output_dir,
                output_prefix=leg,
            )
        else:
            psfgen_script = generate_ligand_psfgen(
                hybrid_pdb=hybrid_pdb,
                topology_files=topology_files,
                output_dir=output_dir,
                output_prefix=leg,
            )

        print(f"    PSF script: {psfgen_script.name}")
        vmd_out = run_vmd_script(psfgen_script, vmd_binary=self.vmd_binary)

        psf_file = output_dir / f"{leg}.psf"
        pdb_file = output_dir / f"{leg}.pdb"

        if not psf_file.exists():
            raise RuntimeError(
                f"VMD psfgen failed for {leg} leg.\n"
                f"VMD output:\n{vmd_out[-2000:]}"
            )

        # --- Solvate + ionize ---
        solvate_script = generate_solvate_ionize(
            psf_file=psf_file,
            pdb_file=pdb_file,
            output_dir=output_dir,
            output_prefix=output_prefix,
            padding=self.padding,
            salt_concentration=self.salt_concentration,
        )

        print(f"    Solvate script: {solvate_script.name}")
        vmd_out = run_vmd_script(solvate_script, vmd_binary=self.vmd_binary)

        ionized_psf = output_dir / f"{output_prefix}.psf"
        ionized_pdb = output_dir / f"{output_prefix}.pdb"

        if not ionized_psf.exists():
            raise RuntimeError(
                f"VMD solvate/ionize failed for {leg} leg.\n"
                f"VMD output:\n{vmd_out[-2000:]}"
            )

        box_info = extract_box_info(vmd_out)
        print(f"    Box: {box_info['box_x']:.1f} × {box_info['box_y']:.1f} "
              f"× {box_info['box_z']:.1f} Å³")

        return {
            "psf": str(ionized_psf),
            "pdb": str(ionized_pdb),
            "box_info": box_info,
            "vmd_log": vmd_out,
        }

    # ------------------------------------------------------------------
    # Step 3: Generate NAMD config files
    # ------------------------------------------------------------------

    def write_namd_configs(
        self,
        box_info: dict,
        *,
        parameter_files: Optional[list[str]] = None,
        nvt_steps: int = 50000,
        nvt_min_steps: int = 5000,
        npt_steps: int = 50000,
        prod_steps: int = 500000,
        output_freq: int = 500,
        fep_out_freq: int = 500,
        alch_equil_steps: int = 50000,
    ) -> dict:
        """Generate all NAMD configuration files."""
        if parameter_files is None:
            parameter_files = DEFAULT_PARAMETER_FILES

        # Prefix parameter files with relative path
        param_paths = [f"../../toppar/{pf}" for pf in parameter_files]

        gen = NAMDConfigGenerator(
            temperature=self.temperature,
            parameter_files=param_paths,
            output_freq=output_freq,
            fep_out_freq=fep_out_freq,
            alch_equil_steps=alch_equil_steps,
            nvt_steps=nvt_steps,
            nvt_min_steps=nvt_min_steps,
            npt_steps=npt_steps,
            prod_steps=prod_steps,
            prod_steps_per_window=prod_steps,
        )

        configs = {}

        for leg in ["complex", "solvent"]:
            leg_dir = self.workdir / leg
            leg_dir.mkdir(parents=True, exist_ok=True)
            leg_info = box_info[leg]
            b = leg_info["box_info"]
            psf = f"ionized.psf" if leg == "complex" else f"ionized.psf"
            pdb = f"{leg}_ionized.pdb" if leg == "complex" else "ionized.pdb"
            fep_file = f"ionized_{leg}.fep"

            # Determine correct PDB/PSF filenames
            _psf = "ionized.psf"
            _pdb = "ionized.pdb"

            # NVT
            nvt_path = gen.write_nvt(
                leg_dir, leg,
                psf_file=_psf, pdb_file=_pdb, fep_file=fep_file,
                output_name="nvt_equil",
                box_info=b,
            )
            configs[f"{leg}_nvt"] = str(nvt_path)

            # NPT
            npt_path = gen.write_npt(
                leg_dir, leg,
                psf_file=_psf, pdb_file=_pdb, fep_file=fep_file,
                output_name="npt_equil",
                prev_coor="nvt_equil.coor",
                prev_vel="nvt_equil.vel",
                prev_xsc="nvt_equil.xsc",
            )
            configs[f"{leg}_npt"] = str(npt_path)

            # Forward production
            fwd_path = gen.write_production(
                leg_dir, leg, "forward",
                psf_file=_psf, pdb_file=_pdb, fep_file=fep_file,
                output_name="md_forward",
                prev_coor="npt_equil.coor",
                prev_vel="npt_equil.vel",
                prev_xsc="npt_equil.xsc",
                lambda_start=0.0, lambda_end=1.0,
                n_windows=self.n_windows,
            )
            configs[f"{leg}_forward"] = str(fwd_path)

            # Backward production
            bwd_path = gen.write_production(
                leg_dir, leg, "backward",
                psf_file=_psf, pdb_file=_pdb, fep_file=fep_file,
                output_name="md_backward",
                prev_coor="md_forward.coor",
                prev_vel="md_forward.vel",
                prev_xsc="md_forward.xsc",
                lambda_start=1.0, lambda_end=0.0,
                n_windows=self.n_windows,
            )
            configs[f"{leg}_backward"] = str(bwd_path)

        return configs

    # ------------------------------------------------------------------
    # Step 4: Generate run script
    # ------------------------------------------------------------------

    def write_run_script(
        self,
        *,
        namd_binary: str = "namd3",
        cpus: int = 1,
        gpu_device: int = 0,
        use_gpu: bool = True,
    ) -> Path:
        """Write a convenience shell script to run all simulations.

        ``cpus`` defaults to 1: in GPU-resident mode throughput scales inversely
        with PE count (see ``namd_config.namd_gpu_command``).
        """
        lines = [
            "#!/bin/bash",
            "# Auto-generated FEP run script",
            "# Generated by fep_pipeline v{}".format(__version__),
            "",
            "set -euo pipefail",
            "",
            f"NAMD={namd_binary}",
            f"CPUS={cpus}",
            f"DEVICE={gpu_device}",
            "",
        ]

        if use_gpu:
            lines.append('# NOTE: CPUS=1 is intentional -- GPU-resident throughput scales')
            lines.append('#       inversely with PE count (+p1 75.8 vs +p8 13.2 ns/day on a V100).')
            lines.append('NAMD_OPTS="+p${CPUS} +setcpuaffinity --GPUresident on +devices ${DEVICE}"')
        else:
            lines.append('NAMD_OPTS="+p${CPUS}"')

        lines.extend([
            "",
            'echo "=========================================="',
            'echo "NAMD-FEP Simulation"',
            'echo "=========================================="',
            "",
        ])

        for leg in ["complex", "solvent"]:
            lines.extend([
                f'echo ""',
                f'echo "=== {leg.upper()} LEG ==="',
                f'echo "--- NVT ---"',
                f'cd {leg}',
                f'$NAMD $NAMD_OPTS nvt_equil.namd | tee nvt_equil.log',
                f'echo "--- NPT ---"',
                f'$NAMD $NAMD_OPTS npt_equil.namd | tee npt_equil.log',
                f'echo "--- Forward ---"',
                f'$NAMD $NAMD_OPTS md_forward.namd | tee md_forward.log',
                f'echo "--- Backward ---"',
                f'$NAMD $NAMD_OPTS md_backward.namd | tee md_backward.log',
                f'cd ..',
                f'echo "=== {leg.upper()} COMPLETE ==="',
                "",
            ])

        lines.extend([
            'echo ""',
            'echo "=========================================="',
            'echo "All simulations complete."',
            'echo "Run ParseFEP in VMD to analyze results."',
            'echo "=========================================="',
        ])

        script_path = self.workdir / "run_fep.sh"
        script_path.write_text("\n".join(lines) + "\n")
        script_path.chmod(0o755)
        return script_path

    # ------------------------------------------------------------------
    # Full pipeline
    # ------------------------------------------------------------------

    def run(
        self,
        protein_pdb: Path,
        ref_pdb: Path,
        mut_pdb: Path,
        *,
        skip_vmd: bool = False,
    ) -> Path:
        """Run the complete pipeline and return path to the status file."""
        self.workdir.mkdir(parents=True, exist_ok=True)

        status = {
            "version": __version__,
            "workdir": str(self.workdir),
            "steps": {},
        }

        print("=" * 60)
        print("NAMD-FEP Preparation Pipeline")
        print("=" * 60)

        # Step 1: Ligand prep
        print("\n[1/4] Preparing ligands (alignment + hybrid PDB)...")
        lig_info = self.prepare_ligands(ref_pdb, mut_pdb)
        status["steps"]["ligands"] = lig_info

        # Copy hybrid PDB to input_files for reference
        hybrid_pdb = Path(lig_info["hybrid_pdb"])

        # Step 2: System building
        print("\n[2/4] Building systems (VMD psfgen + solvate + ionize)...")
        if skip_vmd:
            print("  SKIPPED (--skip-vmd). Using placeholder box info.")
            # Provide default box info for both legs. User must update
            # cellBasisVector/cellOrigin in the generated .namd files
            # with values from vmd_log.txt after running VMD manually.
            default_box = {
                "psf": "ionized.psf", "pdb": "ionized.pdb",
                "box_info": {"box_x": 0.0, "box_y": 0.0, "box_z": 0.0,
                             "box_ox": 0.0, "box_oy": 0.0, "box_oz": 0.0},
                "vmd_log": "",
            }
            box_info = {"complex": dict(default_box), "solvent": dict(default_box)}
        else:
            box_info = self.build_systems(
                protein_pdb=protein_pdb,
                hybrid_pdb=hybrid_pdb,
            )
        status["steps"]["systems"] = box_info

        # Step 3: NAMD configs
        print("\n[3/4] Writing NAMD configuration files...")
        configs = self.write_namd_configs(box_info)
        status["steps"]["namd_configs"] = configs

        # Step 4: Run script
        print("\n[4/4] Writing run script...")
        run_script = self.write_run_script()
        status["steps"]["run_script"] = str(run_script)

        # Write status
        status_path = self.workdir / "fep_status.json"
        status_path.write_text(json.dumps(status, indent=2, default=str))

        print(f"\n{'=' * 60}")
        print(f"Preparation complete!")
        print(f"  Work directory: {self.workdir}")
        print(f"  Run script:     {run_script}")
        print(f"  Status file:    {status_path}")
        print(f"\nNext steps:")
        print(f"  1. Check the generated files in {self.workdir}")
        print(f"  2. Run:  bash {run_script}")
        print(f"  3. Analyze with VMD ParseFEP")
        print(f"{'=' * 60}")

        return status_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="NAMD-FEP Preparation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full preparation:
  prepare_fep.py run --protein 4w53_prepared.pdb \\
      --ligand-a XLC.pdb --ligand-b XLD.pdb \\
      --workdir ./my_fep_run --toppar ./toppar

  # Preparation only (skip VMD if systems already built):
  prepare_fep.py run --protein prot.pdb \\
      --ligand-a ref.pdb --ligand-b mut.pdb \\
      --workdir ./fep --skip-vmd

  # Align ligands only:
  prepare_fep.py align --reference ref.pdb --mobile mut.pdb \\
      --output mut_aligned.pdb

  # Hybrid PDB generation only:
  prepare_fep.py hybrid --ligand-a ref.pdb --ligand-b mut.pdb \\
      --output hybrid.pdb --fep-output ligand.fep
        """,
    )

    parser.add_argument("--version", action="version",
                        version=f"fep_pipeline {__version__}")
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # ---- run ----
    run_parser = subparsers.add_parser("run", help="Run full preparation pipeline")
    run_parser.add_argument("--protein", type=Path, required=True,
                            help="Prepared protein PDB (no hydrogens needed)")
    run_parser.add_argument("--ligand-a", type=Path, required=True,
                            help="Reference ligand PDB (with force field params)")
    run_parser.add_argument("--ligand-b", type=Path, required=True,
                            help="Mutant ligand PDB (with force field params)")
    run_parser.add_argument("--workdir", type=Path, default=Path("./fep_run"),
                            help="Working directory (default: ./fep_run)")
    run_parser.add_argument("--toppar", type=Path,
                            help="Path to CHARMM36 toppar directory")
    run_parser.add_argument("--temperature", type=float, default=300.0)
    run_parser.add_argument("--n-windows", type=int, default=16)
    run_parser.add_argument("--padding", type=float, default=15.0)
    run_parser.add_argument("--salt", type=float, default=None,
                            help="Salt concentration in mol/L")
    run_parser.add_argument("--vmd-binary", default="vmd")
    run_parser.add_argument("--namd-binary", default="namd3")
    run_parser.add_argument("--cpus", type=int, default=4)
    run_parser.add_argument("--gpu-device", type=int, default=0)
    run_parser.add_argument("--no-gpu", action="store_true")
    run_parser.add_argument("--skip-vmd", action="store_true",
                            help="Skip VMD system building (use existing files)")
    run_parser.add_argument("--nvt-steps", type=int, default=50000)
    run_parser.add_argument("--npt-steps", type=int, default=50000)
    run_parser.add_argument("--prod-steps", type=int, default=500000)

    # ---- align ----
    align_parser = subparsers.add_parser("align", help="Align two ligand PDBs")
    align_parser.add_argument("--reference", type=Path, required=True)
    align_parser.add_argument("--mobile", type=Path, required=True)
    align_parser.add_argument("--output", type=Path, required=True)

    # ---- hybrid ----
    hybrid_parser = subparsers.add_parser("hybrid",
                                          help="Generate hybrid PDB + FEP file")
    hybrid_parser.add_argument("--ligand-a", type=Path, required=True)
    hybrid_parser.add_argument("--ligand-b", type=Path, required=True)
    hybrid_parser.add_argument("--output", type=Path, required=True)
    hybrid_parser.add_argument("--fep-output", type=Path, required=True)

    # ---- windows ----
    win_parser = subparsers.add_parser("windows", help="Print lambda schedule")
    win_parser.add_argument("--n-windows", type=int, default=16)

    # ---- command ----
    cmd_parser = subparsers.add_parser("command",
                                       help="Print NAMD command line for a config")
    cmd_parser.add_argument("--config", type=Path, required=True)
    cmd_parser.add_argument("--namd", default="namd3")
    cmd_parser.add_argument("--cpus", type=int, default=4)
    cmd_parser.add_argument("--device", type=int, default=0)
    cmd_parser.add_argument("--no-gpu", action="store_true")

    args = parser.parse_args(argv)

    if args.command == "run":
        pipeline = FEPPipeline(
            workdir=args.workdir,
            toppar_dir=args.toppar,
            vmd_binary=args.vmd_binary,
            temperature=args.temperature,
            n_windows=args.n_windows,
            padding=args.padding,
            salt_concentration=args.salt,
        )
        pipeline.run(
            protein_pdb=args.protein,
            ref_pdb=args.ligand_a,
            mut_pdb=args.ligand_b,
            skip_vmd=args.skip_vmd,
        )
        pipeline.write_run_script(
            namd_binary=args.namd_binary,
            cpus=args.cpus,
            gpu_device=args.gpu_device,
            use_gpu=not args.no_gpu,
        )

    elif args.command == "align":
        from .alignment import write_aligned_pdb
        write_aligned_pdb(args.mobile, args.output, args.reference)
        print(f"Aligned PDB written to {args.output}")

    elif args.command == "hybrid":
        from .alignment import write_hybrid_and_fep
        n_c, n_r, n_m = write_hybrid_and_fep(
            args.ligand_a, args.ligand_b, args.output, args.fep_output
        )
        print(f"Hybrid: {n_c} common, {n_r} ref-only, {n_m} mut-only atoms")
        print(f"Hybrid PDB: {args.output}")
        print(f"FEP ref:    {args.fep_output}")

    elif args.command == "windows":
        for i, (lam, lam2) in enumerate(lambda_windows(args.n_windows)):
            print(f"Window {i:02d}: λ={lam:.6f} → λ={lam2:.6f}")

    elif args.command == "command":
        if args.no_gpu:
            print(namd_cpu_command(args.namd, str(args.config), args.cpus))
        else:
            print(namd_gpu_command(args.namd, str(args.config),
                                   args.cpus, args.device))

    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
