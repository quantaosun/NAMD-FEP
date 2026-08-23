#!/usr/bin/env python3
"""
Build the complex and solvent systems for the 6I5I dual-topology FEP.

Pure VMD + NAMD workflow (no FEPrepare):
  complex  : protein (2 segments) + hybrid ligand -> solvate -> ionize
  solvent  : hybrid ligand alone -> solvate -> ionize

Then marks the solvated system's B-factor column (alchFile .fep) so NAMD
knows which atoms vanish/appear, and writes fep.tcl + the NAMD configs.

Outputs (./complex and ./solvent): psf/pdb, *.fep, *.namd; plus ../fep.tcl.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
HYB = HERE / "hybrid"
TOPPAR = HERE.parent / "toppar"
VMD = "/home/aistudio/vmd-env/bin/vmd"
NAMD = "/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++/namd3"
PROT = HERE / "inputs" / "protein.pdb"

# protein segments (chain A has a gap between 412 and 416)
SEGMENTS = [("P1", 148, 412), ("P2", 416, 483)]


def run_vmd(script: Path, cwd: Path) -> str:
    r = subprocess.run(
        [VMD, "-dispdev", "text", "-e", str(script)],
        capture_output=True, text=True, cwd=str(cwd), timeout=900,
    )
    return r.stdout + "\n" + r.stderr


def split_protein(outdir: Path) -> list[Path]:
    """Write one PDB per segment (residues in the given numbering range)."""
    paths = []
    for seg, lo, hi in SEGMENTS:
        out = []
        for line in PROT.read_text().splitlines():
            if line.startswith("ATOM"):
                rn = int(line[22:26].strip())
                if lo <= rn <= hi:
                    out.append(line)
        out.append("END")
        p = outdir / f"protein_{seg}.pdb"
        p.write_text("\n".join(out) + "\n")
        paths.append(p)
    return paths


def psfgen_complex(outdir: Path, seg_paths: list[Path]) -> Path:
    lines = [
        "package require psfgen",
        f"topology {TOPPAR / 'top_all36_prot.rtf'}",
        f"topology {HYB / 'hybrid.rtf'}",
        "pdbalias atom ILE CD1 CD",
        "pdbalias atom SER HG HG1",
        "pdbalias residue HIS HSD",
    ]
    for (seg, lo, hi), p in zip(SEGMENTS, seg_paths):
        lines += [
            f"segment {seg} {{ first NTER; last CTER; pdb {p} }}",
            f"coordpdb {p} {seg}",
        ]
    lines += [
        f"segment LIG {{ first none; last none; pdb {HYB / 'hybrid.pdb'} }}",
        f"coordpdb {HYB / 'hybrid.pdb'} LIG",
        "guesscoord",
        "writepsf complex.psf",
        "writepdb complex.pdb",
        "exit",
    ]
    p = outdir / "psfgen.tcl"
    p.write_text("\n".join(lines) + "\n")
    return p


def psfgen_solvent(outdir: Path) -> Path:
    lines = [
        "package require psfgen",
        f"topology {HYB / 'hybrid.rtf'}",
        f"segment LIG {{ first none; last none; pdb {HYB / 'hybrid.pdb'} }}",
        f"coordpdb {HYB / 'hybrid.pdb'} LIG",
        "guesscoord",
        "writepsf ligand.psf",
        "writepdb ligand.pdb",
        "exit",
    ]
    p = outdir / "psfgen.tcl"
    p.write_text("\n".join(lines) + "\n")
    return p


def solvate(outdir: Path, psf: str, pdb: str) -> Path:
    lines = [
        "package require solvate",
        "package require autoionize",
        "mol delete all",
        f"mol load psf {psf} pdb {pdb}",
        f"solvate {psf} {pdb} -t 15 -o solvated",
        f"autoionize -psf solvated.psf -pdb solvated.pdb -neutralize -o ionized",
        "mol delete all",
        "mol load psf ionized.psf pdb ionized.pdb",
        "set all [atomselect top all]",
        "puts \"CELL [molinfo top get {a b c}]\"",
        "puts \"ORIGIN [molinfo top get {center}]\"",
        "exit",
    ]
    p = outdir / "solvate.tcl"
    p.write_text("\n".join(lines) + "\n")
    return p


def main() -> None:
    for leg, has_protein in [("complex", True), ("solvent", False)]:
        d = HERE / leg
        d.mkdir(parents=True, exist_ok=True)
        print(f"\n===== {leg.upper()} leg =====")

        # 1. psfgen
        if has_protein:
            segs = split_protein(d)
            script = psfgen_complex(d, segs)
            psf, pdb = "complex.psf", "complex.pdb"
        else:
            script = psfgen_solvent(d)
            psf, pdb = "ligand.psf", "ligand.pdb"
        out = run_vmd(script, d)
        if not (d / psf).exists():
            print(out[-3000:])
            raise SystemExit(f"psfgen failed for {leg}")
        print(f"  psfgen OK -> {psf}")

        # 2. solvate + ionize
        script = solvate(d, psf, pdb)
        out = run_vmd(script, d)
        if not (d / "ionized.psf").exists():
            print(out[-3000:])
            raise SystemExit(f"solvate failed for {leg}")
        print("  solvate/ionize OK -> ionized.psf/pdb")
        for l in out.splitlines():
            if l.startswith(("CELL ", "ORIGIN ")):
                print(f"    {l.strip()}")


if __name__ == "__main__":
    main()
