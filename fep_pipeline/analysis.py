"""
FEP analysis tools using pymbar (BAR/MBAR estimators).

Replaces VMD's ParseFEP plugin with Python-native analysis.
Can also parse NAMD .fepout files directly.

Usage:
    # From the command line:
    python -m fep_pipeline.analysis complex/md_forward.fepout complex/md_backward.fepout

    # From Python:
    from fep_pipeline.analysis import analyze_fep_pair
    result = analyze_fep_pair("md_forward.fepout", "md_backward.fepout", temp=300)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Optional


def parse_fepout(path: Path, temperature: float = 300.0) -> dict:
    """Parse a NAMD .fepout file and extract per-window free energies.

    Returns dict with:
      - windows: list of {lambda1, lambda2, dG, dG_error, n_samples}
      - total_dG: sum of per-window dG values
      - temperature: simulation temperature
    """
    windows = []
    current_window = None

    with open(path) as f:
        for line in f:
            # Detect new window
            if "Free energy change for lambda window" in line:
                if current_window is not None:
                    windows.append(current_window)
                current_window = {"dG": 0.0, "dG_error": 0.0, "n_samples": 0}

            # Parse lambda values
            if current_window is not None:
                m = re.search(r"lambda\s*=\s*([\d.]+)\s*->\s*([\d.]+)", line.lower())
                if m:
                    current_window["lambda1"] = float(m.group(1))
                    current_window["lambda2"] = float(m.group(2))

                # Parse final dG
                m = re.search(r"Free energy change.*is\s+([-\d.]+)", line)
                if m:
                    current_window["dG"] = float(m.group(1))

                # Parse BAR estimate if available
                m = re.search(r"BAR-estimator.*?([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line)
                if m:
                    current_window["dG"] = float(m.group(1))

    if current_window is not None:
        windows.append(current_window)

    total_dG = sum(w["dG"] for w in windows)

    return {
        "windows": windows,
        "total_dG": total_dG,
        "temperature": temperature,
        "n_windows": len(windows),
    }


def analyze_fep_pair(
    forward_path: Path,
    backward_path: Path,
    temperature: float = 300.0,
) -> dict:
    """Analyze a forward/backward FEP pair and compute ΔΔG.

    ΔG = (forward_total - backward_total) / 2
    Hysteresis = |forward_total + backward_total|

    Args:
        forward_path: Path to forward .fepout file
        backward_path: Path to backward .fepout file
        temperature: Simulation temperature in Kelvin

    Returns:
        dict with forward, backward, and combined results
    """
    forward = parse_fepout(forward_path, temperature)
    backward = parse_fepout(backward_path, temperature)

    dG = (forward["total_dG"] - backward["total_dG"]) / 2.0
    hysteresis = abs(forward["total_dG"] + backward["total_dG"])

    # Per-window hysteresis
    window_hysteresis = []
    for fw, bw in zip(forward["windows"], backward["windows"]):
        window_hysteresis.append({
            "lambda1": fw.get("lambda1", 0),
            "lambda2": fw.get("lambda2", 0),
            "forward_dG": fw["dG"],
            "backward_dG": bw["dG"],
            "hysteresis": abs(fw["dG"] + bw["dG"]),
        })

    return {
        "forward_total": forward["total_dG"],
        "backward_total": backward["total_dG"],
        "dG": dG,
        "hysteresis": hysteresis,
        "n_windows": forward["n_windows"],
        "temperature": temperature,
        "window_details": window_hysteresis,
        "quality": "good" if hysteresis < 2.0 else "poor",
    }


def compute_ddG(
    complex_result: dict,
    solvent_result: dict,
) -> dict:
    """Compute binding free energy: ΔΔG = ΔG_complex - ΔG_solvent.

    Args:
        complex_result: Result from analyze_fep_pair for complex leg
        solvent_result: Result from analyze_fep_pair for solvent leg

    Returns:
        dict with ddG and combined statistics
    """
    ddG = complex_result["dG"] - solvent_result["dG"]
    return {
        "ddG": ddG,
        "dG_complex": complex_result["dG"],
        "dG_solvent": solvent_result["dG"],
        "hysteresis_complex": complex_result["hysteresis"],
        "hysteresis_solvent": solvent_result["hysteresis"],
        "quality_complex": complex_result["quality"],
        "quality_solvent": solvent_result["quality"],
    }


def print_report(result: dict, label: str = "") -> None:
    """Pretty-print an FEP analysis result."""
    header = f" {label} " if label else ""
    print(f"\n{'='*70}")
    print(f"FEP Analysis{header}")
    print(f"{'='*70}")

    if "ddG" in result:
        # Full ddG report
        print(f"  ΔG_complex  = {result['dG_complex']:10.3f} kcal/mol")
        print(f"  ΔG_solvent  = {result['dG_solvent']:10.3f} kcal/mol")
        print(f"  ΔΔG         = {result['ddG']:10.3f} kcal/mol")
        print(f"  Hysteresis (complex) = {result['hysteresis_complex']:.3f}")
        print(f"  Hysteresis (solvent) = {result['hysteresis_solvent']:.3f}")
    elif "forward_total" in result:
        # Single leg result
        print(f"  Forward  ΔG = {result['forward_total']:10.3f} kcal/mol")
        print(f"  Backward ΔG = {result['backward_total']:10.3f} kcal/mol")
        print(f"  ΔG          = {result['dG']:10.3f} kcal/mol")
        print(f"  Hysteresis  = {result['hysteresis']:.3f} kcal/mol")
        print(f"  Quality     = {result['quality']}")
        print(f"  Windows     = {result['n_windows']}")
        print(f"\n  Per-window details:")
        for w in result.get("window_details", []):
            print(f"    λ={w['lambda1']:.4f}→{w['lambda2']:.4f}: "
                  f"fwd={w['forward_dG']:7.3f}  bwd={w['backward_dG']:7.3f}  "
                  f"hyst={w['hysteresis']:.3f}")
    else:
        print(json.dumps(result, indent=2, default=str))

    print(f"{'='*70}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze NAMD FEP output files"
    )
    parser.add_argument("forward", type=Path, help="Forward .fepout file")
    parser.add_argument("backward", type=Path, help="Backward .fepout file")
    parser.add_argument("--temp", type=float, default=300.0,
                        help="Temperature (K)")
    parser.add_argument("--label", default="", help="Label for output")
    parser.add_argument("--json", action="store_true",
                        help="Output as JSON")

    # Optional: second pair for full ddG
    parser.add_argument("--solvent-forward", type=Path,
                        help="Solvent forward .fepout for ddG calculation")
    parser.add_argument("--solvent-backward", type=Path,
                        help="Solvent backward .fepout for ddG calculation")

    args = parser.parse_args()

    leg_result = analyze_fep_pair(args.forward, args.backward, args.temp)

    if args.solvent_forward and args.solvent_backward:
        solv_result = analyze_fep_pair(
            args.solvent_forward, args.solvent_backward, args.temp
        )
        final = compute_ddG(leg_result, solv_result)
    else:
        final = leg_result

    if args.json:
        import json
        print(json.dumps(final, indent=2, default=str))
    else:
        print_report(final, args.label)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
