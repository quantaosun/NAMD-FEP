#!/usr/bin/env python3
"""
Analyze 6I5I dual-topology NAMD FEP output directly from .fepout files.

Reads NAMD's raw FepEnergy rows (the *_combined.fepout files built by
fep_run.py assemble). Each FepEnergy row is
    STEP  Elec(l) Elec(l+dl)  vdW(l) vdW(l+dl)  dE  dE_avg  Temp  dG
where dE = U(l+dl) - U(l) is evaluated at the window's base lambda.

Estimators
----------
exp  : one-sided Zwanzig (exponential) averaging of dE, per window.
       USES ONLY THE FORWARD DATA -> gives a preview delta-delta-G while the
       solvent backward leg is still running. This is a BIASED estimator;
       the final number should come from bar once both directions exist.
bar  : two-state Bennett acceptance ratio per window from a fwd/bwd pair.
       Needs both directions of the same leg (pure Python, no pymbar).

Usage
-----
  # EXP preview from the two FORWARD legs (no backward needed)
  python3 analyze_fep.py exp complex/md_forward_combined.fepout \
                              solvent/md_forward_combined.fepout

  # BAR from forward+backward pairs (all legs complete)
  python3 analyze_fep.py bar complex/md_forward_combined.fepout \
                              complex/md_backward_combined.fepout \
                              solvent/md_forward_combined.fepout \
                              solvent/md_backward_combined.fepout

  # single-leg EXP (just the total for one fepout)
  python3 analyze_fep.py exp complex/md_forward_combined.fepout

Temperature defaults to 300 K (match `set temp` in the .namd configs).
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

R = 0.001987204259          # kcal / (mol K)  (CODATA)
WIN_HDR = re.compile(r"LAMBDA SET TO\s+([\d.]+)\s+LAMBDA2\s+([\d.]+)")


def parse_combined(path: Path) -> list[dict]:
    """Split a NAMD fepout into windows of dE samples (kcal/mol).

    Windows are delimited by '#NEW FEP WINDOW: LAMBDA SET TO ... LAMBDA2 ...'
    lines. The first such line is often duplicated inside the copied header
    (canonical_header includes the '#NEW FEP WINDOW: LAMBDA SET TO 0 ...' line),
    so a window that ends up with zero data rows is dropped. Returns windows in
    file order, each {'l1','l2','dE':[...]}.
    """
    windows: list[dict] = []
    cur: dict | None = None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#NEW FEP WINDOW"):
                m = WIN_HDR.search(ln)
                if m:
                    if cur is not None and cur["dE"]:      # close previous
                        windows.append(cur)
                    cur = {"l1": float(m.group(1)), "l2": float(m.group(2)), "dE": []}
                else:
                    cur = None
            elif ln.startswith("FepEnergy") and cur is not None:
                p = ln.split()
                if len(p) >= 7:
                    cur["dE"].append(float(p[6]))           # col 6 = dE
    if cur is not None and cur["dE"]:
        windows.append(cur)
    return windows


def exp_window(dE: list[float], beta: float) -> tuple[float, int]:
    """One-sided Zwanzig dG = -(1/beta) ln <exp(-beta dE)> (log-sum-exp stable)."""
    x = [-beta * v for v in dE]
    m = max(x)
    # dG = -(1/beta) * ln( (1/N) sum exp(x) )
    #    = -(1/beta) * ( m + ln( sum exp(x-m) / N ) )
    return -(1.0 / beta) * (m + math.log(math.fsum(math.exp(v - m) for v in x) / len(x))), len(x)


def exp_leg(path: Path, temp: float) -> dict:
    """EXP (forward-only) analysis of a single leg."""
    beta = 1.0 / (R * temp)
    wins = parse_combined(path)
    total = 0.0
    rows = []
    for w in wins:
        dg, n = exp_window(w["dE"], beta)
        total += dg
        rows.append((w["l1"], w["l2"], dg, n))
    return {"path": str(path), "windows": rows, "total_dG": total, "temp": temp}


def bar_window(fwd: list[float], bwd: list[float], beta: float, iters: int = 200) -> float:
    """Bennett acceptance ratio between states A(low) and B(high).

    fwd: dU(B) - dU(A) sampled at A   (dE from the forward window)
    bwd: dU(B) - dU(A) sampled at B   = -1 * (dE from the backward window,
         because the backward run samples dU(A)-dU(B) at B)
    Solves sum_A s(f - dG) = sum_B s(-(g - dG)) with s(x) = 1/(1+exp(beta x)).
    Returns dG_A->B in kcal/mol.
    """
    f = [beta * v for v in fwd]
    g = [beta * v for v in bwd]

    def resid(dG: float) -> float:
        left = sum(1.0 / (1.0 + math.exp(fj - dG)) for fj in f)
        right = sum(1.0 / (1.0 + math.exp(dG - gj)) for gj in g)  # = 1/(1+exp(beta(g-dG)))
        return left - right

    lo, hi = -100.0, 100.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        # resid is monotonically INCREASING in dG (both sigmoid sums rise with
        # dG), so resid(mid) > 0 means the root lies BELOW mid -> take hi = mid.
        if resid(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) / beta


def bar_leg(fwd_path: Path, bwd_path: Path, temp: float) -> dict:
    """BAR between a forward (0->1) and backward (1->0) leg of one system."""
    beta = 1.0 / (R * temp)
    fwd = parse_combined(fwd_path)
    bwd = parse_combined(bwd_path)
    if len(fwd) != len(bwd):
        print(f"[err] fwd/bwd window count mismatch: {len(fwd)} vs {len(bwd)}", file=sys.stderr)
        sys.exit(1)
    # Backward windows are stored high->low (1.0->0) in file order, so window i
    # of the forward leg shares lambda endpoints with bwd[N-1-i]. Pair by that.
    rows, total = [], 0.0
    for i, wf in enumerate(fwd):
        wb = bwd[len(bwd) - 1 - i]
        dE_fwd = wf["dE"]                 # sampled at A = lambda_i
        # Backward fepout rows report dU(A)-dU(B) sampled at B (NAMD evaluates
        # E(l+dl)-E(l) as lambda decreases), so flip sign to get the B-side
        # forward work dU(B)-dU(A) that bar_window() expects.
        dE_bwd = [-v for v in wb["dE"]]
        dG = bar_window(dE_fwd, dE_bwd, beta)
        total += dG
        rows.append((wf["l1"], wf["l2"], dG))
    return {"fwd_path": str(fwd_path), "bwd_path": str(bwd_path),
            "windows": rows, "total_dG": total, "temp": temp}


def _print_leg(title: str, res: dict, exp: bool) -> None:
    print(f"\n=== {title}  (T = {res['temp']:.0f} K)  ===")
    hdr = f"{'lambda1':>9} {'lambda2':>9} {'dG':>10} {'n':>6}" if exp else \
          f"{'lambda1':>9} {'lambda2':>9} {'dG(BAR)':>10}"
    print(hdr)
    for r in res["windows"]:
        if exp:
            print(f"{r[0]:9.5f} {r[1]:9.5f} {r[2]:10.4f} {r[3]:6d}")
        else:
            print(f"{r[0]:9.5f} {r[1]:9.5f} {r[2]:10.4f}")
    print(f"{'TOTAL dG':>10} = {res['total_dG']:.4f} kcal/mol")


def main() -> int:
    ap = argparse.ArgumentParser(description="Analyze NAMD .fepout files (EXP/BAR)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("exp", help="one-sided Zwanzig from forward data")
    e.add_argument("files", nargs="+", type=Path,
                   help="forward .fepout(s); two = complex + solvent -> ddG")
    e.add_argument("--temp", type=float, default=300.0)

    b = sub.add_parser("bar", help="two-sided BAR from fwd+bwd pairs")
    b.add_argument("files", nargs=4, type=Path,
                   help="complex_fwd complex_bwd solvent_fwd solvent_bwd")
    b.add_argument("--temp", type=float, default=300.0)

    a = ap.parse_args()

    if a.cmd == "exp":
        results = [exp_leg(p, a.temp) for p in a.files]
        if len(results) == 1:
            labels = [f"{results[0]['path'].rsplit('/', 1)[-1]} (EXP)"]
        else:
            labels = ["Complex (forward)", "Solvent (forward)"][:len(results)]
        for lab, res in zip(labels, results):
            _print_leg(lab, res, exp=True)
        if len(results) == 2:
            ddg = results[0]["total_dG"] - results[1]["total_dG"]
            print(f"\n{'='*60}")
            print("EXP (forward-only) PREVIEW — biased, use BAR for the final value")
            print(f"  dG_complex(0->1) = {results[0]['total_dG']:10.4f} kcal/mol")
            print(f"  dG_solvent(0->1) = {results[1]['total_dG']:10.4f} kcal/mol")
            print(f"  ddG = dG_complex - dG_solvent = {ddg:10.4f} kcal/mol")
            print(f"{'='*60}")
    else:
        cplx = bar_leg(a.files[0], a.files[1], a.temp)
        solv = bar_leg(a.files[2], a.files[3], a.temp)
        _print_leg("Complex (BAR)", cplx, exp=False)
        _print_leg("Solvent (BAR)", solv, exp=False)
        ddg = cplx["total_dG"] - solv["total_dG"]
        print(f"\n{'='*60}")
        print(f"  dG_complex = {cplx['total_dG']:10.4f} kcal/mol")
        print(f"  dG_solvent = {solv['total_dG']:10.4f} kcal/mol")
        print(f"  ddG (BAR)  = {ddg:10.4f} kcal/mol")
        print(f"{'='*60}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
