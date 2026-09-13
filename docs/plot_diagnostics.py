#!/usr/bin/env python3
"""Render the convergence-check figure from the 6I5I run's own output.

The banner answers "what is the number". This answers "should you believe it",
which is a different question and needs two different checks, both per window:

  hysteresis   EXP(forward) + EXP(backward) for the same transformation. The two
               directions are independent estimates, so a perfectly reversible
               leg sums to zero. This is what the BAR estimator throws away --
               it *combines* the directions rather than comparing them.
  stationarity BAR(second half) - BAR(first half) within each window. A window
               that is still drifting sums to non-zero; a converged one is flat.

Both are recomputed here through `audit_fep`, matching its sections C and E
exactly -- including the pairing (a backward leg stores its windows high->low,
so forward window i pairs with backward window N-1-i).

    python3 docs/plot_diagnostics.py     # writes docs/diagnostics-{light,dark}.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent / "6I5I_DUAL_FEP"
sys.path.insert(0, str(WORKFLOW))
sys.path.insert(0, str(HERE))

import audit_fep  # noqa: E402
from plot_banner import (GRID_LW, LW, MS, RING, THEMES, load_legs,  # noqa: E402
                         style_axes)

BAR_W = 0.38          # group is 1.0 wide, so two bars + a 2px surface gap


def per_window_checks(legs: dict):
    """Hysteresis and stationarity for every window of every leg."""
    out = {}
    for leg in ("complex", "solvent"):
        hyst, drift = [], []
        for i in range(audit_fep.NWIN):
            # Forward window i pairs with backward window N-1-i.
            wf = legs[leg]["fwd"][i]
            wb = legs[leg]["bwd"][audit_fep.NWIN - 1 - i]
            # Hysteresis: both directions estimate the same transformation.
            hyst.append(audit_fep.exp_1sided(wf["dE"])
                        + audit_fep.exp_1sided(wb["dE"]))
            # Stationarity: BAR on the late half minus BAR on the early half.
            h = len(wf["dE"]) // 2
            drift.append(audit_fep.bar_np(wf["dE"][h:], -wb["dE"][h:])
                         - audit_fep.bar_np(wf["dE"][:h], -wb["dE"][:h]))
        out[leg] = {"hyst": hyst, "drift": drift}
    return out


def panel_bars(ax, t, checks, key, title, subtitle, ylabel):
    x = np.arange(audit_fep.NWIN)
    off = BAR_W / 2 + 0.02      # the 2px surface gap between the two bars
    for sign, leg in ((-1, "complex"), (1, "solvent")):
        ax.bar(x + sign * off, checks[leg][key], width=BAR_W,
               color=t[leg], linewidth=0, zorder=3, label=leg)
    ax.axhline(0, color=t["baseline"], linewidth=GRID_LW, zorder=4)
    # Headroom above the tallest bar so the totals line cannot land on one: the
    # single +0.258 drift bar reaches the top of its axis otherwise.
    vals = [v for leg in ("complex", "solvent") for v in checks[leg][key]]
    m = max(abs(min(vals)), abs(max(vals)))
    ax.set_ylim(min(min(vals), -0.12 * m) - 0.12 * m, max(vals) + 0.45 * m)
    ax.set_xlim(-0.7, audit_fep.NWIN - 0.3)
    ax.set_xticks(np.arange(0, audit_fep.NWIN, 2))
    ax.set_xlabel("$\\lambda$ window index  (0 → 14,  $\\lambda$ 0 → 1)",
                  color=t["secondary"], fontsize=10)
    ax.text(0, 1.13, title, transform=ax.transAxes, va="bottom", ha="left",
            color=t["ink"], fontsize=13, fontweight="semibold")
    ax.text(0, 1.06, subtitle, transform=ax.transAxes, va="bottom", ha="left",
            color=t["secondary"], fontsize=9)
    style_axes(ax, t, ylabel)
    # Totals in the corner: the per-window detail is the point, but the sums are
    # the numbers audit_fep.py reports.
    totals = " · ".join(f"{leg} {sum(checks[leg][key]):+.3f}"
                        for leg in ("complex", "solvent"))
    ax.text(0.98, 0.97, totals, transform=ax.transAxes, ha="right", va="top",
            color=t["ink"], fontsize=9.5, fontweight="semibold")
    ax.legend(loc="upper left", frameon=False, fontsize=9.5,
              labelcolor=t["secondary"], handlelength=1.2, ncol=1)


def render(theme: str, out: Path, checks: dict) -> None:
    t = THEMES[theme]
    fig = plt.figure(figsize=(10, 4.6), dpi=200, facecolor=t["surface"])
    gs = GridSpec(1, 2, wspace=0.20, left=0.085, right=0.975,
                  top=0.72, bottom=0.135)

    fig.text(0.03, 0.945, "6I5I convergence checks", color=t["ink"],
             fontsize=17, fontweight="semibold")
    fig.text(0.03, 0.878,
             "Two per-window checks the headline number does not show: do the two "
             "directions agree, and has each window stopped drifting?",
             color=t["secondary"], fontsize=10.5)

    panel_bars(fig.add_subplot(gs[0, 0]), t, checks, "hyst",
               "Hysteresis per window",
               "EXP(forward) + EXP(backward)   —   0 = perfectly reversible",
               "$\\Delta$G  (kcal/mol)")
    panel_bars(fig.add_subplot(gs[0, 1]), t, checks, "drift",
               "Stationarity per window",
               "BAR(second half) − BAR(first half)   —   0 = equilibrated",
               "$\\Delta\\Delta$G  (kcal/mol)")

    fig.savefig(out, facecolor=t["surface"])
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=WORKFLOW)
    ap.add_argument("--out-dir", type=Path, default=HERE)
    a = ap.parse_args()

    checks = per_window_checks(load_legs(a.root))
    for key, name in (("hyst", "hysteresis"), ("drift", "stationarity")):
        print(f"  {name:12s} " + "  ".join(
            f"{leg} {sum(checks[leg][key]):+.4f}" for leg in ("complex", "solvent")))

    for theme in ("light", "dark"):
        render(theme, a.out_dir / f"diagnostics-{theme}.png", checks)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
