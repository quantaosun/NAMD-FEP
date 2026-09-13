#!/usr/bin/env python3
"""Render the README banner from the 6I5I production run's own output.

The banner is not an illustration: every number in it is recomputed here from
the per-window `.fepout` files, through the *same* extraction the README's
quoted result uses (`audit_fep.collect_leg` + `bar_leg`). If the run changes,
re-running this script changes the picture.

    python3 docs/plot_banner.py            # writes docs/banner-{light,dark}.png

Two files are written because GitHub renders README images in both themes and
an automatic inversion of a light chart is not a dark chart — the dark variant
uses the dark steps of the same palette, validated separately against the dark
surface.

The 95% CI is not recomputed here (the moving-block bootstrap costs minutes and
lives in `audit_fep.py`); it is passed in, defaulting to the audited values.
`python3 6I5I_DUAL_FEP/audit_fep.py` reproduces them.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent / "6I5I_DUAL_FEP"
sys.path.insert(0, str(WORKFLOW))

import audit_fep  # noqa: E402  (needs the sys.path line above)

# The audited point estimate and its bootstrap interval. See the module
# docstring: these come from audit_fep.py section D, not from this script.
DDG_DEFAULT = -0.1062
CI_LO_DEFAULT = -0.3552
CI_HI_DEFAULT = 0.1182

# --------------------------------------------------------------------------
# palette — the reference data-viz instance, one dict per theme so the two
# renders swap in one place. Both were checked with the skill's validator:
#   validate_palette.py "#2a78d6,#eb6834" --mode light   -> all checks pass
#   validate_palette.py "#3987e5,#d95926" --mode dark    -> all checks pass
# worst adjacent CVD dE 24.7 light / 26.8 dark (target >= 8).
# --------------------------------------------------------------------------
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "ink": "#0b0b0b",
        "secondary": "#52514e",
        "muted": "#898781",
        "grid": "#e1e0d9",
        "baseline": "#c3c2b7",
        "complex": "#2a78d6",     # categorical slot 1, light step
        "solvent": "#eb6834",     # categorical slot 2, light step
    },
    "dark": {
        "surface": "#1a1a19",
        "ink": "#ffffff",
        "secondary": "#c3c2b7",
        "muted": "#898781",
        "grid": "#2c2c2a",
        "baseline": "#383835",
        "complex": "#3987e5",     # slot 1, dark step
        "solvent": "#d95926",     # slot 2, dark step
    },
}

# Mark specs from the same reference: 2px line, >=8px markers, hairline solid
# grid. The figure is 15in wide and READMEs show it at ~1500px, so 1pt here is
# 2 screen px and the spec values convert as below.
LW = 1.5              # 2px line
MS = 6.0              # 8px marker
RING = 1.4            # 2px surface ring on markers
GRID_LW = 0.7         # 1px hairline


def load_data(root: Path):
    """Per-window BAR for both legs — the same call path audit_fep.py uses."""
    legs = {}
    for leg in ("complex", "solvent"):
        legs[leg] = {d: audit_fep.collect_leg(
            root, leg, "forward" if d == "fwd" else "backward")
            for d in ("fwd", "bwd")}
    per = {}
    total = {}
    for leg in ("complex", "solvent"):
        total[leg], per[leg] = audit_fep.bar_pair(legs, leg)
    lam = [(w["l1"] + w["l2"]) / 2.0 for w in legs["complex"]["fwd"]]
    return lam, per, total


def style_axes(ax, t, ylabel, xlabel=None):
    ax.set_facecolor(t["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["baseline"])
        ax.spines[side].set_linewidth(GRID_LW)
    ax.tick_params(colors=t["muted"], labelsize=9, length=3, width=GRID_LW)
    ax.set_axisbelow(True)
    # Solid hairlines; dashed grid reads as "threshold" when it is just a grid.
    ax.grid(axis="y", color=t["grid"], linewidth=GRID_LW, linestyle="-")
    ax.set_ylabel(ylabel, color=t["secondary"], fontsize=10)


def panel_profiles(ax, t, lam, per):
    """Per-window dG for both legs: the two halves that almost cancel."""
    for leg, color, label in (("complex", t["complex"], "complex (protein + ligand)"),
                              ("solvent", t["solvent"], "solvent (ligand only)")):
        ax.plot(lam, per[leg], color=color, linewidth=LW, marker="o",
                markersize=MS, markeredgecolor=t["surface"],
                markeredgewidth=RING, label=label, zorder=3)
    ax.axhline(0, color=t["baseline"], linewidth=GRID_LW, zorder=1)
    # The sign flip is not decoration: it sits at alchElecLambdaStart, where the
    # methyl's charge has finished vanishing and the N-H's has begun appearing.
    ax.axvline(0.5, color=t["grid"], linewidth=GRID_LW, zorder=1)
    # Sit the note in the empty lower-right quadrant: the curves are up at
    # 1.5-2.2 by then, so nothing crosses the text there.
    ax.annotate("sign flip at $\\lambda$ = 0.5  (alchElecLambdaStart)",
                xy=(0.53, 0.03), xycoords=("data", "axes fraction"),
                va="bottom", ha="left", color=t["muted"], fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("$\\lambda$ (window midpoint)", color=t["secondary"], fontsize=10)
    ax.set_title("Per-window free energy change", color=t["ink"],
                 fontsize=13, fontweight="semibold", loc="left", pad=10)
    style_axes(ax, t, "$\\Delta$G per window  (kcal/mol)")
    leg = ax.legend(loc="upper left", frameon=False, fontsize=10,
                    labelcolor=t["secondary"], handlelength=1.6)
    leg.set_zorder(5)


def panel_cumulative(ax, t, lam, per, total, ci_lo, ci_hi):
    """Cumulative ddG: the answer is a near-cancellation, and it is zero."""
    diffs = [c - s for c, s in zip(per["complex"], per["solvent"])]
    xs, ys = [0.0], [0.0]
    run = 0.0
    for x, d in zip(lam, diffs):
        run += d
        xs.append(x)
        ys.append(run)
    ddg = total["complex"] - total["solvent"]

    ax.axhline(0, color=t["baseline"], linewidth=GRID_LW, zorder=1)
    ax.plot(xs, ys, color=t["complex"], linewidth=LW, marker="o",
            markersize=MS, markeredgecolor=t["surface"],
            markeredgewidth=RING, zorder=3)
    # Asymmetric CI at the endpoint; the number itself is direct-labelled.
    ax.errorbar([xs[-1]], [ddg], yerr=[[ddg - ci_lo], [ci_hi - ddg]],
                fmt="none", ecolor=t["complex"], elinewidth=LW,
                capsize=4, capthick=LW, zorder=4)
    # The value block goes in the top-left: the curve dives away from it
    # immediately, so that corner stays empty, and the endpoint is left to the
    # error bar itself rather than being crowded by a label.
    ax.text(0.03, 0.97, f"$\\Delta\\Delta$G = {ddg:+.3f} kcal/mol",
            transform=ax.transAxes, ha="left", va="top",
            color=t["ink"], fontsize=12, fontweight="semibold")
    ax.text(0.03, 0.855, f"95% CI [{ci_lo:+.3f}, {ci_hi:+.3f}]",
            transform=ax.transAxes, ha="left", va="top",
            color=t["secondary"], fontsize=10)
    ax.text(0.03, 0.775, "indistinguishable from zero",
            transform=ax.transAxes, ha="left", va="top",
            color=t["muted"], fontsize=9.5)
    ax.set_xlim(-0.02, 1.05)
    ax.set_xlabel("$\\lambda$ (window midpoint)", color=t["secondary"], fontsize=10)
    ax.set_title("Cumulative $\\Delta\\Delta$G", color=t["ink"],
                 fontsize=13, fontweight="semibold", loc="left", pad=10)
    style_axes(ax, t, "$\\Sigma\\,\\Delta\\Delta$G  (kcal/mol)")


def render(theme: str, out: Path, lam, per, total, ci_lo, ci_hi) -> None:
    t = THEMES[theme]
    fig = plt.figure(figsize=(15, 5.9), dpi=200, facecolor=t["surface"])
    gs = GridSpec(1, 2, width_ratios=[1.45, 1.0], wspace=0.20,
                  left=0.055, right=0.975, top=0.755, bottom=0.115)

    fig.text(0.055, 0.935,
             "6I5I CLK1 · dual-topology RBFE — N–CH$_3$  →  N–H (desmethyl)",
             color=t["ink"], fontsize=17, fontweight="semibold")
    fig.text(0.055, 0.875,
             "H3E “12H” vs its desmethyl analogue · 15 $\\lambda$ windows × 500 ps "
             "per direction, forward + backward · Bennett acceptance ratio",
             color=t["secondary"], fontsize=10.5)
    fig.text(0.055, 0.833,
             "plotted directly from this repository's run output "
             "(6I5I_DUAL_FEP/*/md_*.fepout)",
             color=t["muted"], fontsize=9.5)

    panel_profiles(fig.add_subplot(gs[0, 0]), t, lam, per)
    panel_cumulative(fig.add_subplot(gs[0, 1]), t, lam, per, total, ci_lo, ci_hi)

    fig.savefig(out, facecolor=t["surface"])
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=WORKFLOW,
                    help="the 6I5I_DUAL_FEP directory holding the legs")
    ap.add_argument("--out-dir", type=Path, default=HERE)
    ap.add_argument("--ddg", type=float, default=DDG_DEFAULT)
    ap.add_argument("--ci-lo", type=float, default=CI_LO_DEFAULT)
    ap.add_argument("--ci-hi", type=float, default=CI_HI_DEFAULT)
    a = ap.parse_args()

    lam, per, total = load_data(a.root)
    ddg = total["complex"] - total["solvent"]
    print(f"dG_complex = {total['complex']:+.4f}   "
          f"dG_solvent = {total['solvent']:+.4f}   ddG = {ddg:+.4f}")

    for theme in ("light", "dark"):
        render(theme, a.out_dir / f"banner-{theme}.png", lam, per, total,
               a.ci_lo, a.ci_hi)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
