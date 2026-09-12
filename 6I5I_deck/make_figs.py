#!/usr/bin/env python3
"""Generate the schematic figures for the 6I5I FEP slide deck.

Static scientific diagrams (no data-viz dashboard); figures are simple, light,
labelled. Outputs PNGs into ./assets.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm

import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
os.makedirs(OUT, exist_ok=True)

# ---- palette (light surfaces) -------------------------------------------------
INK    = "#232a34"   # primary text / marks
SUB    = "#5b6572"   # secondary text
MUTE   = "#c9cfd6"   # recessive grid / connectors
BLUE   = "#1f6fb2"   # primary accent (sequential-ish / single series)
AMBER  = "#e08a1e"   # second accent
RED    = "#c0392b"   # vanish highlight (has label -> not color-alone)
GREEN  = "#1e8449"   # appear highlight
FILL   = "#eef1f4"   # light fill
BG     = "#ffffff"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "text.color": INK,
    "axes.edgecolor": MUTE,
    "axes.labelcolor": INK,
    "xtick.color": SUB,
    "ytick.color": SUB,
    "figure.facecolor": BG,
    "axes.facecolor": BG,
    "savefig.facecolor": BG,
    "mathtext.fontset": "dejavusans",
})


def eq_png(name, tex, fs=46):
    fig = plt.figure(figsize=(0.1, 0.1))
    t = fig.text(0.5, 0.5, tex, ha="center", va="center", fontsize=fs, color=INK)
    fig.canvas.draw()
    bb = t.get_window_extent()
    w, h = bb.width / fig.dpi, bb.height / fig.dpi
    fig.set_size_inches(w + 0.4, h + 0.4)
    fig.savefig(os.path.join(OUT, name), dpi=220, transparent=True,
                bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


# ==============================================================================
# 1. thermodynamic cycle
# ==============================================================================
def cycle():
    fig, ax = plt.subplots(figsize=(10.4, 7.2), dpi=200)
    ax.set_xlim(0, 10); ax.set_ylim(0, 7.2); ax.axis("off")

    TL = (0.3, 4.6); TR = (6.3, 4.6)          # complex row
    BL = (0.3, 0.6); BR = (6.3, 0.6)          # solvent row
    W, H = 3.4, 1.9

    def box(xy, title, sub, fc=FILL, ec=BLUE):
        x, y = xy
        ax.add_patch(Rectangle((x, y), W, H, fc=fc, ec=ec, lw=1.6, zorder=2,
                               facecolor=fc, edgecolor=ec))
        ax.text(x + W / 2, y + H * 0.60, title, ha="center", va="center",
                fontsize=14.5, fontweight="bold", color=INK)
        ax.text(x + W / 2, y + H * 0.22, sub, ha="center", va="center",
                fontsize=10.5, color=SUB)

    box(TL, "CLK1 · 12H (methyl)", "complex — reference  (λ=0)")
    box(TR, "CLK1 · desmethyl", "complex — mutant  (λ=1)")
    box(BL, "12H (methyl)  aq.", "in water — reference  (λ=0)")
    box(BR, "desmethyl  aq.", "in water — mutant  (λ=1)")

    # alchemical legs (ref -> mut): top and bottom arrows
    top_y = TL[1] + H / 2; bot_y = BL[1] + H / 2
    for y, lab in ((top_y, r"alchemical  $\Delta G^{\mathrm{complex}}$"),
                   (bot_y, r"alchemical  $\Delta G^{\mathrm{solvent}}$")):
        ax.add_patch(FancyArrowPatch((TL[0] + W, y), (TR[0] - 0.15, y),
                      arrowstyle="-|>", mutation_scale=26, lw=2.6, color=BLUE, zorder=3))
        ax.text(4.85, y + 0.28, lab, ha="center", fontsize=13, color=BLUE)

    # binding legs: bottom -> top (solvation -> complex), labeled binding
    for x, lab in ((TL[0] + W / 2, r"$\Delta G_{\mathrm{bind}}$ (12H)"),
                   (TR[0] + W / 2, r"$\Delta G_{\mathrm{bind}}$ (desmethyl)")):
        ax.add_patch(FancyArrowPatch((x, BL[1] + H + 0.12), (x, TL[1] - 0.12),
                      arrowstyle="-|>", mutation_scale=24, lw=2.2, color=AMBER, zorder=3))
        ax.text(x + 0.14, 3.0, lab, ha="left", fontsize=11.5, color=AMBER,
                rotation=90, va="center")

    # headline formula
    ax.text(9.9, 3.6, r"$\Delta\Delta G = \Delta G^{\mathrm{complex}} - \Delta G^{\mathrm{solvent}}$",
            ha="right", fontsize=17, color=INK)
    ax.text(9.9, 3.05, "=  ΔG_bind(desmethyl)  −  ΔG_bind(12H)",
            ha="right", fontsize=11, color=SUB)
    ax.add_patch(FancyBboxPatch((6.0, 2.15), 4.15, 2.0,
                 boxstyle="round,pad=0.18,rounding_size=0.12",
                 fc="#fdf6ec", ec=AMBER, lw=1.2, zorder=1))
    ax.text(9.9, 1.62, "Absolute binding of the two ligands cancels\n"
                       "in the difference → systematic errors reduce.",
            ha="right", fontsize=10.5, color=SUB, zorder=2)

    # corner lambda tags
    ax.text(2.0, 7.05, "Four states, two alchemical legs — no ligand needs to",
            fontsize=11.5, color=SUB, ha="left")
    ax.text(2.0, 6.78, "leaving or entering the pocket during the simulation.",
            fontsize=11.5, color=SUB, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "cycle.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# 2. ligand modification schematic (N4 methyl vs H)
# ==============================================================================
def ligand_swap():
    fig, ax = plt.subplots(figsize=(13, 6.2), dpi=200)
    ax.set_xlim(0, 13); ax.set_ylim(0, 6.2); ax.axis("off")

    def draw_site(xc, yc, r=1.0):
        # generic 5-ring with N4 at the top vertex bearing the R group
        import numpy as np
        n = 6
        ang = np.linspace(90, 90 - 360, n, endpoint=False) * np.pi / 180
        ring = [(xc + r * np.cos(a), yc + r * np.sin(a)) for a in ang]
        # N4 at index 0 (top); ring neighbours idx1 and idx5
        for i in range(n):
            p1, p2 = ring[i], ring[(i + 1) % n]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=INK, lw=2.4,
                    zorder=1, solid_capstyle="round")
        n4 = ring[0]
        # ring neighbours of N4 are scaffold atoms (one of them another ring N);
        # left unlabelled on purpose — this is a schematic of the N4 site only.
        # mark N4
        ax.text(n4[0], n4[1] + 0.10, "N4", ha="center", va="center", fontsize=14,
                fontweight="bold", color=BLUE, zorder=5)
        ax.plot([n4[0]], [n4[1]], "o", ms=12, mfc="white", mec=BLUE, mew=2,
                zorder=4)
        return n4

    # ---- REF ----
    ax.text(3.0, 5.7, "Reference — co-crystal 12H", ha="center", fontsize=15,
            fontweight="bold")
    ax.text(3.0, 5.25, "N4 is N-methylated", ha="center", fontsize=11, color=SUB)
    n4r = draw_site(3.0, 3.0, 0.95)
    # N4 at top of ring: bonds up already drawn to neighbour atoms; add CH3 upward
    ax.plot([n4r[0], n4r[0]], [n4r[1] + 0.0, n4r[1] + 1.15], color=RED, lw=3.0,
            zorder=1)
    ax.text(n4r[0], n4r[1] + 1.45, "CH$_3$", ha="center", fontsize=17,
            fontweight="bold", color=RED)
    ax.text(n4r[0], n4r[1] + 1.95, "C12 · H7–H9", ha="center", fontsize=10.5,
            color=RED)
    ax.add_patch(Rectangle((1.35, 1.55), 3.3, 3.6, fc="none", ec=RED, lw=1.4,
                           ls=(0, (4, 3))))

    # ---- MUT ----
    ax.text(10.0, 5.7, "Mutant — desmethyl analogue", ha="center", fontsize=15,
            fontweight="bold")
    ax.text(10.0, 5.25, "methyl removed → N4–H", ha="center", fontsize=11,
            color=SUB)
    n4m = draw_site(10.0, 3.0, 0.95)
    ax.plot([n4m[0], n4m[0]], [n4m[1] + 0.0, n4m[1] + 0.75], color=GREEN, lw=3.0)
    ax.text(n4m[0], n4m[1] + 1.02, "H", ha="center", fontsize=17, fontweight="bold",
            color=GREEN)
    ax.text(n4m[0], n4m[1] + 1.5, "H17", ha="center", fontsize=10.5, color=GREEN)
    ax.add_patch(Rectangle((8.3, 1.55), 3.4, 3.6, fc="none", ec=GREEN, lw=1.4,
                           ls=(0, (4, 3))))

    # arrow between
    ax.add_patch(FancyArrowPatch((5.0, 3.0), (8.1, 3.0), arrowstyle="-|>",
                 mutation_scale=34, lw=3.0, color=BLUE))
    ax.text(6.55, 3.35, "alchemical λ\n0 → 1", ha="center", fontsize=12.5,
            color=BLUE)
    ax.text(6.55, 0.35, "ligand core (furo[3,2-b]pyridine scaffold) drawn schematically; "
                        "only the N4 substituent changes.", ha="center",
            fontsize=10.5, color=SUB)

    fig.savefig(os.path.join(OUT, "ligand_swap.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# 3. GPU scaling bars  (+p sweep)
# ==============================================================================
def gpu_scaling():
    pes = ["+p1", "+p2", "+p4", "+p8", "+p16", "+p32"]
    vals = [75.8, 59.5, 27.2, 13.2, 6.3, 2.7]
    colors = [BLUE, SUB, SUB, SUB, SUB, SUB]
    fig, ax = plt.subplots(figsize=(9.6, 5.6), dpi=200)
    bars = ax.bar(range(len(pes)), vals, color=colors, width=0.62, zorder=3)
    bars[0].set_edgecolor(BLUE)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.8, f"{v:g}", ha="center", fontsize=12.5,
                fontweight="bold" if i == 0 else "normal",
                color=INK if i == 0 else SUB)
    ax.text(0, 66, "21×", ha="center", fontsize=13, fontweight="bold", color=BLUE)
    ax.set_xticks(range(len(pes))); ax.set_xticklabels(pes, fontsize=12.5)
    ax.set_ylabel("ns / day  (64651-atom complex, plain MD)", fontsize=12.5)
    ax.set_xlabel("CPU threads (PE count, GPU-resident build)", fontsize=12.5)
    ax.set_ylim(0, 88)
    ax.grid(axis="y", color=MUTE, lw=0.7, alpha=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "gpu_scaling.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# 4. lambda schedule
# ==============================================================================
def lambdas():
    lam = [0.0, 0.045, 0.09, 0.14546, 0.22425, 0.30303, 0.38182, 0.46061,
           0.5394, 0.61819, 0.697, 0.77576, 0.85455, 0.91, 0.955, 1.0]
    fig, ax = plt.subplots(figsize=(12, 4.6), dpi=200)
    for i in range(len(lam) - 1):
        w = lam[i + 1] - lam[i]
        ax.add_patch(Rectangle((lam[i], 0.06), w, 0.5, fc=BLUE, ec="none",
                               alpha=0.16 + 0.5 * (1 - w / 0.1), lw=0))
    for i, x in enumerate(lam):
        edge = MUTE if 0 < i < len(lam) - 1 else BLUE
        ax.plot([x, x], [0.02, 0.62], color=edge, lw=1.2)
        if i % 2 == 0 or i in (0, 15):
            ax.text(x, -0.34, f"{x:g}", ha="center", fontsize=10, color=SUB)
        if i % 2 == 0:
            ax.text(x, 0.75, f"{i}", ha="center", fontsize=8.5, color=MUTE)
    ax.text(0.5, 0.88, "15 windows  ·  250 000 steps (0.5 ns) each",
            ha="center", fontsize=12.5, color=INK)
    ax.text(0.5, -0.78, "window i couples λ$_i$ → λ$_{i+1}$; denser at the ends where "
                        "electrostatics / soft-core VdW change fastest",
            ha="center", fontsize=10.5, color=SUB)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.9, 1.0)
    ax.set_xlabel("alchemical coupling parameter λ", fontsize=12.5, labelpad=28)
    ax.set_yticks([]); ax.spines[["left", "right", "top"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "lambda_schedule.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# 5. runtime gantt (projected)
# ==============================================================================
def runtime():
    # (leg, stage, start_h, dur_h)
    rows = [
        ("complex", "nvt + npt", 0.0, 0.10),
        ("complex", "FEP forward", 0.10, 3.80),
        ("complex", "FEP backward", 3.90, 3.80),
        ("solvent", "nvt + npt", 7.70, 0.09),
        ("solvent", "FEP forward", 7.79, 3.10),
        ("solvent", "FEP backward", 10.89, 3.10),
    ]
    colors = {"nvt + npt": FILL, "FEP forward": BLUE, "FEP backward": AMBER}
    fig, ax = plt.subplots(figsize=(12.2, 4.6), dpi=200)
    for r_i, (leg, stage, s, d) in enumerate(rows):
        ax.barh(r_i, d, left=s, height=0.55, color=colors[stage],
                edgecolor=INK if stage == "nvt + npt" else "none", lw=0.8,
                zorder=3)
        ax.text(s + d / 2, r_i, f"{d:.2g} h", ha="center", va="center",
                fontsize=11,
                color=SUB if stage == "nvt + npt" else "white",
                fontweight="bold")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{l}  ·  {st}" for l, st, _, _ in rows], fontsize=11.5)
    ax.set_xlim(0, 15.5)
    ax.set_xticks(range(0, 16, 2))
    ax.set_xlabel("projected wall-clock (hours) on 1× V100", fontsize=12.5)
    ax.grid(axis="x", color=MUTE, lw=0.7, alpha=0.5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.invert_yaxis()
    handles = [Line2D([0], [0], color=BLUE, lw=6, label="FEP production (15 windows)"),
               Line2D([0], [0], color=AMBER, lw=6, label="backward leg (λ 1→0)"),
               Line2D([0], [0], color=INK, lw=0, marker="s", ms=8,
                      mfc=FILL, mec=INK, label="minimize + nvt + npt (~6 min)")]
    ax.legend(handles=handles, loc="lower right", fontsize=10, frameon=False)
    ax.text(15.2, -0.75, "complex 7.7 h + solvent 6.2 h  ≈  14 h  (sequential)",
            ha="right", fontsize=11.5, color=SUB)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "runtime.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# 6. equations
# ==============================================================================
eq_png("eq_zwanzig.png",
       r"$\Delta G_{A\to B} = - k_B T\, \ln \, \langle \exp[\, - (V_B - V_A)\, /\, k_B T\, ] \rangle_A$",
       fs=44)
eq_png("eq_ddg.png",
       r"$\Delta\Delta G = \Delta G^{\mathrm{complex}} - \Delta G^{\mathrm{solvent}} = \Delta G_{\mathrm{bind}}(B) - \Delta G_{\mathrm{bind}}(A)$",
       fs=40)
eq_png("eq_sum.png",
       r"$\Delta G_{0\to 1} = \sum_{i=1}^{M} \Delta G_i , \quad \Delta G_i = - \frac{1}{\beta} \ln \langle e^{-\beta \Delta V_i}\rangle_i , \quad \beta = \frac{1}{k_B T}$",
       fs=40)


if __name__ == "__main__":
    cycle()
    ligand_swap()
    gpu_scaling()
    lambdas()
    runtime()
    print("figures written to", OUT)
    for f in sorted(os.listdir(OUT)):
        print("  ", f, os.path.getsize(os.path.join(OUT, f)), "bytes")
