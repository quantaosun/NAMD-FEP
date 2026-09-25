#!/usr/bin/env python3
"""
Independent audit of the 6I5I dual-topology RBFE.

Written to *cross-check* analyze_fep.py, not to replace it. It deliberately
uses a different parser and a different BAR solver, and — most importantly —
validates the whole dE-extraction chain against **NAMD's own free-energy
output** written into the .fepout.

Two things this script does that analyze_fep.py does not:

1.  **Discards the per-window equilibration samples.**  alchEquilSteps 50000 at
    2 fs = 100 ps of equilibration at each lambda, and alchOutFreq 500 = 1 ps,
    so every window carries 99 pre-equilibration rows ahead of its 401
    production rows. NAMD's own accumulator resets at stepInRun ==
    alchEquilSteps (Controller::outputFepEnergy), so its dE_avg/dG columns are
    built from production samples ONLY. analyze_fep.py averages all 500 —
    including data NAMD itself throws away.

    The production start is located from NAMD's own
    '#<N> STEPS OF EQUILIBRATION AT LAMBDA ... COMPLETED' marker rather than
    hardcoded, so this stays correct if alchEquilSteps/alchOutFreq change.

2.  **Validates** that.  For each window we recompute one-sided EXP from the
    raw dE column and compare against NAMD's own per-window dG (emitted as
    '#Free energy change for lambda window [ l1 l2 ] is <dG>' and in the dG
    column). Agreement to ~1e-6 kcal/mol means the dE column, window
    splitting, and equilibration trim are all correct.

Estimators
----------
exp : one-sided Zwanzig, log-sum-exp stable.
bar : two-state Bennett per window from a fwd/bwd pair (numpy, bisection on a
      residual that is provably monotone in dG).

Usage
-----
    python3 audit_fep.py                     # full audit + block bootstrap
    python3 audit_fep.py --skip-boot         # fast, no error bars
    python3 audit_fep.py --trim 0            # reproduce the old (untrimmed) run
    python3 audit_fep.py --scan              # ddG vs trim sensitivity

Requires numpy. Read-only.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

R = 0.001987204259                      # kcal / (mol K)
TEMP = 300.0
BETA = 1.0 / (R * TEMP)
KT = R * TEMP                           # 0.5962 kcal/mol

# The window count, the legs, the temperature and the equilibration trim all
# come from system.ini.  They used to be literals here -- NWIN = 15 in four
# separate files, and the trim 99 copy-pasted in three -- which is exactly the
# kind of constant that survives a change to the schedule and silently mis-parses
# every window.  `configure()` is called once, from main().
NWIN = 15
LEGS: tuple[str, ...] = ("complex", "solvent")
EQUIL_ROWS = 99
SYSTEM_NAME = "?"


def configure(cfg) -> None:
    """Take NWIN / TEMP / legs / trim from a loaded Config."""
    global NWIN, TEMP, BETA, KT, LEGS, EQUIL_ROWS, SYSTEM_NAME
    NWIN = cfg.nwin
    TEMP = cfg.get("run", "temperature")
    BETA = 1.0 / (R * TEMP)
    KT = R * TEMP
    LEGS = tuple(cfg.legs)
    EQUIL_ROWS = cfg.equil_rows
    SYSTEM_NAME = cfg.name

WIN_HDR = re.compile(r"LAMBDA SET TO\s+([\d.eE+-]+)\s+LAMBDA2\s+([\d.eE+-]+)")
NAMD_DG = re.compile(
    r"#Free energy change for lambda window \[\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s*\]"
    r"\s*is\s+([-\d.eE+]+)\s*;\s*net change until now is\s+([-\d.eE+]+)")
EQUIL_MARK = re.compile(r"STEPS OF EQUILIBRATION AT LAMBDA")


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------
def parse_windows(path: Path) -> list[dict]:
    """Split one .fepout into windows, using NAMD's own markers.

    Returns [{'l1','l2','all','prod_start','namd_dG'}], where 'all' is every dE
    sample in the window and 'prod_start' is the index of the first
    post-equilibration sample.

    Works for both a per-window file (one window) and a single-run file
    (NAMD writes only ONE '#NEW FEP WINDOW' line per process, so window
    boundaries must come from the '#Free energy change ...' marker).
    The '... STEPS OF EQUILIBRATION ... COMPLETED' marker terminates the
    equilibration phase, so the first sample after it is production.
    """
    windows: list[dict] = []
    cur: dict | None = None
    prod = False
    hdr_l1 = hdr_l2 = None

    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#NEW FEP WINDOW"):
                m = WIN_HDR.search(ln)
                if m:
                    hdr_l1, hdr_l2 = float(m.group(1)), float(m.group(2))
                continue

            if ln.startswith("#") and EQUIL_MARK.search(ln):
                prod = True
                continue

            if ln.startswith("#Free energy change for lambda window"):
                m = NAMD_DG.search(ln)
                if m and cur is not None and cur["all"]:
                    cur["l1"] = float(m.group(1))
                    cur["l2"] = float(m.group(2))
                    cur["namd_dG"] = float(m.group(3))
                    windows.append(cur)
                cur = None
                prod = False
                continue

            if ln.startswith("FepEnergy"):
                p = ln.split()
                if len(p) >= 7:
                    if cur is None:
                        cur = {"l1": hdr_l1, "l2": hdr_l2, "all": [],
                               "prod_start": None, "namd_dG": None}
                    if prod and cur["prod_start"] is None:
                        cur["prod_start"] = len(cur["all"])
                    cur["all"].append(float(p[6]))

    if cur is not None and cur["all"]:
        windows.append(cur)

    for w in windows:
        if w["prod_start"] is None:                 # no marker found
            w["prod_start"] = 0
        w["all"] = np.asarray(w["all"], dtype=float)
    return windows


def collect_leg(root: Path, leg: str, direction: str,
                trim: int | None = None) -> list[dict]:
    """Per-window production dE for one leg+direction.

    Prefers the per-window files. complex/forward windows 0-5 have no
    per-window file — that part of the leg was salvaged from the crashed
    single run (complex/md_forward.fepout) — so those come from the single-run
    file via its own markers.

    trim=None keeps exactly the post-equilibration samples (recommended);
    trim=<int> overrides, for the sensitivity scan.
    """
    prefix = f"md_{direction}"
    wins: list[dict | None] = []
    for i in range(NWIN):
        f = root / leg / f"{prefix}_w{i:02d}.fepout"
        if not f.exists():
            wins.append(None)
            continue
        w = parse_windows(f)
        if len(w) != 1:
            raise SystemExit(f"[err] {f}: expected 1 window, got {len(w)}")
        w[0]["src"] = f.name
        wins.append(w[0])

    if any(w is None for w in wins):
        single = root / leg / f"{prefix}.fepout"
        if not single.exists():
            raise SystemExit(f"[err] no per-window files and no {single}")
        salv = parse_windows(single)
        n_missing = sum(w is None for w in wins)
        print(f"  [{leg}/{direction}] {n_missing} windows from the salvaged "
              f"single run {single.name} ({len(salv)} complete windows)")
        for i, w in enumerate(wins):
            if w is None:
                if i >= len(salv):
                    raise SystemExit(f"[err] single run has no window {i}")
                wins[i] = salv[i]
                wins[i]["src"] = f"{single.name}[w{i}]"

    out = []
    for w in wins:
        start = w["prod_start"] if trim is None else trim
        out.append({"l1": w["l1"], "l2": w["l2"],
                    "dE": w["all"][start:], "prod_start": w["prod_start"],
                    "namd_dG": w["namd_dG"], "src": w["src"]})
    return out


# --------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------
def exp_1sided(dE: np.ndarray) -> float:
    """-kT ln <exp(-beta dE)>  (log-sum-exp stable).

    dG = -(1/beta) * [ m + ln( (1/N) sum exp(x - m) ) ],  x = -beta*dE.

    The 1/beta must multiply the WHOLE bracket: m is already in units of
    x = -beta*dE. Leaving it outside rescales m by kT and silently inflates
    every result (it violates Jensen: dG > <dE>).
    """
    x = -BETA * dE
    m = float(np.max(x))
    return float(-(m + math.log(float(np.mean(np.exp(x - m))))) / BETA)


def bar_np(fwd_ab: np.ndarray, bwd_ab: np.ndarray, iters: int = 60) -> float:
    """BAR between A and B.

    fwd_ab : dU(B)-dU(A) sampled at A
    bwd_ab : dU(B)-dU(A) sampled at B   (caller sign-corrects)

    Solves sum_A s(f - dG) = sum_B s(-(g - dG)) with s(x) = 1/(1+exp(x)).
    resid is monotonically INCREASING in dG (the A-sum rises, the B-sum
    falls), so bisection is unambiguous; the bracket is checked explicitly.
    """
    f = BETA * np.asarray(fwd_ab, dtype=float)
    g = BETA * np.asarray(bwd_ab, dtype=float)

    def resid(dG: float) -> float:
        return float(np.sum(1.0 / (1.0 + np.exp(f - dG)))
                     - np.sum(1.0 / (1.0 + np.exp(dG - g))))

    lo, hi = -200.0, 200.0
    if not (resid(lo) < 0 < resid(hi)):
        raise ValueError("BAR bracket failed — residual is not monotone")
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) / BETA


def bar_leg(wins_f: list[dict], wins_b: list[dict]) -> tuple[float, list[float]]:
    """BAR over a leg. A backward leg stores its windows high->low, so window i
    of the forward leg pairs with bwd[N-1-i]."""
    per = []
    for i, wf in enumerate(wins_f):
        wb = wins_b[len(wins_b) - 1 - i]
        # A backward run reports U(low)-U(high) sampled at the HIGH lambda, so
        # negating gives dU(B)-dU(A) sampled at B, which is what BAR wants.
        per.append(bar_np(wf["dE"], -wb["dE"]))
    return float(sum(per)), per


def exp_leg(wins: list[dict]) -> tuple[float, list[float]]:
    per = [exp_1sided(w["dE"]) for w in wins]
    return float(sum(per)), per


def bar_pair(legs: dict, leg: str) -> tuple[float, list[float]]:
    return bar_leg(legs[leg]["fwd"], legs[leg]["bwd"])


# --------------------------------------------------------------------------
# report sections
# --------------------------------------------------------------------------
def section_A(legs: dict) -> None:
    print("\n" + "=" * 74)
    print("A. PARSER VALIDATION — our EXP vs NAMD's own per-window free energy")
    print("=" * 74)
    worst = 0.0
    n = 0
    for key in ("complex", "solvent"):
        rows = []
        for i, w in enumerate(legs[key]["fwd"]):
            if w["namd_dG"] is None:
                continue
            ours = exp_1sided(w["dE"])
            rows.append((i, ours, w["namd_dG"]))
            worst = max(worst, abs(ours - w["namd_dG"]))
            n += 1
        if rows:
            print(f"\n  {key}/forward  (first 3 of {len(rows)} windows)")
            print(f"    {'win':>3} {'ours':>11} {'NAMD':>11} {'diff':>11}")
            for i, ours, ref in rows[:3]:
                print(f"    {i:3d} {ours:11.6f} {ref:11.6f} {ours - ref:11.2e}")

    print(f"\n  {n} windows compared;  max |our EXP - NAMD dG| = {worst:.3e} kcal/mol")
    if worst < 1e-3:
        print("  -> PASS. The dE column, window splitting and the equilibration")
        print("     trim all reproduce NAMD's internal estimator exactly.")
    else:
        print("  -> FAIL. The extraction chain does not match NAMD.")


def section_B(legs: dict) -> tuple[float, float, float]:
    print("\n" + "=" * 74)
    print("B. BAR from the PER-WINDOW files (independent of fep_run.py assemble)")
    print("=" * 74)
    print(f"    {'win':>3} {'lambda1':>9} {'lambda2':>9} "
          f"{'dG_cplx':>10} {'dG_solv':>10} {'diff':>9}")
    tc = ts = 0.0
    for i in range(NWIN):
        cf, cb = legs["complex"]["fwd"][i], legs["complex"]["bwd"][NWIN - 1 - i]
        sf, sb = legs["solvent"]["fwd"][i], legs["solvent"]["bwd"][NWIN - 1 - i]
        dc = bar_np(cf["dE"], -cb["dE"])
        ds = bar_np(sf["dE"], -sb["dE"])
        tc += dc
        ts += ds
        print(f"    {i:3d} {cf['l1']:9.5f} {cf['l2']:9.5f} "
              f"{dc:10.4f} {ds:10.4f} {dc - ds:9.4f}")
    print(f"\n    dG_complex = {tc:8.4f} kcal/mol")
    print(f"    dG_solvent = {ts:8.4f} kcal/mol")
    print(f"    ddG        = {tc - ts:8.4f} kcal/mol")
    return tc, ts, tc - ts


def section_C(legs: dict) -> None:
    print("\n" + "=" * 74)
    print("C. HYSTERESIS  (EXP_fwd + EXP_bwd ; 0 = perfectly reversible)")
    print("=" * 74)
    for key in ("complex", "solvent"):
        ef, _ = exp_leg(legs[key]["fwd"])
        eb, _ = exp_leg(legs[key]["bwd"])
        print(f"  {key:9s} EXP fwd {ef:+8.4f}   EXP bwd {eb:+8.4f}   "
              f"hysteresis {ef + eb:+8.4f} kcal/mol")


def section_F(legs: dict) -> None:
    print("\n" + "=" * 74)
    print("F. OVERLAP — spread of dE within a window, in kT units")
    print("=" * 74)
    print("    sigma(dE)/kT >~ 2-3 means poor overlap; those windows dominate error")
    for key in ("complex", "solvent"):
        rows = sorted(((float(np.std(w["dE"])) / KT, i, w["l1"], w["l2"])
                       for i, w in enumerate(legs[key]["fwd"])), reverse=True)
        print(f"\n  {key}: worst-overlap windows")
        for s, i, l1, l2 in rows[:4]:
            print(f"    win {i:2d}  lambda {l1:.5f}->{l2:.5f}   sigma/kT = {s:5.2f}")


def section_E(legs: dict) -> None:
    print("\n" + "=" * 74)
    print("E. STATIONARITY — BAR on the first vs second half of each window")
    print("=" * 74)
    print("    a converged, equilibrated leg sums to ~0")
    for key in ("complex", "solvent"):
        drift = 0.0
        worst = (0.0, -1)
        for i in range(NWIN):
            wf, wb = legs[key]["fwd"][i], legs[key]["bwd"][NWIN - 1 - i]
            h = len(wf["dE"]) // 2
            d = (bar_np(wf["dE"][h:], -wb["dE"][h:])
                 - bar_np(wf["dE"][:h], -wb["dE"][:h]))
            drift += d
            if abs(d) > abs(worst[0]):
                worst = (d, i)
        print(f"  {key:9s} sum(late - early) = {drift:+8.4f} kcal/mol"
              f"   (largest single window: win {worst[1]}, {worst[0]:+.4f})")


def section_D(legs: dict, blocks: int, nboot: int, seed: int = 20260912) -> None:
    print("\n" + "=" * 74)
    print(f"D. BLOCK BOOTSTRAP ERROR on ddG  (block = {blocks} samples = "
          f"{blocks} ps, {nboot} resamples)")
    print("=" * 74)
    rng = np.random.default_rng(seed)
    n = len(legs["complex"]["fwd"][0]["dE"])
    nb = n // blocks

    def resample(arr: np.ndarray) -> np.ndarray:
        idx = rng.integers(0, nb, nb)
        return np.concatenate([arr[j * blocks:(j + 1) * blocks] for j in idx])

    ddgs = np.empty(nboot)
    for b in range(nboot):
        tc = ts = 0.0
        for i in range(NWIN):
            cf, cb = legs["complex"]["fwd"][i], legs["complex"]["bwd"][NWIN - 1 - i]
            sf, sb = legs["solvent"]["fwd"][i], legs["solvent"]["bwd"][NWIN - 1 - i]
            tc += bar_np(resample(cf["dE"]), -resample(cb["dE"]))
            ts += bar_np(resample(sf["dE"]), -resample(sb["dE"]))
        ddgs[b] = tc - ts
        if (b + 1) % 100 == 0:
            print(f"    ... {b + 1}/{nboot}", flush=True)

    lo, hi = np.percentile(ddgs, [2.5, 97.5])
    print(f"\n  bootstrap mean ddG = {float(np.mean(ddgs)):+.4f} kcal/mol")
    print(f"  bootstrap SD       = {ddgs.std(ddof=1):.4f} kcal/mol")
    print(f"  95% CI             = [{lo:+.4f}, {hi:+.4f}] kcal/mol")


def run_scan(root: Path, trims) -> None:
    print("\n" + "=" * 74)
    print("SENSITIVITY — ddG vs how many per-window samples are discarded")
    print("=" * 74)
    print(f"{'trim':>5} {'n/window':>9} {'dG_cplx':>9} {'dG_solv':>9} "
          f"{'ddG':>9} {'hyst_cplx':>10} {'hyst_solv':>10}")
    for trim in trims:
        legs = {leg: {d: collect_leg(root, leg,
                                     "forward" if d == "fwd" else "backward", trim)
                      for d in ("fwd", "bwd")}
                for leg in ("complex", "solvent")}
        tc, _ = bar_pair(legs, "complex")
        ts, _ = bar_pair(legs, "solvent")
        hc = sum(exp_leg(legs["complex"][d])[0] for d in ("fwd", "bwd"))
        hs = sum(exp_leg(legs["solvent"][d])[0] for d in ("fwd", "bwd"))
        n = len(legs["complex"]["fwd"][0]["dE"])
        print(f"{trim:5d} {n:9d} {tc:9.4f} {ts:9.4f} {tc - ts:9.4f} "
              f"{hc:10.4f} {hs:10.4f}")


# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None,
                    help="directory holding the leg directories "
                         "(default: the system root from system.ini)")
    ap.add_argument("--blocks", type=int, default=25)
    ap.add_argument("--boot", type=int, default=400)
    ap.add_argument("--skip-boot", action="store_true")
    ap.add_argument("--scan", action="store_true",
                    help="ddG vs equilibration-trim sensitivity, then exit")
    ap.add_argument("--trim", type=int, default=None,
                    help="override the auto-detected production start")
    a = ap.parse_args()

    from rbfe.config import find
    cfg = find()
    configure(cfg)
    root: Path = a.root or cfg.root
    print(f"{SYSTEM_NAME} FEP audit — T = {TEMP} K, kT = {KT:.4f} kcal/mol, "
          f"{NWIN} windows")
    print(f"root: {root}")

    if a.scan:
        # Derived from the config's equilibration trim rather than a literal
        # ladder that only made sense for alchEquilSteps 50000 / alchOutFreq 500.
        e = EQUIL_ROWS
        run_scan(root, sorted({0, e // 4, e // 2, 3 * e // 4, e,
                               int(e * 1.25), int(e * 1.5), int(e * 2),
                               int(e * 2.5), int(e * 3)}))
        return 0

    legs = {}
    for leg in LEGS:
        legs[leg] = {}
        for d in ("fwd", "bwd"):
            print(f"  reading {leg}/{d} ...", flush=True)
            legs[leg][d] = collect_leg(root, leg,
                                       "forward" if d == "fwd" else "backward",
                                       a.trim)

    starts = {w["prod_start"] for leg in legs.values()
              for d in leg.values() for w in d}
    print(f"\n  production start index per window (from NAMD's marker): {sorted(starts)}")
    print(f"  samples per window: {len(legs['complex']['fwd'][0]['dE'])}")

    section_A(legs)
    section_B(legs)
    section_C(legs)
    section_F(legs)
    section_E(legs)
    if not a.skip_boot:
        section_D(legs, a.blocks, a.boot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
