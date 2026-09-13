#!/usr/bin/env python3
"""
Checkpoint-safe FEP production runner for the 6I5I dual-topology job.

The stock .namd runs all 15 lambda windows in ONE NAMD process, so a crash
mid-leg makes NAMD forget the current lambda and the whole direction replays
from window 0. This driver instead runs ONE NAMD invocation per lambda window,
each writing its own restart (coor/vel/xsc) + fepout + a '.done' marker, so any
interruption costs at most the window that was in flight, and re-running the
driver simply skips finished windows.

Subcommands
-----------
run <leg> <dir> [--start N] [--from STEM] [--steps N]
    dir in {forward, backward}. Runs windows N..14, seeding window N from
    --from (a restart stem, e.g. window_snapshots/md_forward_w05) or, by
    default, chaining from npt_equil (window 0) / the previous window's own
    restart. Writes markers; skips windows already marked done.

assemble <leg> <stage>
    stage in {md_forward, md_backward}. Builds <stage>_combined.fepout:
      * per-window runs 0..14 -> concatenate in order (drop repeated headers)
      * a completed single-run <stage>.fepout with no per-window runs -> copy
      * a crashed single-run <stage>.fepout (windows 0..k-1 complete) followed
        by per-window k..14 -> salvage + append
"""
from __future__ import annotations
import argparse, json, os, pathlib, re, shutil, subprocess, sys

BASE = pathlib.Path(__file__).resolve().parent
NAMD = os.environ.get(
    "NAMD", "/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3")
# NOTE: "+devices" and "0" must be SEPARATE argv tokens. subprocess.run passes
# each list element verbatim, so "+devices 0" (one element w/ embedded space)
# is NOT parsed by the Charm++ RTS -> NAMD swallows it as the config file and
# dies with "FATAL ERROR: Unknown command-line option +devices 0". (A shell
# would word-split "+p1 +devices 0"; a list does not.)
FLAGS = ["+p1", "+devices", "0"]

# must match the schedule written into fep.tcl
LAM = [0.0, 0.045, 0.09, 0.14546, 0.22425, 0.30303, 0.38182, 0.46061,
       0.5394, 0.61819, 0.697, 0.77576, 0.85455, 0.91, 0.955, 1.0]
NWIN = len(LAM) - 1                     # 15
EXPECT_LINES = 250000 // 500            # FepEnergy rows per complete window (alchOutFreq 500)
DEFAULT_STEPS = 250000
# Of those EXPECT_LINES rows, the first EQUIL_ROWS are the per-window
# equilibration (alchEquilSteps 50000 / alchOutFreq 500, minus the row written
# at the reset step itself). NAMD's own dE_avg/dG columns exclude them, and any
# analysis must too -- see analyze_fep.py / audit_fep.py.
EQUIL_ROWS = 99


def stage_of(dir_: str) -> str:
    return "md_forward" if dir_ == "forward" else "md_backward"


def window_lambdas(dir_: str, i: int):
    if dir_ == "forward":
        return LAM[i], LAM[i + 1]
    return LAM[NWIN - i], LAM[NWIN - 1 - i]


def write_window_config(leg: pathlib.Path, stage: str, i: int, base: str,
                        steps: int) -> pathlib.Path:
    """Render a single-window config from the stock single-run .namd template."""
    cfg = (leg / f"{stage}.namd").read_text()
    out = f"{stage}_w{i:02d}"

    def repl(key: str, val: str) -> None:
        nonlocal cfg
        cfg = re.sub(rf"(?m)^({key})\s+\S+\s*$", rf"\1 {val}", cfg)

    repl("outputname", out)
    repl("restartname", out)
    repl("alchOutFile", f"{out}.fepout")
    repl("bincoordinates", f"{base}.coor")
    repl("binvelocities", f"{base}.vel")
    repl("extendedsystem", f"{base}.xsc")

    keep = [ln for ln in cfg.splitlines()
            if not ln.lstrip().startswith(("source ", "runFEP "))]
    l1, l2 = window_lambdas(stage.replace("md_", ""), i)
    keep += [f"alchLambda             {l1}",
             f"alchLambda2            {l2}",
             f"run                    {steps}"]
    out_path = leg / f"{out}.namd"
    out_path.write_text("\n".join(keep) + "\n")
    return out_path


def run_window(leg: pathlib.Path, stage: str, i: int, base: str,
               steps: int, dry: bool) -> None:
    marker = leg / f".{stage}_w{i:02d}.done"
    if marker.exists():
        print(f"[skip] {stage}_w{i:02d} already done")
        return
    cfg = write_window_config(leg, stage, i, base, steps)
    log = leg / f"{stage}_w{i:02d}.log"
    print(f"[run ] {stage}_w{i:02d}  lambda {window_lambdas(stage.replace('md_',''), i)[0]} "
          f"-> {window_lambdas(stage.replace('md_',''), i)[1]}  restart={base}", flush=True)
    if dry:
        # Deliberately do NOT touch the marker. A 'dry run' that marks work as
        # done poisons the job: every later real run skips those windows, and
        # run_checkpointed.sh sees per_window_done == 15 and declares the leg
        # complete. Must be the only state --dry does not change.
        print(f"       (dry) would run: {NAMD} {' '.join(FLAGS)} {cfg.name}")
        print(f"       (dry) marker {marker.name} NOT written")
        return
    with log.open("w") as fh:
        r = subprocess.run([NAMD, *FLAGS, cfg.name], cwd=leg, stdout=fh,
                           stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"[FAIL] {stage}_w{i:02d} exit={r.returncode} (see {log}) — re-run to resume", flush=True)
        sys.exit(2)
    marker.touch()
    print(f"[done] {stage}_w{i:02d}", flush=True)


def bucket_of(step: int) -> int:
    return (step - 1) // DEFAULT_STEPS          # which lambda window a step belongs to


def _is_column_header(line: str) -> bool:
    """True for the two '#...STEP...Elec...dG' / '#...E(l+dl)-E(l)' headers.

    The 'Elec' test matters: NAMD's equilibration marker is
    '#50000 STEPS OF EQUILIBRATION AT LAMBDA 0.5394 COMPLETED', which also
    contains 'STEP'. Matching on 'STEP' alone silently deletes that marker.
    """
    return (line.startswith("#")
            and ("E(l+dl)-E(l)" in line
                 or ("STEP" in line and "Elec" in line)))


def _keep(line: str) -> bool:
    """Lines worth carrying into the combined file.

    KEEPS NAMD's '#<N> STEPS OF EQUILIBRATION ... COMPLETED' and '#Free energy
    change for lambda window ... is <dG>' comments. They are the only record of
    where production sampling starts and they let the result be validated
    against NAMD's own estimator. An earlier version of this function kept only
    'FepEnergy' lines, which is exactly why an analysis that averaged the
    per-window equilibration rows went unnoticed for three days.
    """
    if _is_column_header(line):
        return False
    if line.startswith("#NEW FEP WINDOW"):
        return False                    # assemble writes a correct one per window
    return True


def per_window_lines(path: pathlib.Path) -> list[str]:
    """Content lines of a per-window .fepout (fresh process, single window)."""
    return [ln for ln in path.read_text().splitlines() if _keep(ln)]


def single_window_data(path: pathlib.Path):
    """Split a (single-run) .fepout into per-window lines by step bucket.

    NAMD emits the '#NEW FEP WINDOW' comment once per process, but steps are
    continuous across windows, so windows are delimited by step range:
    window i covers (i*SPW, (i+1)*SPW]. Returns (idx->[lines], n_complete).

    Comment lines are attached to the window of the most recent FepEnergy row,
    which is where NAMD wrote them (the '<N> STEPS OF EQUILIBRATION' marker for
    window i is emitted at step i*SPW + alchEquilSteps, i.e. mid-window).
    """
    data: dict[int, list[str]] = {}
    if not path.exists():
        return data, 0
    cur: int | None = None
    for line in path.read_text().splitlines():
        if line.startswith("FepEnergy"):
            parts = line.split()
            if len(parts) > 1 and parts[1].isdigit():
                cur = bucket_of(int(parts[1]))
                data.setdefault(cur, []).append(line)
            continue
        if not _keep(line) or not line.startswith("#"):
            continue
        if cur is not None:
            data.setdefault(cur, []).append(line)
    have = set()
    for i, lns in data.items():
        if any(int(ln.split()[1]) == (i + 1) * DEFAULT_STEPS
               for ln in lns if ln.startswith("FepEnergy")):
            have.add(i)                          # window closed on its boundary step
    full = 0
    for i in range(len(data)):
        if i in have:
            full += 1
        else:
            break
    return data, full


def canonical_header() -> list[str]:
    """The two column-header lines for the combined file.

    Deliberately NOT a '#NEW FEP WINDOW' line: assemble writes a correct one
    per window below. The previous version copied the first three lines of the
    single-run .fepout, which injected a bogus 'LAMBDA SET TO 0 LAMBDA2 0.045'
    header into the BACKWARD combined file too.
    """
    return [
        "#            STEP                 Elec                            vdW                    dE           dE_avg             Temp             dG",
        "#                           l             l+dl             l            l+dl         E(l+dl)-E(l)",
    ]


def run_cmd(args: argparse.Namespace) -> None:
    leg = BASE / args.leg
    stage = stage_of(args.dir)
    start = args.start
    base = args.from_stem if (args.from_stem and start is not None) else None
    for i in range(start, NWIN):
        if i == start:
            if base is None:
                base = "npt_equil"
        else:
            base = f"{stage}_w{i-1:02d}"
        run_window(leg, stage, i, base, args.steps, args.dry)


def assemble(args: argparse.Namespace) -> None:
    leg = BASE / args.leg
    stage = args.stage
    direction = stage.replace("md_", "")
    single, full = single_window_data(leg / f"{stage}.fepout")
    out = leg / f"{stage}_combined.fepout"

    # pick data for every window i: a COMPLETE per-window run (has both the
    # .fepout and its .done marker), else the single-run step-bucket.
    lines: dict[int, list[str]] = {}
    for i in range(NWIN):
        wf = leg / f"{stage}_w{i:02d}.fepout"
        marker = leg / f".{stage}_w{i:02d}.done"
        if wf.exists() and marker.exists():
            lines[i] = per_window_lines(wf)
        elif wf.exists():
            # A killed window leaves a partial .fepout behind. Gating on the
            # .fepout alone would silently assemble a TRUNCATED window, so
            # require the marker (which run_window only writes on exit 0).
            print(f"[err] {wf.name} exists but {marker.name} is missing: that "
                  f"window was interrupted and its .fepout is PARTIAL.\n"
                  f"      Delete the partial file and re-run:\n"
                  f"        python3 fep_run.py run {args.leg} {direction} --start {i}")
            sys.exit(1)
        elif i < full:
            lines[i] = single[i]
        else:
            print(f"[err] window {i} has no complete data — run the driver first "
                  f"(python3 fep_run.py run {args.leg} {direction})")
            sys.exit(1)

    # Refuse to assemble a short window: truncation is invisible downstream,
    # it just quietly biases the estimate.
    for i in range(NWIN):
        n = sum(1 for ln in lines[i] if ln.startswith("FepEnergy"))
        if n != EXPECT_LINES:
            print(f"[err] window {i}: {n} FepEnergy rows, expected {EXPECT_LINES} "
                  f"— refusing to assemble a truncated window")
            sys.exit(1)

    out_lines = canonical_header()
    for i in range(NWIN):
        l1, l2 = window_lambdas(direction, i)
        out_lines.append(f"#NEW FEP WINDOW: LAMBDA SET TO {l1} LAMBDA2 {l2}")
        out_lines.extend(lines[i])
    out.write_text("\n".join(out_lines).rstrip("\n") + "\n")
    src = "per-window" if any((leg / f"{stage}_w{i:02d}.fepout").exists()
                              for i in range(NWIN)) else "single-run"
    print(f"[ok] assembled {out} ({NWIN} windows x {EXPECT_LINES} rows, "
          f"source: {src}); NAMD equilibration / free-energy markers preserved")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run")
    r.add_argument("leg", choices=["complex", "solvent"])
    r.add_argument("dir", choices=["forward", "backward"])
    r.add_argument("--start", type=int, default=0)
    r.add_argument("--from", dest="from_stem")
    r.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    r.add_argument("--dry", action="store_true")
    r.set_defaults(fn=run_cmd)

    a = sub.add_parser("assemble")
    a.add_argument("leg", choices=["complex", "solvent"])
    a.add_argument("stage", choices=["md_forward", "md_backward"])
    a.set_defaults(fn=assemble)

    w = sub.add_parser("windows")
    w.add_argument("leg", choices=["complex", "solvent"])
    w.add_argument("stage", choices=["md_forward", "md_backward"])
    w.set_defaults(fn=probe)

    args = ap.parse_args()
    args.fn(args)


def probe(args: argparse.Namespace) -> None:
    """Report how much of a direction is already on disk (JSON)."""
    leg = BASE / args.leg
    stage = args.stage
    _data, full = single_window_data(leg / f"{stage}.fepout")
    per_files = sum((leg / f"{stage}_w{i:02d}.fepout").exists() for i in range(NWIN))
    per_done = sum((leg / f".{stage}_w{i:02d}.done").exists() for i in range(NWIN))
    snapdir = leg / "window_snapshots"
    snap = sorted(p.name for p in snapdir.glob(f"{stage}_w*.coor")) if snapdir.exists() else []
    print(json.dumps({
        "single_full_windows": full,     # complete windows in the single-run .fepout
        "per_window_files": per_files,
        "per_window_done": per_done,
        "snapshots": snap,
    }))


if __name__ == "__main__":
    main()
