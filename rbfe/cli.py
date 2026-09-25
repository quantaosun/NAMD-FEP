"""`rbfe` -- the one command surface, the same for every system.

Every subcommand discovers its system by walking up from the current directory
for a `system.ini`, so none of them needs a path argument.  That is the whole
point: the sequence

    rbfe hybrid ; rbfe build ; rbfe inputs ; rbfe check ; rbfe run ; rbfe audit

is identical for every protein-ligand complex, which is what makes it
documentable once in RUNBOOK.md instead of per system.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from rbfe import config as C
from rbfe.errors import RbfeError

TEMPLATE = """\
# One system, one config. `rbfe init` wrote this; fill in the [system],
# [ligand], [protein] and [mutation] sections. Every other section has a
# working default -- see RUNBOOK.md Appendix A for all of them.

[system]
name        = {name}
description =

[ligand]
dir     = inputs
ref     = ref
mut     = mut
resname = UNL

[protein]
pdb            = inputs/protein.pdb
chain          = A
segment_prefix = P
first          = NTER
last           = CTER

[mutation]
# What changes between the two ligands. `vanish` are atoms only the reference
# has, `appear` only the mutant. An unrecognised strategy is a hard error --
# the engine does not guess at chemistry.
#
#   element_swap  -- one atom becomes a different element (I -> Br)
#   atom_addition -- the mutant has an atom the reference lacks (N-CH3 -> N-H)
#   mcs           -- anything else; the mapping is DERIVED from the two
#                    structures, so omit vanish/appear entirely
strategy      = element_swap
vanish        =
appear        =
placement     = native
missing_terms = synth_from_geometry

[hybrid]
title       = {name}
type_prefix = L
rtf_order   = sorted
prm_order   = source_then_new

[topology]
prot_rtf = {{root}}/../../toppar/top_all36_prot.rtf
params   = {{root}}/../../toppar/par_all36m_prot.prm
           {{root}}/../../toppar/par_all36_na.prm
           {{root}}/../../toppar/par_all36_carb.prm
           {{root}}/../../toppar/par_all36_lipid.prm
           {{root}}/../../toppar/par_all36_cgenff.prm
           {{root}}/../../toppar/par_water_ions_clean.prm
           {{hybrid}}

[run]
lambdas           = 0.0 0.045 0.09 0.14546 0.22425 0.30303 0.38182 0.46061
                    0.5394 0.61819 0.697 0.77576 0.85455 0.91 0.955 1.0
steps_per_window  = 500000

[binaries]
namd = ${{NAMD:-/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3}}
vmd  = ${{VMD:-/home/aistudio/vmd-env/bin/vmd}}
"""


def cmd_init(a) -> int:
    d = Path(a.dir).resolve()
    ini = d / "system.ini"
    if ini.exists() and not a.force:
        raise RbfeError(f"{ini} already exists (use --force to overwrite)")
    (d / "inputs").mkdir(parents=True, exist_ok=True)
    ini.write_text(TEMPLATE.format(name=d.name))
    print(f"wrote {ini}")
    print(f"wrote {d/'inputs'}/")
    print(f"""
Next:
  1. copy your inputs in:
       cp <protein>.pdb            {d.name}/inputs/protein.pdb
       cp <ref>.{{pdb,rtf,prm}}    {d.name}/inputs/
       cp <mut>.{{pdb,rtf,prm}}    {d.name}/inputs/
  2. edit {ini} -- the ligand resname, the [mutation] block, the chain
  3. cd {d.name} && rbfe hybrid""")
    return 0


def cmd_hybrid(a) -> int:
    from rbfe import hybrid as H
    from rbfe.mutations import get_strategy, spec_from_section
    cfg = C.find()
    spec = spec_from_section(cfg.mutation_section())
    strategy = get_strategy(spec.strategy).from_spec(spec)
    lig = cfg.resolve(cfg.get("ligand", "dir"))
    H.build(lig / cfg.get("ligand", "ref"), lig / cfg.get("ligand", "mut"),
            out_dir(cfg) / "hybrid",
            resname=cfg.get("ligand", "resname"), title=cfg.get("hybrid", "title"),
            type_prefix=cfg.get("hybrid", "type_prefix"),
            rtf_order=cfg.get("hybrid", "rtf_order"),
            prm_order=cfg.get("hybrid", "prm_order"),
            missing_terms=spec.missing_terms, strategy=strategy,
            similarity_threshold=(cfg.get("similarity", "threshold")
                                  if cfg.get("similarity", "enabled") else None))
    return 0


def cmd_similarity(a) -> int:
    """Score the ligand pair for RBFE suitability. Builds nothing, writes nothing.

    Always exits 0 even when the pair is below the threshold: the threshold is a
    heuristic, the operator may know better, and this reports rather than gates.
    Read the output, not the exit code.
    """
    from rbfe import similarity as S
    from rbfe.charmm import load_ligand
    cfg = C.find()
    if not cfg.get("similarity", "enabled"):
        print("similarity check disabled by [similarity] enabled = no")
        return 0
    lig = cfg.resolve(cfg.get("ligand", "dir"))
    ref = load_ligand(lig / cfg.get("ligand", "ref"))
    mut = load_ligand(lig / cfg.get("ligand", "mut"))
    result = S.check(ref, mut, threshold=cfg.get("similarity", "threshold"),
                     printer=(lambda *x, **k: None) if a.json else print)
    if a.json:
        print(json.dumps(result))
    return 0


def cmd_build(a) -> int:
    from rbfe import build as B
    B.build(C.find(), out_dir(C.find()), a.legs)
    return 0


def cmd_inputs(a) -> int:
    from rbfe import inputs as I
    cfg = C.find()
    I.write_all(cfg, out_dir(cfg), a.legs)
    return 0


def cmd_check(a) -> int:
    """Pre-flight: parse every config, then optionally time a short run.

    `--dryrun` makes NAMD read and validate the config without integrating, so
    a mistake in the generated file costs seconds rather than a GPU-hour.  The
    throughput figure comes from TIMING:, not from NAMD's cumulative
    PERFORMANCE: line, which is an average over the whole run and reads high
    early on.
    """
    from rbfe import inputs as I
    import re as _re
    cfg = C.find()
    out = out_dir(cfg)
    namd = str(cfg.resolve(cfg.get("binaries", "namd")))
    flags = list(cfg.get("binaries", "flags"))

    if not (out / "fep.tcl").exists():
        raise RbfeError(f"{out/'fep.tcl'} missing -- run `rbfe inputs` first.")
    if (out / "fep.tcl").read_text() != I.render_fep_tcl(cfg):
        raise RbfeError(
            f"{out/'fep.tcl'} does not match [run] lambdas in system.ini -- "
            f"re-run `rbfe inputs`.")

    rc = 0
    for leg in cfg.legs:
        for stage in ("nvt_equil", "npt_equil", "md_forward", "md_backward"):
            f = out / leg / f"{stage}.namd"
            if not f.exists():
                print(f"  MISSING {leg}/{stage}.namd")
                rc = 1
                continue
            r = subprocess.run([namd, *flags, "--dryrun", f"{stage}.namd"],
                               cwd=str(out / leg), capture_output=True, text=True)
            ok = r.returncode == 0
            print(f"  {'ok  ' if ok else 'FAIL'} {leg}/{stage}.namd")
            if not ok:
                print((r.stdout + r.stderr)[-1500:])
                rc = 1

    if a.smoke and rc == 0:
        leg = cfg.legs[0]
        print(f"\nsmoke run: {leg}/nvt_equil ({cfg.get('run','min_steps')} "
              f"min + a few hundred steps)")
        r = subprocess.run([namd, *flags, "nvt_equil.namd"],
                           cwd=str(out / leg), capture_output=True, text=True,
                           timeout=a.smoke_timeout)
        m = _re.findall(r"^TIMING:.*?([\d.]+)/step", r.stdout + r.stderr, _re.M)
        if not m:
            print("  no TIMING: line -- the run did not get going")
            print((r.stdout + r.stderr)[-1500:])
            return 1
        s = float(m[-1])
        ns_day = 86400 * cfg.get("run", "timestep") * 1e-6 / s
        print(f"  {s:.6f} sec/step -> {ns_day:.1f} ns/day")
        hours = cfg.nwin * 2 * len(cfg.legs) * cfg.steps_per_window * s / 3600
        print(f"  projected full job: {hours:.1f} h "
              f"({cfg.nwin} windows x 2 directions x {len(cfg.legs)} legs)")
        if cfg.get("run", "max_hours") and hours > cfg.get("run", "max_hours"):
            print(f"  NOTE: exceeds [run] max_hours = {cfg.get('run','max_hours')}")
    return rc


def _run_argv(argv: list[str]) -> int:
    """`rbfe.run` parses sys.argv itself (it was a standalone script)."""
    from rbfe import run as R
    sys.argv = ["rbfe", *argv]
    return R.main()


def cmd_run(a) -> int:
    cfg = C.find()
    if a.respawn:
        return _spawn(cfg, ["bash", str(out_dir(cfg) / "respawn_controller.sh")],
                      "controller_run.log")
    if a.watchdog:
        return _spawn(cfg, ["bash", str(out_dir(cfg) / "checkpoint_status.sh"),
                            str(a.watchdog)], "checkpoint_status.log")
    from rbfe import run as R
    R.configure(cfg)
    R.check_fep_tcl(cfg)

    # The leg order is fixed and is not the same as cfg.legs: both directions of
    # the complex must finish before the solvent starts, because the solvent leg
    # is the cheap one and a crash before it costs least.
    order = ([(a.leg, a.direction)] if a.leg
             else [(leg, d) for leg in cfg.legs for d in ("forward", "backward")])
    for leg, direction in order:
        rc = _run_argv(["run", leg, direction])
        if rc:
            return rc
    return 0


def _spawn(cfg, cmd, logname) -> int:
    out = out_dir(cfg)
    log = out / logname
    with log.open("a") as fh:
        p = subprocess.Popen(cmd, cwd=str(out), stdout=fh, stderr=subprocess.STDOUT,
                             start_new_session=True)
    print(f"started pid {p.pid}, logging to {log}")
    return 0


def cmd_status(a) -> int:
    """How much of every direction is already on disk.

    The completion test is the per-window `.done` markers, never the existence
    of a `.fepout` -- a killed window leaves a partial `.fepout` behind, so a
    file-exists test would call the leg complete and assemble a truncated one.
    """
    import json as _json
    from types import SimpleNamespace
    from rbfe import run as R
    cfg = C.find()
    R.configure(cfg)

    rows, total_done = [], 0
    for leg in cfg.legs:
        for direction in ("forward", "backward"):
            st = SimpleNamespace(leg=leg, stage=f"md_{direction}")
            import io
            import contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                R.probe(st)
            d = _json.loads(buf.getvalue())
            total_done += d["per_window_done"]
            rows.append((leg, direction, d))

    if a.json:
        print(_json.dumps({"system": cfg.name, "nwin": cfg.nwin,
                           "windows": [{"leg": l, "direction": di, **d}
                                       for l, di, d in rows]}, indent=2))
        return 0

    print(f"{cfg.name}: {cfg.nwin} windows per direction")
    for leg, direction, d in rows:
        done = d["per_window_done"]
        bar = "#" * done + "." * (cfg.nwin - done)
        note = "" if done == cfg.nwin else f"  <- next: w{done:02d}"
        print(f"  {leg:8s} {direction:9s} [{bar}] {done}/{cfg.nwin}{note}")
    if total_done == cfg.nwin * len(cfg.legs) * 2:
        print("\nall directions complete -- run `rbfe audit`")
    return 0


def cmd_windows(a) -> int:
    from rbfe import run as R
    cfg = C.find()
    R.configure(cfg)
    return _run_argv(["windows", a.leg, a.stage])


def cmd_assemble(a) -> int:
    from rbfe import run as R
    cfg = C.find()
    R.configure(cfg)
    return _run_argv(["assemble", a.leg, a.stage])


def cmd_audit(a) -> int:
    from rbfe import analysis as A
    cfg = C.find()
    A.configure(cfg)
    argv = []
    if a.scan:
        argv.append("--scan")
    if a.skip_boot:
        argv.append("--skip-boot")
    if a.trim is not None:
        argv += ["--trim", str(a.trim)]
    if a.blocks:
        argv += ["--blocks", str(a.blocks)]
    if a.root:
        argv += ["--root", a.root]
    sys.argv = ["rbfe audit", *argv]
    return A.main()


def cmd_verify(a) -> int:
    from verify.verify import run as verify_run
    return verify_run(a.case, deep=a.deep, scratch=a.scratch)


def out_dir(cfg: C.Config) -> Path:
    """Where generated artefacts go. `--out` overrides the system root."""
    return Path(_ARGS.out).resolve() if _ARGS.out else cfg.root


_ARGS = None


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="rbfe", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="write generated files here instead of the "
                                  "system directory (used by the verifier)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add(name, fn, help_):
        p = sub.add_parser(name, help=help_)
        p.set_defaults(fn=fn)
        return p

    p = add("init", cmd_init, "scaffold a new system directory")
    p.add_argument("dir", help="directory to create")
    p.add_argument("--force", action="store_true")

    add("hybrid", cmd_hybrid, "build the dual-topology hybrid ligand")
    p = add("similarity", cmd_similarity,
            "score the ligand pair for RBFE suitability (warns, never blocks)")
    p.add_argument("--json", action="store_true")
    p = add("build", cmd_build, "psfgen + solvate + ionize the systems")
    p.add_argument("--legs", nargs="+")
    p = add("inputs", cmd_inputs, "write fep.tcl, ionized.fep and the .namd files")
    p.add_argument("--legs", nargs="+")

    p = add("check", cmd_check, "validate every config; optionally time a smoke run")
    p.add_argument("--smoke", action="store_true", help="also run a short nvt")
    p.add_argument("--smoke-timeout", type=int, default=900)

    p = add("run", cmd_run, "run production (resumable, per-window)")
    p.add_argument("leg", nargs="?", choices=["complex", "solvent"],
                   help="omit to run every leg and direction in order")
    p.add_argument("direction", nargs="?", default="forward",
                   choices=["forward", "backward"])
    p.add_argument("--respawn", action="store_true",
                   help="start the self-healing wrapper instead of running directly")
    p.add_argument("--watchdog", nargs="?", const="300", metavar="SECONDS",
                   help="start the read-only status watchdog")

    p = add("status", cmd_status, "where the run is")
    p.add_argument("--json", action="store_true")
    p = add("windows", cmd_windows, "how much of one direction is on disk")
    p.add_argument("leg", choices=["complex", "solvent"])
    p.add_argument("stage", choices=["md_forward", "md_backward"])
    p = add("assemble", cmd_assemble, "concatenate the per-window .fepout files")
    p.add_argument("leg", choices=["complex", "solvent"])
    p.add_argument("stage", choices=["md_forward", "md_backward"])

    p = add("audit", cmd_audit, "independent BAR audit with bootstrap error bars")
    p.add_argument("--scan", action="store_true")
    p.add_argument("--skip-boot", action="store_true")
    p.add_argument("--blocks", type=int)
    p.add_argument("--trim", type=int)
    p.add_argument("--root")

    p = add("verify", cmd_verify, "reproduce a frozen system, byte for byte")
    p.add_argument("--case", action="append", default=None,
                   help="case name from verify/cases (repeatable)")
    p.add_argument("--deep", action="store_true", help="also re-run VMD (slow)")
    p.add_argument("--scratch", help="write here (default: a temp dir)")
    return ap


def main(argv=None) -> int:
    global _ARGS
    ap = build_parser()
    _ARGS = ap.parse_args(argv)
    try:
        return _ARGS.fn(_ARGS)
    except RbfeError as e:
        print(str(e), file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130
