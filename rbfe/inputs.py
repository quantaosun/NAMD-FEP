"""Write the NAMD runtime inputs for a dual-topology FEP.

Produces, for each leg:

  ionized.fep   the built system's PDB with the alchCol B-factor column marked
                (vanishing = -1, common = 0, appearing = +1)
  nvt_equil.namd, npt_equil.namd, md_forward.namd, md_backward.namd

and once per system:

  fep.tcl       the runFEP Tcl proc, sourced by the single-run configs

Ported from `4YLJ/write_fep_inputs.py`, which was already the generalized
version of this step: it reads *which* atoms are alchemical from
`hybrid/hybrid.pdb`'s B-factor column (written by the hybrid builder) and the
periodic cell from `<leg>/box.json` (written by the build step).  The 6I5I copy
hardcoded both, which is why a new system could be silently simulated in the
wrong cell.

What is new here is that *everything* comes from `system.ini`: the lambda
schedule (which used to exist in three places per system), the step counts, the
force-field constants, and the parameter-file list.  `fep.tcl` is generated from
the config rather than kept as a hand-copied file, and carries a hash so the
driver can refuse to run against a schedule that has drifted from the config.
"""

from __future__ import annotations

import hashlib
import os
import json
from pathlib import Path

from rbfe.config import Config
from rbfe.errors import ConfigError

HERE = Path(__file__).resolve().parent

# Rows a NAMD FEP window writes per outfreq interval.

LAMBDA_PER_LINE = 8


# ---------------------------------------------------------------------------
# fep.tcl
# ---------------------------------------------------------------------------

def render_fep_tcl(cfg: Config) -> str:
    """The `runFEP` proc, generated from `[run] lambdas`.

    The single-run path (`run_all.sh`-style) sources this and lets NAMD step
    through all windows in one process.  The per-window driver does not use it
    -- it substitutes explicit `alchLambda`/`alchLambda2` -- which is why the
    schedule existing in two places was survivable but not safe.
    """
    toks = cfg.lambda_tokens
    lines = []
    for i in range(0, len(toks), LAMBDA_PER_LINE):
        chunk = " ".join(toks[i:i + LAMBDA_PER_LINE])
        cont = " \\" if i + LAMBDA_PER_LINE < len(toks) else ""
        pad = "                 " if i else ""
        lines.append(f"{pad}{chunk}{cont}")
    block = "\n".join(lines)
    return f"""\
# fep.tcl -- runFEP helper for NAMD alchemical FEP (dual topology)
# {len(toks)} lambda values ({len(toks) - 1} windows), denser near the endpoints where dE/dlambda
# changes fastest (electrostatics decoupling).
# GENERATED from [run] lambdas in system.ini by `rbfe inputs` -- edit the config,
# not this file.  schedule-hash: {schedule_hash(cfg)}
set fep_lambdas {{{block}}}

proc runFEP {{nsteps {{direction forward}}}} {{
    global fep_lambdas
    set n [llength $fep_lambdas]
    for {{set i 0}} {{$i < $n - 1}} {{incr i}} {{
        if {{$direction eq "backward"}} {{
            set lambda1 [lindex $fep_lambdas [expr {{$n - 1 - $i}}]]
            set lambda2 [lindex $fep_lambdas [expr {{$n - 2 - $i}}]]
        }} else {{
            set lambda1 [lindex $fep_lambdas $i]
            set lambda2 [lindex $fep_lambdas [expr {{$i + 1}}]]
        }}
        alchLambda  $lambda1
        alchLambda2 $lambda2
        print "FEP window $i: lambda $lambda1 -> $lambda2"
        run $nsteps
    }}
}}
"""


def schedule_hash(cfg: Config) -> str:
    """A short digest of the lambda schedule, for the drift guard."""
    return hashlib.sha256(" ".join(cfg.lambda_tokens).encode()).hexdigest()[:12]


# ---------------------------------------------------------------------------
# ionized.fep
# ---------------------------------------------------------------------------

def atom_key(line: str) -> str:
    """The atom name as a case-insensitive key.

    psfgen upper-cases atom names when it builds a segment, so an element-swap
    ligand comes back from the build with `BR1` where the hybrid builder wrote
    `Br1`.  Matching the two by exact string silently finds no appearing atom
    and writes an alchFile that mutates in one direction only -- the run looks
    fine, burns the GPU time, and yields a meaningless ddG.  Chemically
    meaningful ligand atom names are not distinguished by case alone, so folding
    case is safe; every marked name is still checked against what the build
    actually produced (see `verify_marked`).
    """
    return line[12:16].strip().upper()


def alchemical_atoms(hybrid_pdb: Path) -> tuple[set[str], set[str]]:
    """Read the -1 / +1 markers straight out of the hybrid ligand's B-factors."""
    if not hybrid_pdb.exists():
        raise ConfigError(
            f"{hybrid_pdb} not found -- run `rbfe hybrid` first. The alchemical "
            f"atoms are read from its B-factor column rather than repeated in "
            f"the config, so the two cannot disagree.")
    vanish, appear = set(), set()
    for line in hybrid_pdb.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        b = float(line[60:66])
        name = atom_key(line)
        if b < 0:
            vanish.add(name)
        elif b > 0:
            appear.add(name)
    if not vanish and not appear:
        raise ConfigError(
            f"{hybrid_pdb} marks no alchemical atoms (B-factor column all zero). "
            f"A FEP with nothing alchemical would compute a ddG of zero.")
    return vanish, appear


def write_fep(leg_dir: Path, vanish: set[str], appear: set[str],
              resname: str) -> tuple[set[str], set[str]]:
    """Mark the B-factor column of ionized.pdb and write ionized.fep."""
    src = leg_dir / "ionized.pdb"
    if not src.exists():
        raise ConfigError(f"{src} not found -- run `rbfe build` first.")
    out = []
    marked_vanish, marked_appear = set(), set()
    for line in src.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            rn = line[17:20].strip()
            name = atom_key(line)
            if rn == resname:
                if name in vanish:
                    b = -1.0
                    marked_vanish.add(name)
                elif name in appear:
                    b = 1.0
                    marked_appear.add(name)
                else:
                    b = 0.0
            else:
                b = 0.0
            line = f"{line[:60]}{b:6.2f}{line[66:]}"
        out.append(line)
    (leg_dir / "ionized.fep").write_text("\n".join(out) + "\n")
    return marked_vanish, marked_appear


def verify_marked(leg: str, resname: str, vanish: set[str], appear: set[str],
                  marked_vanish: set[str], marked_appear: set[str]) -> None:
    """Fail loudly if the hybrid and the built system disagree about the ligand.

    This is the check that was missing when psfgen's upper-casing turned `Br1`
    into `BR1`: the run started, and nothing anywhere said the appearing atom
    had never been marked.

    Only the *declared* side has to be non-empty.  A pure addition (benzene ->
    toluene: a methyl appears, nothing vanishes) is a legitimate dual-topology
    FEP and an earlier version of this check refused it; what must never pass is
    *neither* side marked, which means no perturbation at all.
    """
    missing = (vanish - marked_vanish) | (appear - marked_appear)
    if missing:
        raise ConfigError(
            f"{leg}/ionized.pdb carries no alchemical atom named "
            f"{', '.join(sorted(missing))} in residue {resname} -- the hybrid "
            f"says it should.\n"
            f"  marked   vanish {sorted(marked_vanish)}  appear {sorted(marked_appear)}\n"
            f"  expected vanish {sorted(vanish)}  appear {sorted(appear)}\n"
            f"Refusing to write an alchFile that mutates in one direction only.")
    if not marked_vanish and not marked_appear:
        raise ConfigError(
            f"{leg}: no alchemical atom was marked at all. The hybrid marks "
            f"{len(vanish)} vanish / {len(appear)} appear; none of them reached "
            f"residue {resname} in the built system.")


# ---------------------------------------------------------------------------
# *.namd
# ---------------------------------------------------------------------------

def _param(p: str, cfg: Config, hybrid_param: str) -> str:
    """One `parameters` line, with the config's placeholders resolved.

    `{root}` marks a path as *derived*, so it is expanded and lexically
    normalised to an absolute path -- which makes the emitted config correct no
    matter how deep the system directory sits, and happens to reproduce 4YLJ's
    frozen absolute paths exactly.

    A value *without* `{root}` is emitted verbatim.  That is deliberate: it is
    the byte-compatibility escape hatch for reproducing a legacy config that
    used relative paths, and it is why `[topology] params` is documented as
    "what the .namd will say" rather than "a path to resolve".
    """
    s = p.replace("{hybrid}", hybrid_param)
    if "{root}" in p:
        return os.path.normpath(s.replace("{root}", str(cfg.root)))
    return s


def render_namd(cfg: Config, leg: str, kind: str, box: dict,
                hybrid_param: str, namd_path: str) -> str:
    """Render one NAMD config.  `kind` in {nvt, npt, forward, backward}."""
    a, b, c = box["cell"]
    ox, oy, oz = box["origin"]
    temp = cfg.get("run", "temperature")
    outfreq = cfg.get("run", "outfreq")
    flags = cfg.get("binaries", "flags")

    params = "\n".join(
        f"parameters              {_param(p, cfg, hybrid_param)}"
        for p in cfg.get("topology", "params"))
    fepout = {"nvt": "nvt_equil", "npt": "npt_equil",
              "forward": "md_forward", "backward": "md_backward"}[kind]

    common = f"""\
###################################################
# {leg.upper()} leg -- {kind}
###################################################
set temp               {temp}
paraTypeCharmm         on
{params}

exclude                scaled1-4
1-4scaling             1.0

structure              ionized.psf
coordinates            ionized.pdb

# NOTE: do NOT set `temperature` in this common header. NAMD errors with
# "Cannot specify both an initial temperature and a velocity file" whenever a
# stage that reads `binvelocities` (npt/forward/backward) also sets an initial
# temperature. Only `nvt` seeds velocities (via `reinitvels`); it sets
# `temperature` in its own block below.

outputenergies         {outfreq}
outputtiming           {outfreq}
outputpressure         {outfreq}
restartfreq            {outfreq}
XSTFreq                {outfreq}
dcdfreq                {outfreq}
outputname             {fepout}
restartname            {fepout}
binaryoutput           yes
binaryrestart          yes

langevin               on
langevinTemp           $temp
langevinDamping        1.0

PME                    yes
PMEGridSpacing         1.0

cellBasisVector1       {a:.3f} 0.0 0.0
cellBasisVector2       0.0 {b:.3f} 0.0
cellBasisVector3       0.0 0.0 {c:.3f}
cellOrigin             {ox:.3f} {oy:.3f} {oz:.3f}
wrapAll                on
wrapWater              on

stepspercycle          20
switching              on
switchdist             {cfg.get('run', 'switchdist'):.1f}
cutoff                 {cfg.get('run', 'cutoff'):.1f}
pairlistdist           {cfg.get('run', 'pairlistdist'):.1f}

timestep               {cfg.get('run', 'timestep'):.1f}
rigidbonds             all
rigidtolerance         0.000001
rigiditerations        400
nonbondedFreq          1
fullElectFrequency     2

ComMotion              no

# GPU-resident integration -- requires the --with-single-node-cuda build at
# {namd_path}
# (~21x faster). Run it with {flags[0] if flags else '+p1'}, NOT +p8: throughput scales inversely with
# PE count in this mode. See "The one flag that matters most" in ../README.md.
GPUresident            {cfg.get('binaries', 'gpu_resident')}

source                 ../fep.tcl
alch                   on
alchType               fep
alchFile               ionized.fep
alchCol                B
alchOutFile            {fepout}.fepout
alchOutFreq            {outfreq}
alchVdwLambdaEnd       {cfg.get('run', 'vdw_lambda_end'):.1f}
alchElecLambdaStart    {cfg.get('run', 'elec_lambda_start'):.1f}
alchVdWShiftCoeff      {cfg.get('run', 'vdw_shift_coeff'):.1f}
alchDecouple           {cfg.get('run', 'decouple')}
alchEquilSteps         {cfg.get('run', 'alch_equil_steps')}
"""

    if kind == "nvt":
        return common + f"""\
temperature            $temp
set numSteps           {cfg.get('run', 'equil_steps')}
set numMinSteps        {cfg.get('run', 'min_steps')}
alchLambda             0.0
alchLambda2            0.0
minimize               $numMinSteps
reinitvels             $temp
run                    $numSteps
"""
    if kind == "npt":
        return common + f"""\
bincoordinates        nvt_equil.coor
binvelocities         nvt_equil.vel
extendedsystem        nvt_equil.xsc

langevinPiston         on
langevinPistonTarget   1.01325
langevinPistonPeriod   100
langevinPistonDecay    50
langevinPistonTemp     $temp
useGroupPressure       yes
useFlexibleCell        no

alchLambda             0.0
alchLambda2            0.0
run                    {cfg.get('run', 'equil_steps')}
"""
    direction = "forward" if kind == "forward" else "backward"
    return common + f"""\
bincoordinates        npt_equil.coor
binvelocities         npt_equil.vel
extendedsystem        npt_equil.xsc

langevinPiston         on
langevinPistonTarget   1.01325
langevinPistonPeriod   100
langevinPistonDecay    50
langevinPistonTemp     $temp
useGroupPressure       yes
useFlexibleCell        no

runFEP                {cfg.steps_per_window} {direction}
"""


# ---------------------------------------------------------------------------

def write_all(cfg: Config, out: Path, leg_dirs: list[str] | None = None) -> None:
    """Generate fep.tcl and every per-leg config."""
    legs = leg_dirs or cfg.legs
    hybrid_pdb = out / "hybrid" / "hybrid.pdb"
    hybrid_prm = out / "hybrid" / "hybrid.prm"
    namd_path = str(cfg.resolve(cfg.get("binaries", "namd")))
    resname = cfg.get("ligand", "resname")

    vanish, appear = alchemical_atoms(hybrid_pdb)
    print(f"alchemical atoms from {hybrid_pdb}:")
    print(f"  vanish (-1): {sorted(vanish)}")
    print(f"  appear (+1): {sorted(appear)}")

    tcl = render_fep_tcl(cfg)
    (out / "fep.tcl").write_text(tcl)
    print(f"wrote fep.tcl   ({cfg.nwin} windows, schedule-hash {schedule_hash(cfg)})")

    for leg in legs:
        d = out / leg
        boxfile = d / "box.json"
        if not boxfile.exists():
            raise ConfigError(
                f"{boxfile} missing -- run `rbfe build` first. The cell is read "
                f"from there rather than from the config, so a system cannot be "
                f"simulated in another system's box.")
        box = json.loads(boxfile.read_text())
        mv, ma = write_fep(d, vanish, appear, resname)
        verify_marked(leg, resname, vanish, appear, mv, ma)
        print(f"\n{leg}: box {[round(x, 3) for x in box['cell']]} "
              f"origin {[round(x, 3) for x in box['origin']]}  (from box.json)")
        print(f"  wrote {leg}/ionized.fep  "
              f"(marked -1 {sorted(mv)}, +1 {sorted(ma)})")
        for kind, name in (("nvt", "nvt_equil"), ("npt", "npt_equil"),
                           ("forward", "md_forward"), ("backward", "md_backward")):
            cfg_text = render_namd(cfg, leg, kind, box, str(hybrid_prm), namd_path)
            (d / f"{name}.namd").write_text(cfg_text)
            print(f"  wrote {leg}/{name}.namd")
