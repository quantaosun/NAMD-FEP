#!/usr/bin/env python3
"""
Write the NAMD FEP runtime inputs for the 6I5I dual-topology FEP:

  *.fep      solvated system PDB with the alchFile B-factor column set
             (methyl = -1, common = 0, N-H = +1; environment = 0)
  fep.tcl    the runFEP / runFEPmin Tcl procs (sourced by every .namd)
  *.namd     nvt / npt / md_forward / md_backward for complex + solvent
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
TOPPAR = HERE.parent / "toppar"

# atoms that vanish (-1) / appear (+1) in the UNL ligand
VANISH = {"C12", "H7", "H8", "H9"}   # the N-methyl
APPEAR = {"H17"}                     # the N-H

# CHARMM36 + ligand parameter files (paths relative to each leg directory)
PARAM_FILES = [
    "../../toppar/par_all36m_prot.prm",
    "../../toppar/par_all36_na.prm",
    "../../toppar/par_all36_carb.prm",
    "../../toppar/par_all36_lipid.prm",
    "../../toppar/par_all36_cgenff.prm",
    "../../toppar/par_water_ions_clean.prm",
    "../hybrid/hybrid.prm",
]

FEP_TCL = """\
# fep.tcl -- runFEP helper for NAMD alchemical FEP (dual topology)
# 16 lambda values (15 windows), denser near the endpoints where dE/dlambda
# changes fastest (electrostatics decoupling).
set fep_lambdas {0.0 0.045 0.09 0.14546 0.22425 0.30303 0.38182 0.46061 \\
                 0.5394 0.61819 0.697 0.77576 0.85455 0.91 0.955 1.0}

proc runFEP {nsteps {direction forward}} {
    global fep_lambdas
    set n [llength $fep_lambdas]
    for {set i 0} {$i < $n - 1} {incr i} {
        if {$direction eq "backward"} {
            set lambda1 [lindex $fep_lambdas [expr {$n - 1 - $i}]]
            set lambda2 [lindex $fep_lambdas [expr {$n - 2 - $i}]]
        } else {
            set lambda1 [lindex $fep_lambdas $i]
            set lambda2 [lindex $fep_lambdas [expr {$i + 1}]]
        }
        alchLambda  $lambda1
        alchLambda2 $lambda2
        print "FEP window $i: lambda $lambda1 -> $lambda2"
        run $nsteps
    }
}
"""


def write_fep(leg: Path) -> None:
    """Mark the B-factor column of ionized.pdb and write ionized.fep."""
    out = []
    for line in (leg / "ionized.pdb").read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            resname = line[17:20].strip()
            name = line[12:16].strip()
            if resname == "UNL":
                b = -1.0 if name in VANISH else (1.0 if name in APPEAR else 0.0)
            else:
                b = 0.0
            line = f"{line[:60]}{b:6.2f}{line[66:]}"
        out.append(line)
    (leg / "ionized.fep").write_text("\n".join(out) + "\n")


def namd_config(leg: str, kind: str, box, steps: dict) -> str:
    """Render one NAMD config. kind in {nvt, npt, forward, backward}."""
    a, b, c = box["cell"]
    ox, oy, oz = box["origin"]   # box centre (matches repo recipe)
    params = "\n".join(f"parameters              {p}" for p in PARAM_FILES)
    temp = 300.0
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
temperature            $temp

outputenergies         {steps['outfreq']}
outputtiming           {steps['outfreq']}
outputpressure         {steps['outfreq']}
restartfreq            {steps['outfreq']}
XSTFreq                {steps['outfreq']}
dcdfreq                {steps['outfreq']}
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
switchdist             10.0
cutoff                 12.0
pairlistdist           14.0

timestep               2.0
rigidbonds             all
rigidtolerance         0.000001
rigiditerations        400
nonbondedFreq          1
fullElectFrequency     2

ComMotion              no

source                 ../fep.tcl
alch                   on
alchType               fep
alchFile               ionized.fep
alchCol                B
alchOutFile            {fepout}.fepout
alchOutFreq            {steps['outfreq']}
alchVdwLambdaEnd       1.0
alchElecLambdaStart    0.5
alchVdWShiftCoeff      5.0
alchDecouple           off
alchEquilSteps         {steps['alch_equil']}
"""

    if kind == "nvt":
        return common + f"""\
set numSteps           {steps['nvt']}
set numMinSteps        {steps['min']}
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
run                    {steps['npt']}
"""
    # production (forward / backward) over the 16-lambda schedule
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

runFEP                {steps['prod']} {direction}
"""


def main() -> None:
    steps = {
        "outfreq": 500,
        "alch_equil": 50000,
        "min": 5000,
        "nvt": 50000,       # 0.1 ns equilibration
        "npt": 50000,       # 0.1 ns equilibration
        "prod": 250000,     # 0.5 ns per lambda window
    }

    (HERE / "fep.tcl").write_text(FEP_TCL)
    print("wrote fep.tcl")

    boxes = {
        "complex": {"cell": (89.193, 85.407, 90.624),
                    "origin": (12.069, 43.533, 21.667)},
        "solvent": {"cell": (35.881, 38.220, 40.158),
                    "origin": (19.362, 44.809, 12.024)},
    }

    for leg, box in boxes.items():
        d = HERE / leg
        write_fep(d)
        print(f"wrote {leg}/ionized.fep")
        for kind in ["nvt", "npt", "forward", "backward"]:
            cfg = namd_config(leg, kind, box, steps)
            name = {"nvt": "nvt_equil", "npt": "npt_equil",
                    "forward": "md_forward", "backward": "md_backward"}[kind]
            (d / f"{name}.namd").write_text(cfg)
            print(f"wrote {leg}/{name}.namd")


if __name__ == "__main__":
    main()
