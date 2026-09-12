#!/bin/bash
# Run the full 6I5I dual-topology FEP: nvt -> npt -> forward -> backward
# for both complex and solvent legs, on the V100 GPU.
set -euo pipefail

# GPU-resident build (--with-single-node-cuda). ~21x faster than the old
# Linux-x86_64-g++ build; see CLAUDE.md "NAMD GPU-resident".
NAMD=/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3

# +p1 is INTENTIONAL. In GPU-resident mode throughput scales inversely with PE
# count -- measured on this box (64651-atom complex, plain MD):
#   +p1 75.8 | +p2 59.5 | +p4 27.2 | +p8 13.2 | +p16 6.3 | +p32 2.7 ns/day
# Do not "optimize" this back to +p8.
FLAGS="+p1 +devices 0"
HERE="$(cd "$(dirname "$0")" && pwd)"

run_leg() {
    local leg="$1"
    echo ""
    echo "#################################################"
    echo "##  $leg leg"
    echo "#################################################"
    cd "$HERE/$leg"
    for stage in nvt_equil npt_equil md_forward md_backward; do
        echo "===== [$leg] $stage  ($(date '+%H:%M:%S')) ====="
        $NAMD $FLAGS "$stage.namd" 2>&1 | tee "$stage.log"
    done
    cd "$HERE"
}

run_leg complex
run_leg solvent

echo ""
echo "ALL LEGS DONE at $(date)"
