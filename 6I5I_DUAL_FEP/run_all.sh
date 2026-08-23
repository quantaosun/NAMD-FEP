#!/bin/bash
# Run the full 6I5I dual-topology FEP: nvt -> npt -> forward -> backward
# for both complex and solvent legs, on the V100 GPU.
set -euo pipefail

NAMD=/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++/namd3
FLAGS="+p8 +devices 0"
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
