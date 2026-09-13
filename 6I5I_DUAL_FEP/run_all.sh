#!/bin/bash
# Run the full 6I5I dual-topology FEP: nvt -> npt -> forward -> backward
# for both complex and solvent legs, on the V100 GPU.
#
# NOTE: this is the ORIGINAL, NON-RESUMABLE launcher. For anything other than a
# completely fresh run, use run_checkpointed.sh instead -- it is marker-guarded
# and resumes at the interrupted window instead of replaying from nvt_equil.
# The guard below refuses to run this script when production state already
# exists; see the comment there.
set -euo pipefail

# GPU-resident build (--with-single-node-cuda). ~21x faster than the old
# Linux-x86_64-g++ build; see ../README.md section 8 "Performance reference".
NAMD=/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3

# +p1 is INTENTIONAL. In GPU-resident mode throughput scales inversely with PE
# count -- measured on this box (64651-atom complex, plain MD):
#   +p1 75.8 | +p2 59.5 | +p4 27.2 | +p8 13.2 | +p16 6.3 | +p32 2.7 ns/day
# Do not "optimize" this back to +p8.
FLAGS="+p1 +devices 0"
HERE="$(cd "$(dirname "$0")" && pwd)"

# ---------------------------------------------------------------------------
# SAFETY GUARD -- do not remove.
#
# This script replays EVERY stage from nvt_equil and overwrites
# md_forward/md_backward. It has no notion of what is already finished, so
# running it on a completed job silently redoes ~14 h of GPU work and clobbers
# the .fepout files the published ΔΔG was computed from.
#
# This file used to be the documented entry point while the same docs said
# "Do NOT relaunch run_all.sh" -- so make the script enforce what the prose
# only warned about.
# ---------------------------------------------------------------------------
if [ "${1:-}" != "--force" ]; then
    existing="$(find "$HERE" -maxdepth 2 \
                  \( -name '.*.done' -o -name '*_combined.fepout' \) \
                  -print -quit 2>/dev/null)"
    if [ -n "$existing" ]; then
        cat >&2 <<MSG
REFUSING TO RUN: production state already exists.

    $existing

run_all.sh is NOT resumable -- it replays every stage from nvt_equil and
overwrites md_forward/md_backward. Use the resumable controller instead:

    bash run_checkpointed.sh

Deliberately overriding (this WILL redo and overwrite finished work):

    bash run_all.sh --force
MSG
        exit 1
    fi
fi

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
