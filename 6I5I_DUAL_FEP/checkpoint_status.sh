#!/bin/bash
# checkpoint_status.sh [interval_seconds]
# Lightweight watchdog: every tick appends one line to checkpoint_status.log
# recording the live state (active window, step, restart-file freshness, and
# how many per-window .done markers exist). Purpose: after an abrupt server
# kill (this box dies ~4.5 h from launch), the last lines show exactly where
# the job stopped and that checkpoints were still advancing. Read-only; safe
# to run alongside the controller.
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
INTERVAL="${1:-300}"
LOG=checkpoint_status.log
while true; do
  # active NAMD stage/window (the *.namd basename of the running namd3)
  namd=$(ps -eo args | grep 'namd3 +p1' | grep -v grep \
         | grep -oE '(md_forward|md_backward|nvt_equil|npt_equil)_w?[0-9]*\.namd' | head -1)
  base="-"; step="-"; fresh="-"
  if [ -n "$namd" ]; then
    base="${namd%.namd}"
    logf=$(find complex solvent -name "${base}.log" ! -path '*window_snapshots*' 2>/dev/null | head -1)
    coorf=$(find complex solvent -name "${base}.coor" ! -path '*window_snapshots*' 2>/dev/null | head -1)
    step=$(awk '/^TIMING:/{s=$2} END{print s+0}' "$logf" 2>/dev/null)
    if [ -n "$coorf" ]; then
      age=$(( $(date +%s) - $(stat -c %Y "$coorf") ))
      fresh="${age}s"
    fi
  fi
  donec=$(find complex solvent -name '*.done' 2>/dev/null | wc -l)
  ctrl=$(pgrep -fc 'run_checkpointed.sh' 2>/dev/null || echo 0)
  echo "[$(date '+%m-%d %H:%M:%S')] ctrl=${ctrl} active=${base} step=${step} restart_age=${fresh} done_markers=${donec}" >> "$LOG"
  sleep "$INTERVAL"
done
