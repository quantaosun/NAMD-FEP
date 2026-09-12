#!/bin/bash
# Self-healing relauncher for run_checkpointed.sh (the 6I5I FEP controller).
#
# This box kills long jobs ~4.5 h after launch (see checkpoint_status.sh header).
# run_checkpointed.sh is resumable (marker-guarded, per-window restarts), so the
# cheapest fix is to relaunch it whenever it exits before "ALL LEGS DONE".
#
# Stop conditions:
#   * controller printed ALL LEGS DONE        -> exit 0
#   * controller exits in <120 s twice in a row (real bug, no progress possible)
#                                            -> exit 1
# Environment kills (run tens of minutes .. hours) always relaunch.
#
# Run detached:  setsid nohup bash respawn_controller.sh >/dev/null 2>&1 &
set -u
cd "$(dirname "$0")"
LOG=controller_run.log
log(){ echo "[respawn $(date '+%m-%d %H:%M:%S')] $*" >> "$LOG"; }

log "self-healing wrapper started (pid $$)"
quick=0
while :; do
  start=$SECONDS
  bash run_checkpointed.sh >> "$LOG" 2>&1
  rc=$?
  elapsed=$(( SECONDS - start ))
  if grep -q "ALL LEGS DONE" "$LOG"; then
    log "controller reports ALL LEGS DONE -> wrapper exits"
    exit 0
  fi
  if [ "$elapsed" -lt 120 ]; then quick=$(( quick + 1 )); else quick=0; fi
  if [ "$quick" -ge 2 ]; then
    log "controller exited rc=$rc after ${elapsed}s twice in a row with no progress -> aborting (real bug?)"
    exit 1
  fi
  log "controller exited rc=$rc after ${elapsed}s -> relaunching in 30 s (quick_strikes=$quick)"
  sleep 30
done
