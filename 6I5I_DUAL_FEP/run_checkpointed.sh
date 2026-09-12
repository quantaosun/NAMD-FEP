#!/bin/bash
# Checkpoint-safe controller for the 6I5I FEP production job.
#
# Replaces run_all.sh from this point on. Strategy:
#   * complex md_forward is still running as a single NAMD process (launched
#     before this controller); a watcher (watch_windows.sh) snapshots a clean
#     restart at every lambda-window boundary so a crash there resumes at the
#     current window, not window 0.
#   * every subsequent production leg (complex backward, solvent forward/
#     backward) runs ONE NAMD invocation per lambda window via fep_run.py,
#     each writing its own restart + marker -> re-running this script after
#     any interruption skips finished windows and continues.
#   * short equilibration stages run as plain single NAMD jobs with a .done
#     marker.
#
# Safe to re-run any number of times.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NAMD=/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3
FLAGS="+p1 +devices 0"
# fep_run.py must be invoked as a proper command (a function), NOT via a
# variable holding "python3 fep_run.py" as one word -> "command not found".
run_fep(){ python3 fep_run.py "$@"; }
# Single GPU (V100): never let a second namd3 start while any namd3 is still
# running — splitting the GPU makes the whole job SLOWER. Gate every launch.
wait_gpu_free(){  # blocks until no namd3 process is alive anywhere
  while pgrep -x namd3 >/dev/null 2>&1; do
    log "GPU busy: a namd3 is still running — waiting 15 s (no GPU splitting)"
    sleep 15
  done
}
log(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

namd_in_leg(){  # $1=leg dir name  $2=stage basename  (namd running in that leg?)
  local p cmd cwd
  for p in $(pgrep -x namd3 2>/dev/null); do
    cmd=$(tr '\0' ' ' < "/proc/$p/cmdline" 2>/dev/null) || continue
    case "$cmd" in *"$2.namd"*) ;; *) continue;; esac
    cwd=$(readlink -f "/proc/$p/cwd" 2>/dev/null) || continue
    [ "${cwd%/}" = "${HERE%/}/$1" ] && return 0
  done
  return 1
}

run_stage(){  # plain single-run equil stage, marker-guarded  ($1=leg $2=basename)
  if [ -f "$1/.$2.done" ]; then log "skip $1/$2 (done)"; return; fi
  log "run  $1/$2"
  wait_gpu_free
  ( cd "$1" && "$NAMD" $FLAGS "$2.namd" > "$2.log" 2>&1 ) || { log "FAIL $1/$2"; exit 3; }
  touch "$1/.$2.done"; log "done $1/$2"
}

run_production(){  # $1=leg $2=forward|backward
  local leg=$1 dir=$2 stage=md_$2 st full per start from
  # if the original single-run NAMD for this leg/stage is still alive, wait
  while namd_in_leg "$leg" "$stage"; do log "wait: $leg/$stage single-run in progress"; sleep 60; done
  st=$(run_fep windows "$leg" "$stage")
  full=$(echo "$st" | python3 -c "import sys,json;print(json.load(sys.stdin)['single_full_windows'])")
  # NOTE: judge completion by the per-window .done MARKERS, not by whether the
  # per-window .fepout FILE exists. A killed window leaves a partial .fepout
  # behind (e.g. solvent md_forward_w14 at 146500/250000), so a file-exists
  # test would declare the leg complete and assemble a TRUNCATED window.
  per_done=$(echo "$st" | python3 -c "import sys,json;print(json.load(sys.stdin)['per_window_done'])")
  if [ "$full" -ge 15 ]; then log "$leg/$stage: single-run already complete"; run_fep assemble "$leg" "$stage" || exit 4; return; fi
  if [ "$per_done" -ge 15 ]; then log "$leg/$stage: per-window already complete"; run_fep assemble "$leg" "$stage" || exit 4; return; fi

  start=$full; from=""
  if [ "$full" -gt 0 ] && [ "$full" -lt 15 ]; then
    w=$(printf '%02d' "$full")
    if [ -f "$leg/window_snapshots/${stage}_w${w}.coor" ]; then
      from="window_snapshots/${stage}_w${w}"          # resume at current window
    else
      log "warn: no snapshot at window $full -> restarting $dir from npt_equil (replays <= $full windows)"
      start=0
    fi
  fi
  log "$leg/$dir: (re)run windows $start..14"
  wait_gpu_free
  if [ -n "$from" ]; then run_fep run "$leg" "$dir" --start "$start" --from "$from"; else run_fep run "$leg" "$dir" --start "$start"; fi \
    || { log "FAIL $leg/$dir (per-window run) -- re-run me to resume"; exit 4; }
  run_fep assemble "$leg" "$stage" || { log "FAIL $leg/$dir (assemble)"; exit 4; }
  log "$leg/$dir: done + assembled"
}

log "checkpointed controller started (resumable; re-run me after any interruption)"
run_production complex forward
run_production complex backward
run_stage      solvent nvt_equil
run_stage      solvent npt_equil
run_production solvent forward
run_production solvent backward
log "ALL LEGS DONE"
