#!/bin/bash
# Watch a running NAMD FEP production stage and snapshot a clean restart
# (coor/vel/xsc) at every lambda-window boundary. On a crash these snapshots
# let the resume driver restart from the current window instead of window 0.
# Usage: watch_windows.sh <legdir> <outname> <steps_per_window>
set -euo pipefail
LEG="$1"; OUT="$2"; SPW="$3"
DIR="$(cd "$LEG" && pwd)"
SNAP="$DIR/window_snapshots"
mkdir -p "$SNAP"
SEEN_FILE="$SNAP/.seen"
# initialise with the current window index so we only snapshot *future* boundaries
last=$(tail -1 "$DIR/$OUT.xst" 2>/dev/null | awk '{print $1}' || echo 0)
[ -f "$SEEN_FILE" ] || echo $(( last / SPW )) > "$SEEN_FILE"
seen=$(cat "$SEEN_FILE")
echo "[watch] $DIR/$OUT  steps/window=$SPW  starting from window index $seen"
# wait until the namd process is gone (we watch the newest namd3)
while pgrep -f "$OUT.namd" >/dev/null 2>&1; do
  last=$(tail -1 "$DIR/$OUT.xst" 2>/dev/null | awk '{print $1}' || echo 0)
  cur=$(( last / SPW ))
  if [ "$cur" -gt "$seen" ]; then
    # crossed boundary/boundaries: snapshot start-of-window states
    for b in $(seq $((seen+1)) $cur); do
      for ext in coor vel xsc; do
        if [ -f "$DIR/$OUT.$ext" ]; then cp "$DIR/$OUT.$ext" "$SNAP/${OUT}_w$(printf '%02d' "$b").$ext"; fi
      done
      echo "[watch] $(date '+%H:%M:%S') snapshot start-of-window $b  (step ~$(( b * SPW )))"
    done
    echo "$cur" > "$SEEN_FILE"
    seen="$cur"
  fi
  sleep 1
done
echo "[watch] process ended; last window boundary seen = $seen"
