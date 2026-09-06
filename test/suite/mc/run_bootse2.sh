#!/bin/bash
# run_bootse2.sh - F058 bootstrap-calibration shards. Two concurrent Stata shards at
# most (project cap: never launch while 6+ "-e do" instances are already running).
# usage: ./run_bootse2.sh <shard-name> <job> [<job> ...]     job = CFG:FIRST:LAST
set -u
MC="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$MC/../../.." && pwd)"
ADO="$REPO/wooldidinstall/wooldid.ado"
OUT="$MC/out2/bootstrap"
STATA="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
REPS="${REPS:-199}"
mkdir -p "$OUT/logs"

slot_wait() {
  local n
  while : ; do
    n=$(pgrep -fl stata-mp | grep -c "^[0-9]* /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do " || true)
    [ "$n" -lt 6 ] && break
    echo "$(date +%H:%M:%S) waiting: $n stata instances" >> "$OUT/logs/$SHARD.sched"
    sleep 60
  done
}

SHARD="$1"; shift
for job in "$@"; do
  CFG="${job%%:*}"; rest="${job#*:}"; FIRST="${rest%%:*}"; LAST="${rest##*:}"
  slot_wait
  d="$OUT/logs/${SHARD}_${CFG}"; mkdir -p "$d"
  echo "$(date +%H:%M:%S) start $CFG $FIRST..$LAST reps=$REPS" >> "$OUT/logs/$SHARD.sched"
  ( cd "$d" && "$STATA" -e do "$MC/mc_bootse2.do" "$ADO" "$CFG" "$FIRST" "$LAST" "$OUT" "$REPS" >/dev/null 2>&1 )
  echo "$(date +%H:%M:%S) end   $CFG" >> "$OUT/logs/$SHARD.sched"
done
