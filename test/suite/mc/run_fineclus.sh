#!/bin/bash
# run_fineclus.sh - F3 shards. Modeled on test/suite/mc/run_bootse2.sh.
# usage: ./run_fineclus.sh <reps> <job> [<job> ...]     job = CFG:FIRST:LAST
# Launches every job CONCURRENTLY, one Stata each. The caller is responsible for the
# concurrency cap: the project ceiling is 8 concurrent Stata instances, and the F3
# readout was produced with 4 (FINE:1:100 FINE:101:200 NEST:1:100 NEST:101:200 at
# reps 199). Outputs land in out_fineclus/ and are deliberately not tracked.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ADO="${ADO:-$(cd "$HERE/../../.." && pwd)/wooldidinstall/wooldid.ado}"
STATA="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
OUT="$HERE/out_fineclus"
mkdir -p "$OUT"

REPS="$1"; shift
for job in "$@"; do
  CFG="${job%%:*}"; rest="${job#*:}"; FIRST="${rest%%:*}"; LAST="${rest##*:}"
  d="$OUT/${CFG}_${FIRST}_${LAST}"; mkdir -p "$d"
  cp "$HERE/mc_fineclus.do" "$HERE/mc_fineclus_dgp.doh" "$d/"
  ( cd "$d" && "$STATA" -e do mc_fineclus.do "$ADO" "$CFG" "$FIRST" "$LAST" "$OUT" "$REPS" >/dev/null 2>&1 ) &
  echo "launched $CFG $FIRST..$LAST reps=$REPS -> $d"
done
wait
echo "ALL_SHARDS_DONE"
