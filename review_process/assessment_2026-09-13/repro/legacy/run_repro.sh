#!/bin/zsh
# Runs each reproduction on the OLD build and then on the NEW build, in SEPARATE Stata
# processes, strictly one at a time (shared cap: wait while >= 8 stata-mp are running),
# with a 300 s watchdog per leg. Then prints the decisive lines and cross-build checks.
# Roots come from the environment so this runs from any checkout:
#   REPO=/path/to/wooldid_dev WD=/path/to/scratch/legacy zsh run_repro.sh [A1A2 A12TTB]
: "${REPO:?set REPO=<wooldid_dev checkout root>}"
: "${WD:?set WD=<scratch workdir holding old/ and new/>}"
STATA="${STATA:-/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$WD" || exit 1
ids=(A1A2 A12TTB)
[ $# -gt 0 ] && ids=("$@")
waitcap() {
  while [ "$(ps -eo args= | grep -c '^/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')" -ge 8 ]; do
    echo "cap reached, waiting"; sleep 20
  done
}
for id in $ids; do
  for b in old new; do
    waitcap
    rm -f "repro_${id}.log" "repro_${id}_${b}.log"
    echo "== $id $b start $(date +%T)"
    "$STATA" -e do "$HERE/repro_${id}.do" "$b" "$REPO" "$WD" > /dev/null 2>&1 &
    pid=$!
    n=0
    while kill -0 $pid 2>/dev/null; do
      sleep 2; n=$((n+2))
      if [ $n -ge 300 ]; then echo "WATCHDOG: killing $id $b after 300 s"; kill $pid; break; fi
    done
    mv -f "repro_${id}.log" "repro_${id}_${b}.batch.log" 2>/dev/null
    echo "== $id $b end $(date +%T) rc-lines:"; grep -h "_RC=" "repro_${id}_${b}.log" 2>/dev/null
  done
  echo "---- $id decisive lines ----"
  grep -h "^REPRO_${id}" "repro_${id}_old.log" "repro_${id}_new.log" 2>/dev/null
done
WOOLDID_REPO="$REPO" WOOLDID_SCRATCH="$WD" python3 "$HERE/repro_cross.py" $ids
