#!/bin/bash
# run_amese.sh — MC-truth runs for the AME-SE diagnosis grid (26 cases x 256 reps).
# Structure-dead cases (d3_*, and the within-on-controls duplicates) are documented
# from the deterministic structure pass and excluded here.
set -u
MC="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$MC/../../.." && pwd)"
ADO="$REPO/wooldidinstall/wooldid.ado"
OUT="$MC/out"
STATA_BIN="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
JOBS="${STATA_JOBS:-8}"
mkdir -p "$OUT/logs"

pool_wait() { while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done; }

CASES="d1_basic_ci d1_basic_cu d1_basic_rb d1_bycoh_ci d1_bycoh_cu d1_bycoh_rb \
d1_bytime_ci d1_bytime_cu d1_bytime_rb \
d2_basic_ci d2_basic_cu d2_basic_rb d2_bycoh_ci d2_bycoh_cu d2_bycoh_rb \
d2_bytime_ci d2_bytime_cu d2_bytime_rb \
d4_within_ci d4_within_cu d4_within_rb d5_within_ci d5_within_cu d5_within_rb \
d1c_basic_ci d1c_basic_rb"

ci=0
for c in $CASES; do
  ci=$((ci + 1))
  base=$((70000000 + ci * 100000))
  for s in 1 2 3 4 5 6 7 8; do
    pool_wait
    (
      d="$OUT/logs/${c}_s${s}"
      mkdir -p "$d"
      cd "$d" && "$STATA_BIN" -e do "$MC/mc_amese.do" "$ADO" "$c" "$s" 32 "$base" "$OUT" >/dev/null 2>&1
    ) &
  done
done
wait
echo "== amese MC complete =="
fail=0
for c in $CASES; do
  n=$(cat "$OUT"/${c}_s*.csv 2>/dev/null | grep -vc "^case," || true)
  echo "$c: $n rows"
  [ "$n" -eq 256 ] || fail=1
done
exit $fail
