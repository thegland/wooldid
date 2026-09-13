#!/bin/bash
# run_advb.sh — adversarial-boundary (AB) MC orchestrator.
# PREREG_ADVBOUND_2026-08-25: maps where the adversarial line test is correctly sized.
#
# The grid lives in advb_grid.json (committed with the preregistration), so this
# runner and the driver are fixed before any rep fires. usage:
#   ./run_advb.sh              run every config in the grid, then analyze
#   ./run_advb.sh AB1 AB2      run only these configs, then analyze
#   ./run_advb.sh analyze      re-judge existing CSVs without running
#   STATA_JOBS=8               pool size (default 8 = the project hard cap)
#   ADVB_OUT=/some/dir         output tree (default out/ = the 2026-08-26 record; see below)
set -u
MC="$(cd "$(dirname "$0")" && pwd)"
SUITE="$(cd "$MC/.." && pwd)"
REPO="$(cd "$SUITE/../.." && pwd)"
ADO="$REPO/wooldidinstall/wooldid.ado"
# 2026-09-13 (audit F6): ADVB_OUT=<dir> re-runs into a fresh tree. mc_advb.do REPLACES
# <outdir>/<cfg>_s<n>.csv and the analyzer writes <outdir>/advb_report.md, so a re-run into
# the default out/ overwrites the 2026-08-26 record it should be compared against.
OUT="${ADVB_OUT:-$MC/out}"
GRID="$MC/advb_grid.json"
STATA_BIN="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
JOBS="${STATA_JOBS:-8}"
mkdir -p "$OUT/logs"
LOG="$OUT/advb_console.log"
note() { echo "== $*" | tee -a "$LOG"; }

pool_wait() { while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done; }

# emit one "cfg|shards|reps|base|G|U|T|vce|err|mode|extra" line per requested config
plan() {
  python3 - "$GRID" "$@" <<'PY'
import json, sys
grid = json.load(open(sys.argv[1]))
want = set(sys.argv[2:])
for c in grid["configs"]:
    if want and c["config"] not in want:
        continue
    print("|".join(str(x) for x in [
        c["config"], c["shards"], c["reps_per_shard"], c["base_seed"],
        c["G"], c["U"], c["T"], c["vce"], c["err"], c["mode"],
        c.get("stag", 0), c.get("scs", 1), c.get("pwr", 0),
        c.get("extra", "")]))
PY
}

FAIL=0
ARGS=("$@")
ANALYZE_ONLY=0
if [ "${1:-}" = "analyze" ]; then ANALYZE_ONLY=1; ARGS=(); fi

if [ "$ANALYZE_ONLY" -eq 0 ]; then
  note "AB campaign start (pool=$JOBS, ado=$ADO)"
  [ -f "$GRID" ] || { note "missing $GRID — the preregistered grid must exist first"; exit 2; }
  while IFS='|' read -r cfg ns reps base G U T vce err mode stag scs pwr extra; do
    [ -z "$cfg" ] && continue
    note "$cfg: G=$G U=$U T=$T vce=$vce err=$err mode=$mode stag=$stag scs=$scs pwr=$pwr extra='$extra' (${ns}x${reps} reps)"
    for s in $(seq 1 "$ns"); do
      pool_wait
      (
        d="$OUT/logs/${cfg}_s${s}"; mkdir -p "$d"
        # -b (not -e): batch log is named from the do-file/args rather than requiring a
        # non-numeric final token; a trailing numeric arg makes -e exit silently
        cd "$d" && "$STATA_BIN" -b do "$MC/mc_advb.do" "$ADO" "$cfg" "$s" "$reps" \
          "$base" "$OUT" "$G" "$U" "$T" "$vce" "$err" "$mode" \
          "$stag" "$scs" "$pwr" $extra >/dev/null 2>&1
      ) &
    done
  done < <(plan "${ARGS[@]+"${ARGS[@]}"}")
  wait
  note "all shards finished; checking completeness"
  while IFS='|' read -r cfg ns reps base G U T vce err mode stag scs pwr extra; do
    [ -z "$cfg" ] && continue
    ok=1
    for s in $(seq 1 "$ns"); do
      f="$OUT/${cfg}_s${s}.csv"
      if [ ! -f "$f" ]; then note "MISSING shard $f"; ok=0; continue; fi
      n=$(($(wc -l < "$f") - 1))
      [ "$n" -eq "$reps" ] || { note "SHORT shard $f: $n of $reps rows"; ok=0; }
    done
    [ $ok -eq 1 ] && note "$cfg: complete" || FAIL=1
  done < <(plan "${ARGS[@]+"${ARGS[@]}"}")
fi

note "analyzing"
python3 "$MC/analyze_advb.py" "$OUT" "$GRID" | tee -a "$LOG" || FAIL=1
[ $FAIL -eq 0 ] && note "AB campaign complete" || note "AB campaign: PROBLEMS PRESENT (see above)"
exit $FAIL
