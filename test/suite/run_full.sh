#!/bin/bash
# wooldid comprehensive testing & assessment suite — orchestrator.
#
# Usage: ./run_full.sh [component ...]
#   components: data oldnew cert newfeat battery coverage all   (default: all)
#   env: STATA_BIN (default: the direct StataMP binary path — the PATH symlink is broken)
#        STATA_JOBS (default 6; hard project cap is 8 concurrent Stata instances)
#        ONLY=<comma-list of oldnew components, e.g. core,trig> to subset the dual-build run
#
# The full run is HEFTY (hours) — meant for campaign-level assessment, not per-fix checks.
# Individual components (and ONLY= subsets) are the per-fix tools. See README.md.
set -u
SUITE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SUITE/../.." && pwd)"
STATA_BIN="${STATA_BIN:-/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp}"
STATA_JOBS="${STATA_JOBS:-6}"
OUT="$SUITE/out"
REPORT="$OUT/SUITE_REPORT.md"
mkdir -p "$OUT"

note() { echo "== $*" | tee -a "$REPORT"; }

# ---------- job pool with per-shard timeout (macOS has no `timeout`) ----------
pool_wait_slot() {
  while [ "$(ls "$OUT"/.running.* 2>/dev/null | wc -l | tr -d ' ')" -ge "$STATA_JOBS" ]; do
    sleep 10
  done
}
run_shard() { # $1=workdir $2=dofile-abs $3=timeout-seconds $4=tag
  local dir="$1" dofile="$2" tmo="$3" tag="$4"
  pool_wait_slot
  (
    marker="$OUT/.running.$tag"
    : > "$marker"
    cd "$dir" || { echo "FAIL cd $dir" >> "$OUT/pool.log"; rm -f "$marker"; exit 1; }
    ( exec "$STATA_BIN" -e do "$dofile" ) >/dev/null 2>&1 &
    pid=$!
    start=$(date +%s)
    while kill -0 "$pid" 2>/dev/null; do
      sleep 15
      if [ $(( $(date +%s) - start )) -gt "$tmo" ]; then
        kill "$pid" 2>/dev/null
        echo "TIMEOUT $tag ($dofile after ${tmo}s)" >> "$OUT/pool.log"
        break
      fi
    done
    wait "$pid" 2>/dev/null
    log="$(basename "${dofile%.do}").log"
    if grep -q "SHARD_DONE" "$log" 2>/dev/null; then
      echo "OK $tag" >> "$OUT/pool.log"
    else
      echo "INCOMPLETE $tag (no SHARD_DONE in $log)" >> "$OUT/pool.log"
    fi
    rm -f "$marker"
  ) &
}
pool_drain() { wait; }

# ---------- components ----------
do_data() {
  note "data: building suite datasets into $SUITE/data"
  mkdir -p "$SUITE/data" "$OUT/datagen"
  ( cd "$OUT/datagen" && "$STATA_BIN" -e do "$SUITE/mkdata_suite.do" "$SUITE/data" ) \
    && grep -q "MKDATA_SUITE_DONE" "$OUT/datagen/mkdata_suite.log" \
    && note "data: OK" || note "data: FAILED (see $OUT/datagen/mkdata_suite.log)"
  ( cd "$OUT/datagen" && "$STATA_BIN" -e do "$SUITE/extfix/build_all.do" "$SUITE/data" ) \
    && grep -q "WDBT_FIXTURES_BUILT" "$OUT/datagen/build_all.log" \
    && note "extfix data: OK" || note "extfix data: FAILED (see $OUT/datagen/build_all.log)"
}

do_oldnew() {
  local only="${ONLY:-}"
  note "oldnew: staging builds"
  mkdir -p "$OUT/build_old" "$OUT/build_new" "$OUT/run_old" "$OUT/run_new" \
           "$OUT/shards_old" "$OUT/shards_new"
  git -C "$REPO" show main:wooldidinstall/wooldid.ado > "$OUT/build_old/wooldid.ado"
  cp "$REPO/wooldidinstall/wooldid.ado" "$OUT/build_new/wooldid.ado"
  local sha_old sha_new
  sha_old=$(git -C "$REPO" rev-parse --short main)
  sha_new=$(git -C "$REPO" rev-parse --short HEAD)
  local extcfg=""
  [ -f "$SUITE/configs_ext.csv" ] && extcfg="--configs $SUITE/configs_ext.csv"
  local comps=""
  [ -n "$only" ] && comps="--components $only"
  for side in old new; do
    local sha; sha=$([ $side = old ] && echo "$sha_old" || echo "$sha_new")
    # shellcheck disable=SC2086
    python3 "$SUITE/mkrun.py" --configs "$SUITE/configs.csv" $extcfg \
      --preps "$SUITE/preps.json" --ado "$OUT/build_$side/wooldid.ado" \
      --label "$side" --sha "$sha" --data "$SUITE/data" --data "$REPO/test/data" \
      --outdir "$OUT/shards_$side" $comps || { note "oldnew: mkrun($side) FAILED"; return 1; }
  done
  note "oldnew: running shards (pool=$STATA_JOBS)"
  : > "$OUT/pool.log"
  for side in old new; do
    while IFS='|' read -r dofile n speed; do
      [ -z "$dofile" ] && continue
      tmo=2700; [ "$speed" = "slow" ] && tmo=5400
      run_shard "$OUT/run_$side" "$OUT/shards_$side/$dofile" "$tmo" "${side}_${dofile%.do}"
    done < "$OUT/shards_$side/shards.manifest"
  done
  pool_drain
  note "oldnew: pool results: $(sort "$OUT/pool.log" | uniq -c | tr '\n' ';' )"
  note "oldnew: comparing"
  local staleflag=""
  [ -n "$only" ] && staleflag="--no-stale-check"
  # shellcheck disable=SC2086
  python3 "$SUITE/compare_oldnew.py" "$OUT/run_old" "$OUT/run_new" --rtol 1e-10 \
    --ledger "$SUITE/ledger.csv" --configs "$SUITE/configs.csv" $staleflag \
    > "$OUT/oldnew_compare.txt" 2>&1
  local rc=$?
  tail -3 "$OUT/oldnew_compare.txt" | tee -a "$REPORT"
  [ $rc -eq 0 ] && note "oldnew: PASS" || note "oldnew: DIFFERENCES/STALE — see $OUT/oldnew_compare.txt"
}

do_cert() {
  note "cert: HEAD certification suites (test/test_*.do + check_validate_sync.py)"
  : > "$OUT/cert.log"
  for f in "$REPO"/test/test_*.do; do
    b=$(basename "$f" .do)
    run_shard "$REPO/test" "$f" 3600 "cert_$b"
  done
  pool_drain
  for f in "$REPO"/test/test_*.do; do
    b=$(basename "$f" .do)
    up=$(echo "${b#test_}" | tr '[:lower:]' '[:upper:]')
    if grep -q "_ALL_PASS" "$REPO/test/$b.log" 2>/dev/null; then
      echo "PASS $b" >> "$OUT/cert.log"
    else
      echo "FAIL $b (no *_ALL_PASS sentinel; grep $up in $REPO/test/$b.log)" >> "$OUT/cert.log"
    fi
  done
  python3 "$REPO/test/check_validate_sync.py" >/dev/null 2>&1 \
    && echo "PASS check_validate_sync" >> "$OUT/cert.log" \
    || echo "FAIL check_validate_sync" >> "$OUT/cert.log"
  note "cert: $(grep -c '^PASS' "$OUT/cert.log")/$(wc -l < "$OUT/cert.log" | tr -d ' ') passing"
  grep '^FAIL' "$OUT/cert.log" | tee -a "$REPORT" || true
}

do_newfeat() {
  note "newfeat: cross-check components (HEAD only)"
  mkdir -p "$OUT/newfeat"
  : > "$OUT/newfeat.log"
  for f in "$SUITE"/newfeat/*.do; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .do)
    mkdir -p "$OUT/newfeat/$b"
    (
      pool_wait_slot
      marker="$OUT/.running.nf_$b"; : > "$marker"
      cd "$OUT/newfeat/$b" && \
        "$STATA_BIN" -e do "$f" "$REPO/wooldidinstall/wooldid.ado" "$SUITE/data" "$REPO/test/data" \
        >/dev/null 2>&1
      if grep -q "_ALL_PASS" "$b.log" 2>/dev/null; then
        echo "PASS $b" >> "$OUT/newfeat.log"
      else
        echo "FAIL $b (see $OUT/newfeat/$b/$b.log)" >> "$OUT/newfeat.log"
      fi
      rm -f "$marker"
    ) &
  done
  pool_drain
  note "newfeat: $(grep -c '^PASS' "$OUT/newfeat.log" 2>/dev/null || echo 0) passing"
  grep '^FAIL' "$OUT/newfeat.log" | tee -a "$REPORT" || true
}

do_battery() {
  note "battery: HEAD bit-exact regression vs committed baseline"
  mkdir -p "$REPO/test/run_suitebatt"
  ( cd "$REPO/test/run_suitebatt" && "$STATA_BIN" -e do ../runbattery.do ) >/dev/null 2>&1
  if python3 "$REPO/test/compare.py" "$REPO/test/baseline" "$REPO/test/run_suitebatt" --rtol 0 \
       > "$OUT/battery_compare.txt" 2>&1; then
    note "battery: IDENTICAL vs baseline"
  else
    note "battery: DIFFERENCES vs baseline — see $OUT/battery_compare.txt"
  fi
}

do_coverage() {
  note "coverage: option-surface assertion"
  local extcfg=""
  [ -f "$SUITE/configs_ext.csv" ] && extcfg="--configs $SUITE/configs_ext.csv"
  # shellcheck disable=SC2086
  python3 "$SUITE/coverage_check.py" --ado "$REPO/wooldidinstall/wooldid.ado" \
    --configs "$SUITE/configs.csv" $extcfg \
    --allow-zero update,replace,verbose,cleanmatrices,oldsyntax,regtol \
    > "$OUT/coverage.txt" 2>&1 \
    && note "coverage: OK" || note "coverage: GAPS — see $OUT/coverage.txt"
  tail -2 "$OUT/coverage.txt" | tee -a "$REPORT"
}

# ---------- main ----------
[ $# -eq 0 ] && set -- all
echo "# wooldid suite run $(date '+%Y-%m-%d %H:%M')" >> "$REPORT"
for c in "$@"; do
  case "$c" in
    data) do_data ;;
    oldnew) do_oldnew ;;
    cert) do_cert ;;
    newfeat) do_newfeat ;;
    battery) do_battery ;;
    coverage) do_coverage ;;
    all) do_data; do_oldnew; do_cert; do_newfeat; do_battery; do_coverage ;;
    *) echo "unknown component: $c"; exit 2 ;;
  esac
done
note "run complete — report: $REPORT"
