#!/bin/bash
# wooldid comprehensive testing & assessment suite — orchestrator.
#
# Usage: ./run_full.sh [component ...]
#   components: data oldnew compare cert newfeat battery coverage all   (default: all)
#   env: STATA_BIN (default: the direct StataMP binary path — the PATH symlink is broken)
#        STATA_JOBS (default 6; hard project cap is 8 concurrent Stata instances)
#        ONLY=<comma-list of oldnew components, e.g. core,trig> to subset the dual-build
#        run (subset runs judge ONLY those components' configs and skip stale checks)
#
# Cost calibration (measured on this machine, 2026-08-18): the full `all` run is ~20-30
# minutes at pool=6 — cheap enough to run at every campaign checkpoint. Subset components
# run in seconds-to-minutes and are the per-fix tools. See README.md.
#
# EXIT STATUS: nonzero if ANY component failed. The cumulative record is out/SUITE_REPORT.md.
set -u
SUITE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SUITE/../.." && pwd)"
STATA_BIN="${STATA_BIN:-/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp}"
STATA_JOBS="${STATA_JOBS:-6}"
OUT="$SUITE/out"
REPORT="$OUT/SUITE_REPORT.md"
mkdir -p "$OUT"
FAILRC=0

note() { echo "== $*" | tee -a "$REPORT"; }

# ---------- job pool with per-shard timeout (macOS has no `timeout`) ----------
# Markers are created by the PARENT before spawning (no TOCTOU overshoot); stale
# markers from interrupted runs are swept at startup.
rm -f "$OUT"/.running.* 2>/dev/null
pool_wait_slot() {
  while [ "$(ls "$OUT"/.running.* 2>/dev/null | wc -l | tr -d ' ')" -ge "$STATA_JOBS" ]; do
    sleep 5
  done
}
run_stata_job() { # $1=workdir $2=dofile-abs $3=timeout-s $4=tag $5=sentinel-regex ('' = none)
  local dir="$1" dofile="$2" tmo="$3" tag="$4" sentinel="$5"
  pool_wait_slot
  local marker="$OUT/.running.$tag"
  : > "$marker"
  (
    started=$(date +%s)
    cd "$dir" || { echo "FAIL $tag (cd $dir)" >> "$OUT/pool.log"; rm -f "$marker"; exit 1; }
    ( exec "$STATA_BIN" -e do "$dofile" ) >/dev/null 2>&1 &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
      sleep 10
      if [ $(( $(date +%s) - started )) -gt "$tmo" ]; then
        kill "$pid" 2>/dev/null
        sleep 10
        kill -9 "$pid" 2>/dev/null
        echo "TIMEOUT $tag ($dofile after ${tmo}s)" >> "$OUT/pool.log"
        break
      fi
    done
    wait "$pid" 2>/dev/null
    dur=$(( $(date +%s) - started ))
    log="$(basename "${dofile%.do}").log"
    if [ -z "$sentinel" ] || grep -qE "$sentinel" "$log" 2>/dev/null; then
      echo "OK $tag (${dur}s)" >> "$OUT/pool.log"
    else
      echo "INCOMPLETE $tag (no /$sentinel/ in $log, ${dur}s)" >> "$OUT/pool.log"
    fi
    rm -f "$marker"
  ) &
}
pool_drain() { wait; }

# grade a Stata log produced by a self-checking suite: the result line must START the
# line (batch echoes of source lines begin '. ' or '> ', so column-1 anchoring rejects
# them) and no column-1 FAILS line may exist. The log must be newer than $2.
suite_log_pass() { # $1=logfile $2=epoch-floor
  [ -f "$1" ] || return 1
  local mt
  mt=$(stat -f %m "$1" 2>/dev/null || echo 0)
  [ "$mt" -ge "$2" ] || return 1
  grep -qE '^[A-Z0-9_]+_ALL_PASS' "$1" || return 1
  grep -qE '^[A-Z0-9_]+_FAILS=' "$1" && return 1
  return 0
}

# ---------- components ----------
do_data() {
  if [ "$(ls "$SUITE/data"/*.dta 2>/dev/null | wc -l | tr -d ' ')" -ge 23 ]; then
    note "data: frozen datasets present — skipping regeneration"
    return 0
  fi
  note "data: building suite datasets into $SUITE/data"
  mkdir -p "$SUITE/data" "$OUT/datagen"
  local ok=0
  ( cd "$OUT/datagen" && "$STATA_BIN" -e do "$SUITE/mkdata_suite.do" "$SUITE/data" ) \
    && grep -q "MKDATA_SUITE_SELFTEST_PASS" "$OUT/datagen/mkdata_suite.log" \
    && note "data: mkdata_suite OK" || { note "data: mkdata_suite FAILED"; ok=1; }
  ( cd "$OUT/datagen" && "$STATA_BIN" -e do "$SUITE/extfix/build_all.do" "$SUITE/data" ) \
    && grep -q "WDBT_FIXTURES_BUILT" "$OUT/datagen/build_all.log" \
    && note "data: extfix OK" || { note "data: extfix FAILED"; ok=1; }
  ( cd "$OUT/datagen" && "$STATA_BIN" -e do "$SUITE/mkdata_admin.do" "$SUITE/data" ) \
    && grep -q "MKDATA_ADMIN_SELFTEST_PASS" "$OUT/datagen/data.log" \
    && note "data: mkdata_admin OK" || { note "data: mkdata_admin FAILED"; ok=1; }
  return $ok
}

_compare_cmd() { # writes compare output to $1; extra args appended
  local outfile="$1"; shift
  local extcfg=()
  [ -f "$SUITE/configs_ext.csv" ] && extcfg=(--configs "$SUITE/configs_ext.csv")
  python3 "$SUITE/compare_oldnew.py" "$OUT/run_old" "$OUT/run_new" --rtol 1e-10 \
    --ledger "$SUITE/ledger.csv" --configs "$SUITE/configs.csv" \
    "${extcfg[@]+"${extcfg[@]}"}" "$@" > "$outfile" 2>&1
}

do_compare() {
  note "compare: re-judging existing dumps"
  local flags=()
  [ -n "${ONLY:-}" ] && flags=(--only "$ONLY" --no-stale-check)
  _compare_cmd "$OUT/oldnew_compare.txt" "${flags[@]+"${flags[@]}"}"
  local rc=$?
  tail -4 "$OUT/oldnew_compare.txt" | tee -a "$REPORT"
  [ $rc -eq 0 ] && note "compare: PASS" || note "compare: FAILURES — see $OUT/oldnew_compare.txt"
  return $rc
}

do_oldnew() {
  local only="${ONLY:-}"
  note "oldnew: staging builds"
  mkdir -p "$OUT/build_old" "$OUT/build_new"
  git -C "$REPO" show main:wooldidinstall/wooldid.ado > "$OUT/build_old/wooldid.ado"
  # stage the NEW side from git HEAD (never the working tree) so the prov sha is honest
  git -C "$REPO" show HEAD:wooldidinstall/wooldid.ado > "$OUT/build_new/wooldid.ado"
  if ! git -C "$REPO" diff --quiet HEAD -- wooldidinstall/wooldid.ado; then
    note "oldnew: WARNING — working-tree wooldid.ado differs from HEAD; testing HEAD"
  fi
  local sha_old sha_new
  sha_old=$(git -C "$REPO" rev-parse --short main)
  sha_new=$(git -C "$REPO" rev-parse --short HEAD)
  local extcfg=()
  [ -f "$SUITE/configs_ext.csv" ] && extcfg=(--configs "$SUITE/configs_ext.csv")
  local comps=()
  [ -n "$only" ] && comps=(--components "$only")
  for side in old new; do
    if [ -z "$only" ]; then
      rm -rf "$OUT/run_$side" "$OUT/shards_$side" "$OUT/wd_$side"  # no mixed-vintage dumps on full runs
    fi
    mkdir -p "$OUT/run_$side" "$OUT/shards_$side"
    local sha; sha=$([ $side = old ] && echo "$sha_old" || echo "$sha_new")
    python3 "$SUITE/mkrun.py" --configs "$SUITE/configs.csv" \
      "${extcfg[@]+"${extcfg[@]}"}" \
      --preps "$SUITE/preps.json" --ado "$OUT/build_$side/wooldid.ado" \
      --label "$side" --sha "$sha" --data "$SUITE/data" --data "$REPO/test/data" \
      --outdir "$OUT/shards_$side" --dumpdir "$OUT/run_$side" \
      "${comps[@]+"${comps[@]}"}" \
      || { note "oldnew: mkrun($side) FAILED"; return 1; }
  done
  note "oldnew: running shards (pool=$STATA_JOBS)"
  : > "$OUT/pool.log"
  for side in old new; do
    while IFS='|' read -r dofile n speed; do
      [ -z "$dofile" ] && continue
      tmo=600; [ "$speed" = "slow" ] && tmo=1800
      # SI-21: every shard gets its OWN, FRESH working directory. dumpres globs the whole
      # cwd for saveplots output and then erases what it listed, so shards sharing a cwd
      # could record each other's files AND delete them before the owner dumped - a race
      # that produced a false FAILURE on a clean build and could equally hide a real
      # saveplots regression. Dumps still land together in run_$side via --dumpdir. The
      # directory is recreated on every launch (not just full runs) so a stale svp* file
      # from a previous subset run can never be attributed to the next run's first config.
      wd="$OUT/wd_$side/${dofile%.do}"
      rm -rf "$wd"; mkdir -p "$wd"
      run_stata_job "$wd" "$OUT/shards_$side/$dofile" "$tmo" \
        "${side}_${dofile%.do}" "^SHARD_DONE"
    done < "$OUT/shards_$side/shards.manifest"
  done
  pool_drain
  local bad
  bad=$(grep -cE '^(TIMEOUT|INCOMPLETE|FAIL)' "$OUT/pool.log" || true)
  note "oldnew: shards $(grep -c '^OK' "$OUT/pool.log") OK, $bad problem(s)"
  [ "$bad" -gt 0 ] && grep -E '^(TIMEOUT|INCOMPLETE|FAIL)' "$OUT/pool.log" | tee -a "$REPORT"
  note "oldnew: comparing"
  local flags=()
  [ -n "$only" ] && flags=(--only "$only" --no-stale-check)
  _compare_cmd "$OUT/oldnew_compare.txt" "${flags[@]+"${flags[@]}"}"
  local rc=$?
  tail -4 "$OUT/oldnew_compare.txt" | tee -a "$REPORT"
  if [ $rc -eq 0 ] && [ "$bad" -eq 0 ]; then
    note "oldnew: PASS"
    return 0
  fi
  note "oldnew: FAILURES — see $OUT/oldnew_compare.txt and $OUT/pool.log"
  return 1
}

do_cert() {
  note "cert: HEAD certification suites (test/test_*.do + check_validate_sync.py)"
  : > "$OUT/cert.log"
  local floor
  floor=$(date +%s)
  for f in "$REPO"/test/test_*.do; do
    b=$(basename "$f" .do)
    rm -f "$REPO/test/$b.log"   # a stale log must never grade a failed launch as PASS
    run_stata_job "$REPO/test" "$f" 3600 "cert_$b" '^[A-Z0-9_]+_ALL_PASS'
  done
  pool_drain
  local ok=0
  for f in "$REPO"/test/test_*.do; do
    b=$(basename "$f" .do)
    if suite_log_pass "$REPO/test/$b.log" "$floor"; then
      echo "PASS $b" >> "$OUT/cert.log"
    else
      echo "FAIL $b (see $REPO/test/$b.log)" >> "$OUT/cert.log"
      ok=1
    fi
  done
  python3 "$REPO/test/check_validate_sync.py" >/dev/null 2>&1 \
    && echo "PASS check_validate_sync" >> "$OUT/cert.log" \
    || { echo "FAIL check_validate_sync" >> "$OUT/cert.log"; ok=1; }
  note "cert: $(grep -c '^PASS' "$OUT/cert.log")/$(wc -l < "$OUT/cert.log" | tr -d ' ') passing"
  grep '^FAIL' "$OUT/cert.log" | tee -a "$REPORT" || true
  return $ok
}

do_newfeat() {
  note "newfeat: cross-check components (HEAD only)"
  : > "$OUT/newfeat.log"
  local floor
  floor=$(date +%s)
  for f in "$SUITE"/newfeat/*.do; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .do)
    rm -rf "$OUT/newfeat/$b"
    mkdir -p "$OUT/newfeat/$b"
    pool_wait_slot
    marker="$OUT/.running.nf_$b"
    : > "$marker"
    (
      cd "$OUT/newfeat/$b" && \
        "$STATA_BIN" -e do "$f" "$REPO/wooldidinstall/wooldid.ado" "$REPO/test/data" \
        >/dev/null 2>&1
      # Stata names the batch log after the LAST cmdline token — check every log here
      pass=1
      for lg in ./*.log; do
        [ -e "$lg" ] && suite_log_pass "$lg" "$floor" && pass=0
      done
      if [ $pass -eq 0 ]; then
        echo "PASS $b" >> "$OUT/newfeat.log"
      else
        echo "FAIL $b (see $OUT/newfeat/$b/)" >> "$OUT/newfeat.log"
      fi
      rm -f "$marker"
    ) &
  done
  pool_drain
  local nfail
  nfail=$(grep -c '^FAIL' "$OUT/newfeat.log" 2>/dev/null || true)
  note "newfeat: $(grep -c '^PASS' "$OUT/newfeat.log" 2>/dev/null || echo 0) passing, $nfail failing"
  grep '^FAIL' "$OUT/newfeat.log" | tee -a "$REPORT" || true
  [ "$nfail" -eq 0 ]
}

do_battery() {
  note "battery: HEAD bit-exact regression vs committed baseline"
  rm -rf "$REPO/test/run_suitebatt"
  mkdir -p "$REPO/test/run_suitebatt"
  ( cd "$REPO/test/run_suitebatt" && "$STATA_BIN" -e do ../runbattery.do ) >/dev/null 2>&1
  if python3 "$REPO/test/compare.py" "$REPO/test/baseline" "$REPO/test/run_suitebatt" --rtol 0 \
       > "$OUT/battery_compare.txt" 2>&1; then
    note "battery: IDENTICAL vs baseline"
    return 0
  fi
  note "battery: DIFFERENCES vs baseline — see $OUT/battery_compare.txt"
  return 1
}

do_coverage() {
  note "coverage: option-surface assertion"
  local extcfg=()
  [ -f "$SUITE/configs_ext.csv" ] && extcfg=(--configs "$SUITE/configs_ext.csv")
  local dclass=()
  [ -f "$SUITE/data_classes.csv" ] && dclass=(--data-classes "$SUITE/data_classes.csv")
  if python3 "$SUITE/coverage_check.py" --ado "$REPO/wooldidinstall/wooldid.ado" \
      --configs "$SUITE/configs.csv" "${extcfg[@]+"${extcfg[@]}"}" \
      --allow-zero update,permprobs "${dclass[@]+"${dclass[@]}"}" \
      > "$OUT/coverage.txt" 2>&1; then
    note "coverage: OK"
    tail -2 "$OUT/coverage.txt" | tee -a "$REPORT"
    return 0
  fi
  note "coverage: GAPS — see $OUT/coverage.txt"
  tail -4 "$OUT/coverage.txt" | tee -a "$REPORT"
  return 1
}

# ---------- main ----------
[ $# -eq 0 ] && set -- all
echo "# wooldid suite run $(date '+%Y-%m-%d %H:%M') (components: $*)" >> "$REPORT"
run_component() {
  "do_$1" || { FAILRC=1; note "COMPONENT FAILED: $1"; }
}
for c in "$@"; do
  case "$c" in
    data|oldnew|compare|cert|newfeat|battery|coverage) run_component "$c" ;;
    all) for cc in data oldnew cert newfeat battery coverage; do run_component "$cc"; done ;;
    *) echo "unknown component: $c"; exit 2 ;;
  esac
done
if [ $FAILRC -eq 0 ]; then
  note "run complete — ALL COMPONENTS PASS — report: $REPORT"
else
  note "run complete — FAILURES PRESENT — report: $REPORT"
fi
exit $FAILRC
