#!/bin/bash
# wooldid comprehensive testing & assessment suite — orchestrator.
#
# Usage: ./run_full.sh [component ...]
#   components: data oldnew compare cert newfeat battery coverage degeneracy all   (default: all)
#   env: STATA_BIN (default: the direct StataMP binary path — the PATH symlink is broken)
#        STATA_JOBS (default 6; hard project cap is 8 concurrent Stata instances. Stata-MP
#          claims several cores per instance, so a pool above the machine's core count
#          slows every shard instead of finishing sooner - measured, the slowest bootstrap
#          shard went from 330s at pool 6 to over 600s at pool 10 on ten cores)
#        ONLY=<comma-list of oldnew components, e.g. core,trig> to subset the dual-build
#        run (subset runs judge ONLY those components' configs and skip stale checks)
#
# Cost calibration (measured on this machine, 2026-09-11): the full `all` run is ~60
# minutes at pool=6 (oldnew ~25, cert ~23 — test_refmargin alone is 22 — newfeat ~7,
# battery ~5-8) — run it at every campaign checkpoint. Subset components run in
# seconds-to-minutes and are the per-fix tools. See README.md.
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
run_stata_job() { # $1=workdir $2=dofile-abs $3=timeout-s $4=tag $5=sentinel-regex ('' = none) [$6.. = do-file arguments]
  local dir="$1" dofile="$2" tmo="$3" tag="$4" sentinel="$5"
  shift 5
  # Stata names the batch log after the LAST cmdline token with its extension replaced
  # by .log: the do-file when it takes no arguments, else the last argument (the newfeat
  # scripts and the data builders end in .../data, so their log is data.log).
  local logtok="$dofile"
  [ $# -gt 0 ] && logtok="${@: -1}"
  pool_wait_slot
  local marker="$OUT/.running.$tag"
  : > "$marker"
  (
    started=$(date +%s)
    timedout=0
    cd "$dir" || { echo "FAIL $tag (cd $dir)" >> "$OUT/pool.log"; rm -f "$marker"; exit 1; }
    ( exec "$STATA_BIN" -e do "$dofile" "$@" ) >/dev/null 2>&1 &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
      sleep 10
      if [ $(( $(date +%s) - started )) -gt "$tmo" ]; then
        kill "$pid" 2>/dev/null
        sleep 10
        kill -9 "$pid" 2>/dev/null
        echo "TIMEOUT $tag ($dofile after ${tmo}s)" >> "$OUT/pool.log"
        timedout=1
        break
      fi
    done
    wait "$pid" 2>/dev/null
    dur=$(( $(date +%s) - started ))
    log="$(basename "$logtok")"; log="${log%.*}.log"
    if [ "$timedout" -eq 1 ]; then
      :   # one problem line per job: TIMEOUT is already recorded, no INCOMPLETE on top
    elif [ -z "$sentinel" ] || grep -qE "$sentinel" "$log" 2>/dev/null; then
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
  # 25 = 13 mkdata_suite + 11 extfix/build_all + 1 mkdata_admin. The three 2026-09-13
  # realistic fixtures (dta_s_stateyear, dta_s_biennial, dta_s_contcell) are inside the
  # mkdata_suite 13, self-tests included, and continuous_untreated_zero_compact was dropped
  # from extfix/build_all.do as a cf-identical duplicate of continuous_control_compact, so a
  # regeneration from an empty data/ still produces all 25.
  # Raise this WITH every new fixture: a stale threshold makes do_data skip regeneration
  # and every config on a new dataset dies r(601) on BOTH builds, which reads as a
  # parity pass.
  if [ "$(ls "$SUITE/data"/*.dta 2>/dev/null | wc -l | tr -d ' ')" -ge 25 ]; then
    note "data: frozen datasets present — skipping regeneration"
    return 0
  fi
  note "data: building suite datasets into $SUITE/data"
  mkdir -p "$SUITE/data"
  local ok=0
  # Stata names the batch log after the LAST cmdline token, so all three builders write
  # data.log (their argument is $SUITE/data). The graders used to look for
  # mkdata_suite.log / build_all.log, which Stata never writes (harness audit H-3), so
  # every regeneration reported those two steps FAILED. Each step now runs in its own
  # subdirectory (all three logs survive for triage), through the pool so it gets a
  # timeout (H-7), and is graded on its column-1 sentinel in that subdirectory's
  # data.log. Stata's shell exit status is 0 even when the do-file errors: the grep is
  # the gate. Steps stay sequential (the builders share $SUITE/data).
  local step sub dofile sentinel label
  for step in "suite|$SUITE/mkdata_suite.do|MKDATA_SUITE_SELFTEST_PASS|mkdata_suite" \
              "extfix|$SUITE/extfix/build_all.do|WDBT_FIXTURES_BUILT|extfix" \
              "admin|$SUITE/mkdata_admin.do|MKDATA_ADMIN_SELFTEST_PASS|mkdata_admin"; do
    IFS='|' read -r sub dofile sentinel label <<< "$step"
    rm -rf "$OUT/datagen/$sub"; mkdir -p "$OUT/datagen/$sub"
    run_stata_job "$OUT/datagen/$sub" "$dofile" 1800 "data_$sub" "^$sentinel" "$SUITE/data"
    pool_drain
    if grep -q "^$sentinel" "$OUT/datagen/$sub/data.log" 2>/dev/null; then
      note "data: $label OK"
    else
      note "data: $label FAILED — see $OUT/datagen/$sub/data.log"; ok=1
    fi
  done
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
  # the comparer checks that the new-side stamps are uniform, not that they are HEAD: a
  # re-judge after new commits would otherwise print PASS for an older build (H-11).
  # -xF on purpose: BSD grep treats `$' as an anchor even mid-pattern.
  local headsha; headsha=$(git -C "$REPO" rev-parse --short HEAD)
  grep -qxF "prov[new] = new|$headsha" "$OUT/oldnew_compare.txt" \
    || note "compare: WARNING — new-side dumps are not from HEAD ($headsha); this re-judges an older build"
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
      # Tier budgets are a HANG detector, not a runtime bound, so they carry real headroom
      # over the measured worst case at the default pool: the slowest fast shard
      # (shard_boot_00) runs 330s and the slowest slow shard (shard_boot_slow_00) 981s, so
      # these sit at roughly 3.5x each. That margin is what lets STATA_JOBS rise above the
      # default without a shard tripping its budget: Stata-MP claims several cores per
      # instance, so raising the pool past the core count slows each shard rather than
      # speeding the set, and a 600s fast tier turned that into a TIMEOUT at pool 10 on a
      # 10-core machine. Raise these with the pool, or leave the pool at the default.
      tmo=1200; [ "$speed" = "slow" ] && tmo=3600; [ "$speed" = "vslow" ] && tmo=7200
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
  note "cert: working-tree certification suites (test/test_*.do + check_validate_sync.py; only oldnew stages HEAD)"
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
  note "newfeat: cross-check components (working tree; only oldnew stages HEAD)"
  : > "$OUT/newfeat.log"
  local floor
  floor=$(date +%s)
  for f in "$SUITE"/newfeat/*.do; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .do)
    rm -rf "$OUT/newfeat/$b"
    mkdir -p "$OUT/newfeat/$b"
    # through the pool so a hung script cannot hang the suite (H-7). The pool's own
    # sentinel check reads data.log (the last cmdline token); the grade below reads
    # every log in the directory, as before.
    run_stata_job "$OUT/newfeat/$b" "$f" 1800 "nf_$b" '^[A-Z0-9_]+_ALL_PASS' \
      "$REPO/wooldidinstall/wooldid.ado" "$REPO/test/data"
  done
  pool_drain
  for f in "$SUITE"/newfeat/*.do; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .do)
    pass=1
    for lg in "$OUT/newfeat/$b"/*.log; do
      [ -e "$lg" ] && suite_log_pass "$lg" "$floor" && pass=0
    done
    if [ $pass -eq 0 ]; then
      echo "PASS $b" >> "$OUT/newfeat.log"
    else
      echo "FAIL $b (see $OUT/newfeat/$b/)" >> "$OUT/newfeat.log"
    fi
  done
  local nfail
  nfail=$(grep -c '^FAIL' "$OUT/newfeat.log" 2>/dev/null || true)
  note "newfeat: $(grep -c '^PASS' "$OUT/newfeat.log" 2>/dev/null || echo 0) passing, $nfail failing"
  grep '^FAIL' "$OUT/newfeat.log" | tee -a "$REPORT" || true
  [ "$nfail" -eq 0 ]
}

do_battery() {
  note "battery: working-tree bit-exact regression vs committed baseline"
  rm -rf "$REPO/test/run_suitebatt"
  mkdir -p "$REPO/test/run_suitebatt"
  # through the pool for the timeout and the BATTERY_DONE sentinel (H-7; the battery has
  # no kill loop of its own and a hung cell hung the suite). pool.log is truncated only
  # by do_oldnew, so a battery-only run would otherwise read the previous run's lines:
  # the tag carries this invocation's PID (pool.log lines are `STATUS tag (...)`).
  # runbattery.do resolves ../.. and ../data from the cwd, which the pool preserves.
  local btag="battery_$$"
  run_stata_job "$REPO/test/run_suitebatt" "$REPO/test/runbattery.do" 3600 "$btag" '^BATTERY_DONE'
  pool_drain
  if grep -qE "^(TIMEOUT|INCOMPLETE|FAIL) $btag " "$OUT/pool.log"; then
    note "battery: run did not complete — see $OUT/pool.log"
    return 1
  fi
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
      --allow-zero update "${dclass[@]+"${dclass[@]}"}" \
      > "$OUT/coverage.txt" 2>&1; then
    note "coverage: OK"
    # every WARN class reaches SUITE_REPORT.md, not just the last two lines (rc-only, A-5)
    grep -E '^(COVERAGE_WARN|info:|COVERAGE_OK)' "$OUT/coverage.txt" | tee -a "$REPORT"
    return 0
  fi
  note "coverage: GAPS — see $OUT/coverage.txt"
  tail -4 "$OUT/coverage.txt" | tee -a "$REPORT"
  return 1
}

do_degeneracy() {
  note "degeneracy: hidden-coverage assertion (degenerate results on the surface a config targets)"
  # CONTRACT (harness audit F12 + critique A2): degeneracy_scan.py --expect asserts every
  # clean-rc wellbehaved/stress config's degeneracy-flag stems against
  # degeneracy_expected.csv (NEW stem = fail, STALE expectation = fail), and a config
  # whose dump is MISSING must be a DEGEN_FAIL, never a skipped row: after a shard
  # timeout run_component continues into this component, and DEGEN_CHECK_OK over absent
  # dumps would be vacuous. Belt on the harness side: refuse when the comparer reported a
  # missing dump.
  #   2026-09-14 (author ruling): this guard USED to also refuse on a TIMEOUT/INCOMPLETE/FAIL
  #   line in pool.log. That condition was wrong, not merely redundant. pool.log is
  #   APPEND-ONLY, so once a shard has timed out the line is there forever - including after
  #   the shard is re-run by hand to its own SHARD_DONE sentinel and every dump is recovered.
  #   The guard then refused a scan that the scanner itself passes at nodump=0, and the only
  #   way to close a run was to invoke degeneracy_scan.py directly, which is exactly the
  #   bypass a guard is supposed to make unnecessary. The comparer's MISSING report is the
  #   right condition on its own: it is recomputed from the dumps that are on disk NOW, so it
  #   clears when a repair actually lands and stays tripped when one does not. The stale-dump
  #   case the pool.log condition was guarding against is also ruled out by that re-judge.
  # Runs on run_new only
  # (v1.3 dumps trip era flags); --dumps-old only annotates the old side's rc/properties.
  [ -f "$SUITE/degeneracy_scan.py" ] || { note "degeneracy: FAIL — $SUITE/degeneracy_scan.py is not present"; return 1; }
  [ -f "$SUITE/degeneracy_expected.csv" ] || { note "degeneracy: FAIL — $SUITE/degeneracy_expected.csv is not present"; return 1; }
  [ -d "$OUT/run_new" ] || { note "degeneracy: FAIL — no new-side dumps in $OUT/run_new (run oldnew first)"; return 1; }
  if grep -qE '^\[[^]]+\] MISSING dump' "$OUT/oldnew_compare.txt" 2>/dev/null; then
    note "degeneracy: FAIL — the last oldnew comparison reported a missing dump (see $OUT/oldnew_compare.txt); re-run oldnew first"
    return 1
  fi
  local extcfg=() logs=() old=()
  [ -f "$SUITE/configs_ext.csv" ] && extcfg=(--configs "$SUITE/configs_ext.csv")
  [ -d "$OUT/wd_new" ] && logs=(--logs "$OUT/wd_new")        # shard logs -> console-caution digest
  [ -d "$OUT/run_old" ] && old=(--dumps-old "$OUT/run_old")   # old-side rc/properties annotation
  if python3 "$SUITE/degeneracy_scan.py" --ado "$REPO/wooldidinstall/wooldid.ado" \
      --configs "$SUITE/configs.csv" "${extcfg[@]+"${extcfg[@]}"}" \
      --data-classes "$SUITE/data_classes.csv" --ledger "$SUITE/ledger.csv" \
      --dumps "$OUT/run_new" "${old[@]+"${old[@]}"}" "${logs[@]+"${logs[@]}"}" \
      --out "$OUT/degeneracy" --expect "$SUITE/degeneracy_expected.csv" \
      > "$OUT/degeneracy.txt" 2>&1; then
    note "degeneracy: OK"
    grep -E '^(DEGEN_WARN|DEGEN_CHECK_OK)' "$OUT/degeneracy.txt" | tee -a "$REPORT"
    return 0
  fi
  note "degeneracy: FAIL — see $OUT/degeneracy.txt"
  tail -8 "$OUT/degeneracy.txt" | tee -a "$REPORT"
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
    data|oldnew|compare|cert|newfeat|battery|coverage|degeneracy) run_component "$c" ;;
    all) for cc in data oldnew cert newfeat battery coverage degeneracy; do run_component "$cc"; done ;;
    *) echo "unknown component: $c"; exit 2 ;;
  esac
done
if [ $FAILRC -eq 0 ]; then
  note "run complete — ALL COMPONENTS PASS — report: $REPORT"
else
  note "run complete — FAILURES PRESENT — report: $REPORT"
fi
exit $FAILRC
