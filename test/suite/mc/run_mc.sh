#!/bin/bash
# run_mc.sh — MC statistical-calibration orchestrator (PREREG_NEWFEAT_2026-08-24).
# A DEEP-VALIDATION component, deliberately NOT part of run_full.sh all (cost: ~1.5h).
# usage: ./run_mc.sh [pretest] [ri] [adv] [pilot] [arms] [analyze]     or: ./run_mc.sh all
# STATA_JOBS caps the Stata pool (default 8 = the project hard cap). MC_OUT=<dir> sets the
# output tree (default out/ = the 2026-08-24 record).
set -u
MC="$(cd "$(dirname "$0")" && pwd)"
SUITE="$(cd "$MC/.." && pwd)"
REPO="$(cd "$SUITE/../.." && pwd)"
ADO="$REPO/wooldidinstall/wooldid.ado"
# 2026-09-13 (audit): MC_OUT=<dir> re-runs into a fresh tree. The drivers REPLACE
# <outdir>/<cfg>_s<n>.csv and analyze_mc.py writes <outdir>/mc_report.md, so a re-run into the
# default out/ overwrites the 2026-08-24 record (MC_REPORT_2026-08-24.md) it is compared against.
OUT="${MC_OUT:-$MC/out}"
STATA_BIN="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
JOBS="${STATA_JOBS:-8}"
mkdir -p "$OUT/logs"
note() { echo "== $*" | tee -a "$OUT/mc_console.log"; }

pool_wait() { while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 1; done; }

run_shards() { # driver cfg nshards reps base [gamma]
  local drv=$1 cfg=$2 ns=$3 rps=$4 base=$5; shift 5
  local s
  for s in $(seq 1 "$ns"); do
    pool_wait
    (
      d="$OUT/logs/${cfg}_s${s}"; mkdir -p "$d"
      cd "$d" && "$STATA_BIN" -e do "$MC/$drv" "$ADO" "$cfg" "$s" "$rps" "$base" "$OUT" "$@" >/dev/null 2>&1
    ) &
  done
}

check_shards() { # cfgprefix nshards expectedreps
  local cfg=$1 ns=$2 exp=$3 ok=1 s f n
  for s in $(seq 1 "$ns"); do
    f="$OUT/${cfg}_s${s}.csv"
    if [ ! -f "$f" ]; then note "MISSING shard: $f"; ok=0; continue; fi
    n=$(($(wc -l < "$f") - 1))
    [ "$n" -eq "$exp" ] || { note "SHORT shard $f: $n of $exp rows"; ok=0; }
  done
  [ $ok -eq 1 ] && note "$cfg: all $ns shards complete ($exp rows each)"
  return $((1 - ok))
}

FAIL=0

do_pretest() {
  note "pretest: null-fidelity (planted-line fidelity of the 9-period shape)"
  d="$OUT/logs/pretest"; mkdir -p "$d"
  ( cd "$d" && "$STATA_BIN" -e do "$MC/mc_pretest.do" "$ADO" >/dev/null 2>&1 )
  # Stata names the batch log after the LAST cmdline token (here the ado path)
  if grep -q '^NULLFID_ALL_PASS' "$d"/*.log 2>/dev/null; then
    note "pretest: PASS"
    grep -h '^NULLFID info' "$d"/*.log | tee -a "$OUT/mc_console.log"
  else
    note "pretest: FAIL — planted-line MC must not fire (see $d/mc_pretest.log)"
    FAIL=1
  fi
}

do_ri() {
  note "ri: R1-R7 (M=400 each, B=99; ~45 min at pool=$JOBS)"
  run_shards mc_ri.do R1 8 50 11000000
  run_shards mc_ri.do R2 8 50 12000000
  run_shards mc_ri.do R3 8 50 13000000
  run_shards mc_ri.do R4 8 50 14000000
  run_shards mc_ri.do R5 8 50 15000000
  run_shards mc_ri.do R6 8 50 16000000
  run_shards mc_ri.do R7 8 50 17000000
  wait
  for c in R1 R2 R3 R4 R5 R6 R7; do check_shards "$c" 8 50 || FAIL=1; done
}

do_adv() {
  note "adv: A1 A2 A3a A3b A6 (M=1000) + A5 (M=500; ~20 min at pool=$JOBS)"
  run_shards mc_adv.do A1 8 125 21000000
  run_shards mc_adv.do A2 8 125 22000000
  run_shards mc_adv.do A3a 8 125 23000000
  run_shards mc_adv.do A3b 8 125 24000000
  run_shards mc_adv.do A5 8 63 25000000
  run_shards mc_adv.do A6 8 125 26000000
  wait
  for c in A1 A2 A3a A3b A6; do check_shards "$c" 8 125 || FAIL=1; done
  check_shards A5 8 63 || FAIL=1
}

do_adv2() {
  note "adv2: A4 continuous adversarial under robust VCE (M=504; AMENDMENT A4)"
  run_shards mc_adv2.do A4 8 63 27000000
  wait
  check_shards A4 8 63 || FAIL=1
}

do_pilot() {
  note "pilot: ptp gamma=0 (104 reps)"
  run_shards mc_ptp.do pilot 8 13 31000000 0
  wait
  check_shards ptp_pilot 8 13 || FAIL=1
}

do_arms() {
  GW=$(python3 "$MC/pilot_median.py" "$OUT" mdt_wald) || { note "arms: pilot median (wald) unavailable"; FAIL=1; return; }
  GN=$(python3 "$MC/pilot_median.py" "$OUT" mdt_nis) || { note "arms: pilot median (nis) unavailable"; FAIL=1; return; }
  note "arms: W gamma=$GW, N gamma=$GN, Z gamma=0 (~10 min)"
  run_shards mc_ptp.do W 8 50 32000000 "$GW"
  run_shards mc_ptp.do N 8 50 33000000 "$GN"
  run_shards mc_ptp.do Z 8 38 34000000 0
  wait
  check_shards ptp_W 8 50 || FAIL=1
  check_shards ptp_N 8 50 || FAIL=1
  check_shards ptp_Z 8 38 || FAIL=1
}

do_analyze() {
  note "analyze: judging against the preregistered bands"
  if python3 "$MC/analyze_mc.py" "$OUT" | tee -a "$OUT/mc_console.log"; then
    note "analyze: complete — report at $OUT/mc_report.md"
  else
    note "analyze: PREREG FAIL verdicts present — see $OUT/mc_report.md"
    FAIL=1
  fi
}

[ $# -eq 0 ] && set -- all
for c in "$@"; do
  case "$c" in
    pretest) do_pretest ;;
    ri) do_ri ;;
    adv) do_adv ;;
    adv2) do_adv2 ;;
    pilot) do_pilot ;;
    arms) do_arms ;;
    analyze) do_analyze ;;
    all) do_pretest; [ $FAIL -eq 0 ] && { do_ri; do_adv; do_adv2; do_pilot; do_arms; do_analyze; } ;;
    *) echo "unknown component: $c"; exit 2 ;;
  esac
done
exit $FAIL
