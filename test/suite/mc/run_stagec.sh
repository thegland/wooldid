#!/bin/bash
# Shard launcher for the three 2026-09-13 Stage C scaffolds. They share one positional
# contract, so one runner covers all three rather than three near-identical scripts:
#
#   <ado> <shard> <nreps> <seedbase> <outdir> <B> <arg7> <arg8> [<arg9> <arg10> <arg11>]
#
#   driver          arg7                 arg8 (the phase/arm/geometry list, + or , separated)
#   mc_bands.do     TAU (planted effect) ARMS   an+bs+bt
#   mc_ri_ext.do    TAU                  PHASES E1+E2+E3+E4+E5
#   mc_adv_bs12.do  SLOPE (0 = null)     GEOMS  AB1+AB13+AB12 (default); AB1C is the
#                                                continuous-treatment rung, non-default, an-arm only
#
# B is the replication count the driver hands to bootstrap()/permute(); wooldid's floor for
# permreps() is 20, so B < 20 makes the permutation arms refuse rather than run.
# Each shard writes <outdir>/<driver stem>_s<shard>.csv, one row per (rep, arm).
#
# Usage:
#   ./run_stagec.sh bands  <outdir> <nreps> <B> <tau>   <arms|->   <shard>...
#   ./run_stagec.sh riext  <outdir> <nreps> <B> <tau>   <phases|-> <shard>...
#   ./run_stagec.sh bs12   <outdir> <nreps> <B> <slope> <geoms|->  <shard>...
# "-" takes the driver's own default for that list. Example (a 200-rep null size ladder):
#   ./run_stagec.sh bs12 out3/bs12 200 999 0 AB1+AB13+AB12 1 2 3 4
#
# mc_bands.do alone accepts three FURTHER positional arguments (arg9..11 = G U T, the cohort
# count, units per cohort and period count; default 30 6 12 = the AB1 geometry). This runner
# always takes those defaults. To vary the geometry, call Stata directly:
#   "$STATA" -e do mc_bands.do <ado> <shard> <nreps> <seedbase> <outdir> <B> <tau> <arms> 40 8 16
#
# Like run_mcri2.sh this honours the machine's Stata instance cap by counting processes whose
# ARGV[0] is the Stata binary (a pgrep on " -e do " also matches monitoring pipelines and
# deadlocks on a busy box), and runs at most 2 of its own shards at a time.
set -u
STATA=/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp
REPO="$(cd "$(dirname "$0")/../../.." && pwd)"
ADO=$REPO/wooldidinstall/wooldid.ado
MC="$(cd "$(dirname "$0")" && pwd)"

[ $# -ge 7 ] || { sed -n '2,26p' "$0"; exit 2; }
case "$1" in
  bands) DO=$MC/mc_bands.do ;;
  riext) DO=$MC/mc_ri_ext.do ;;
  bs12)  DO=$MC/mc_adv_bs12.do ;;
  *) echo "run_stagec.sh: unknown driver '$1' (bands|riext|bs12)" >&2; exit 2 ;;
esac
[ -f "$DO" ] || { echo "run_stagec.sh: $DO is missing" >&2; exit 2; }
[ -f "$ADO" ] || { echo "run_stagec.sh: $ADO is missing" >&2; exit 2; }
shift
OUT=$1; NREPS=$2; B=$3; A7=$4; A8=$5; shift 5
SHARDS=("$@")
[ "$A8" = "-" ] && A8=""
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

nstata () { ps -eo args= | grep -c "^${STATA} -e do" || true; }
waitslot () {
  local n
  while :; do
    n=$(nstata)
    [ "$n" -lt 6 ] && break
    echo "$(date +%H:%M:%S) cap-wait: $n real stata -e do jobs running"
    sleep 60
  done
}

for S in "${SHARDS[@]}"; do
  waitslot
  echo "$(date +%H:%M:%S) launching shard $S (stata jobs before launch: $(nstata))"
  # seed base 900000 + shard keeps these campaigns off every earlier campaign's stream
  ( cd "$OUT" && "$STATA" -e do "$DO" "$ADO" "$S" "$NREPS" $((900000 + S)) "$OUT" "$B" "$A7" "$A8" ) &
  sleep 25
done
wait
echo "$(date +%H:%M:%S) all shards done; rows: $(cat "$OUT"/*_s*.csv 2>/dev/null | wc -l)"
