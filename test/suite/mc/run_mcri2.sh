#!/bin/zsh
# F058 PERMUTATION shard launcher. Honours the machine's Stata instance cap: it never
# launches while 6 or more REAL "stata-mp -e do" jobs are running, and runs at most 2 of
# my own shards at once.
#
# NOTE on the counter. The brief's literal probe
#     pgrep -fl stata-mp | grep -c " -e do "
# also matches other agents' shell pipelines (their own monitoring commands contain both
# strings), so on a busy box it reports 12-16 when only 3-4 Stata jobs exist and would
# deadlock forever. This script counts processes whose ARGV[0] is the Stata binary, which
# is what the cap means.
#
# Usage: run_mcri2.sh <outdir> <bD1> <bD2> <bD3> <bD4> <maxreps> <verify> <permreps> <shards...>
set -u
STATA=/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp
REPO=/Users/hegland/Documents/Research/wooldid_dev
ADO=$REPO/wooldidinstall/wooldid.ado
DO=$REPO/test/suite/mc/mc_ri2.do

OUT=$1; BD1=$2; BD2=$3; BD3=$4; BD4=$5; MAXR=$6; VER=$7; PR=$8; shift 8
SHARDS=("$@")
mkdir -p "$OUT"

nstata () {
  ps -eo args= | grep -c "^${STATA} -e do" || true
}

waitslot () {
  local n
  while :; do
    n=$(nstata)
    if [ "$n" -lt 6 ]; then break; fi
    echo "$(date +%H:%M:%S) cap-wait: $n real stata -e do jobs running"
    sleep 60
  done
}

for S in "${SHARDS[@]}"; do
  waitslot
  echo "$(date +%H:%M:%S) launching shard $S (real stata jobs before launch: $(nstata))"
  ( cd "$OUT" && "$STATA" -e do "$DO" "$ADO" "$S" "$OUT" "$BD1" "$BD2" "$BD3" "$BD4" 900000 "$MAXR" "$VER" "$PR" ) &
  sleep 25
done
wait
echo "$(date +%H:%M:%S) all shards done"
