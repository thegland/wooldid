#!/bin/zsh
# Driver: runs the old-build and new-build reproductions in SEPARATE Stata processes (sequentially),
# respecting the shared cap (<8 stata-mp instances), then summarizes.
# Usage:  REPO=/path/to/wooldid_dev W=/path/to/scratch/contsg zsh repro_run.sh
: "${REPO:?set REPO=<wooldid_dev checkout root>}"
: "${W:?set W=<scratch workdir holding old/, new/ and repro_helpers.doh>}"
STATA="${STATA:-/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# repro_old.do/repro_new.do `include "$W/repro_helpers.doh"`; ship it into the workdir:
cp "$HERE/repro_helpers.doh" "$W/repro_helpers.doh" || exit 1
cd "$W" || exit 1
for leg in old new; do
  while [ "$(ps -eo args= | grep -c '^/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')" -ge 8 ]; do
    echo "cap reached; waiting 30s"; sleep 30
  done
  echo "launch $leg $(date +%H:%M:%S) (stata count before: $(ps -eo args= | grep -c '^/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp'))"
  "$STATA" -e do "$HERE/repro_$leg.do" "$REPO" "$W"
  echo "done $leg rc=$? $(date +%H:%M:%S)"
  grep -c "REPRO_ALL_${leg:u}_DONE" repro_$leg.log
done
WOOLDID_REPO="$REPO" WOOLDID_SCRATCH="$W" python3 "$HERE/repro_summarize.py" > repro_result.txt
cat repro_result.txt
