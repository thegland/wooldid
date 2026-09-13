#!/bin/bash
set -u
MC="$(cd "$(dirname "$0")" && pwd)"
STATA="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
# smoke2/verify2.do and smoke2/smoke2.do are tracked as of 2026-09-13; the 2026-09-13 critic
# C2 inventory note (runner present, target untracked, so a fresh checkout had nothing to run)
# is closed by tracking them rather than by deleting this runner. smoke2/*.log stays ignored.
# Keep this guard: it makes a future re-break loud instead of silent.
[ -f "$MC/smoke2/verify2.do" ] || { echo "run_verify2.sh: $MC/smoke2/verify2.do is missing" >&2; exit 2; }
while : ; do
  n=$(pgrep -fl stata-mp | grep -c "^[0-9]* /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do " || true)
  [ "$n" -lt 6 ] && break
  sleep 60
done
cd "$MC/smoke2" && "$STATA" -e do verify2.do >/dev/null 2>&1
