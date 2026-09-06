#!/bin/bash
set -u
MC="$(cd "$(dirname "$0")" && pwd)"
STATA="/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
while : ; do
  n=$(pgrep -fl stata-mp | grep -c "^[0-9]* /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do " || true)
  [ "$n" -lt 6 ] && break
  sleep 60
done
cd "$MC/smoke2" && "$STATA" -e do verify2.do >/dev/null 2>&1
