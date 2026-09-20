clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid
local dN 50
local dT 10
local donsets "2004 2006"
include "`W'/mkdata.do"

* NO level() option -- the level comes from the session only
foreach L in 99.9 99.99 98.9 99.98 92.5 {
  set level `L'
  display as result "---- session set level `L' (no level() option) ----"
  capture noisily wooldid y i t ttre, att robust
  display "SESSION `L' rc=" _rc " e(level)=" e(level) " cols=" "`: colnames e(wdidmainresults)'"
  set level 95
}
* and the same under bootstrap(), where the engines have their own range check
set level 99.9
capture noisily wooldid y i t ttre, att cluster(i) bootstrap(standard) reps(20) rseed(1) dots(0)
display "SESSIONBOOT 99.9 rc=" _rc
set level 95
display "P12_DONE"
