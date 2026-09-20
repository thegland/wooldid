clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

* Does Stata itself accept set level 99.99 ?
cap noi set level 99.99
dis "SETLEVEL_9999_RC=" _rc "  c(level)=" c(level)
set level 95

* what does syntax Level(real) leave in the local?
cap prog drop SHOWLEV
prog def SHOWLEV
  syntax , Level(real 0)
  dis "LOCAL_AS_STRING=[`level']"
  dis "HEX=" %21x `level'
  cap noi set level `level'
  dis "SETLEV_RC=" _rc "  c(level)=" c(level)
  set level 95
end
foreach L in 99.99 99.98 99.9 97.5 33.33 66.67 10 10.01 {
  dis as result "---- syntax-real round trip for `L' ----"
  cap noi SHOWLEV , level(`L')
}

local dN 60
local dT 10
local donsets "2004 2006"
include "`W'/mkdata.do"
foreach L in 99.99 99.98 99.9 10.01 {
  cap noi wooldid y i t ttre, att robust level(`L')
  dis "WOOLDID_`L'_RC=" _rc "  e(level)=" e(level)
}
* and with the engines, whose own range check names 99.99 as legal
cap noi wooldid y i t ttre, att cluster(i) bootstrap(standard) reps(20) rseed(1) dots(0) level(99.99)
dis "BOOT_9999_RC=" _rc
dis "P9_DONE"
