clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

cap prog drop WDLEVSHOW
program define WDLEVSHOW
  syntax , LVL(real 0)
  display "RENDER [" "`lvl'" "] len=" strlen("`lvl'")
  capture noisily set level `lvl'
  display "  setlevel rc=" _rc " c(level)=" c(level)
  set level 95
end
foreach L in 99.99 99.9 99.98 99.95 99.5 92.5 10.01 88.8 79.9 44.4 {
  display as result "-- `L' --"
  capture noisily WDLEVSHOW , lvl(`L')
  display "  progrc=" _rc
}

local dN 50
local dT 10
local donsets "2004 2006"
include "`W'/mkdata.do"
foreach L in 99.9 99.91 99.92 99.93 99.94 99.96 99.97 99.99 98.9 88.8 79.9 44.4 12.3 {
  capture noisily wooldid y i t ttre, att robust level(`L')
  display "WLEV `L' rc=" _rc
}
display "P11_DONE"
