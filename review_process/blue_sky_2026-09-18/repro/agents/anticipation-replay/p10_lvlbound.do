clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

* ---- 1. does Stata itself take these? --------------------------------------
foreach L in 99.99 99.98 99.95 99.9 99.5 92.5 10.01 {
  cap noi set level `L'
  dis "STATA_SETLEVEL `L' rc=" _rc " c(level)=" c(level)
  set level 95
}

* ---- 2. what does `syntax , Level(real)' leave in the macro? ----------------
cap prog drop WDLEVSHOW
program define WDLEVSHOW
  syntax , Level(real 0)
  local s "`level'"
  display "MACRO_FOR `s' -> [" "`level'" "] strlen=" strlen("`level'")
  capture noisily set level `level'
  display "  set level rc=" _rc " c(level)=" c(level)
  set level 95
end
foreach L in 99.99 99.98 99.95 99.9 99.5 92.5 10.01 {
  capture noisily WDLEVSHOW , level(`L')
  display "  WDLEVSHOW rc=" _rc
}

* ---- 3. wooldid itself -----------------------------------------------------
local dN 60
local dT 10
local donsets "2004 2006"
include "`W'/mkdata.do"
foreach L in 99.99 99.98 99.95 99.9 99.5 92.5 10.01 {
  capture noisily wooldid y i t ttre, att robust level(`L')
  display "WOOLDID level(`L') rc=" _rc " e(level)=" e(level) " cols=" "`: colnames e(wdidmainresults)'"
}
display "P10_DONE"
