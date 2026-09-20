clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- 1. map the two-decimal level grid with STRING-built levels ----
cap prog drop lvlecho
prog def lvlecho
  syntax , [Level(real 0)]
  c_local echoed `level'
end
local nbad 0
local badlist
forvalues k = 1000/9999 {
  local hi = floor(`k'/100)
  local lo = substr("`=mod(`k',100)+100'",2,2)
  local L "`hi'.`lo'"
  lvlecho , level(`L')
  if "`echoed'" != "`L'" {
    local nbad = `nbad' + 1
    if `nbad' <= 60 local badlist `badlist' `L'
  }
}
dis "C_NBAD_OF_9000=" `nbad'
dis "C_BADLIST_FIRST60= `badlist'"
* one-decimal grid too
local n1bad 0
local b1list
forvalues k = 100/999 {
  local hi = floor(`k'/10)
  local lo = mod(`k',10)
  local L "`hi'.`lo'"
  lvlecho , level(`L')
  if "`echoed'" != "`L'" {
    local n1bad = `n1bad' + 1
    local b1list `b1list' `L'
  }
}
dis "C_NBAD1_OF_900=" `n1bad'
dis "C_BAD1LIST= `b1list'"

* ---- 2. FUZZ-06 with a legal permute() ----
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz-refute/gen_mypanel.do"
local mN 60
local mseed 5150
include "`GEN'"
tempfile base
qui save "`base'", replace
cap noi wooldid y i t ttre , att cluster(i) permute(shuffle) permreps(24) rseed(9) dots(-5)
dis "C_PERM_NEG_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) permute(shuffle) permreps(24) rseed(9) dots(0)
dis "C_PERM_ZERO_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) permute(shuffle) permreps(24) rseed(9) dots(4)
dis "C_PERM_POS_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) bootstrap(standard) reps(24) rseed(9) dots(-5)
dis "C_BOOT_NEG_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) permute(shuffle) permreps(24) rseed(9) dots(-5.5)
dis "C_PERM_FRAC_RC=" _rc

* ---- 3. trailing blank from a CLEAN state ----
use "`base'", clear
macro drop WOOLDID_LASTRC
cap noi wooldid y i t ttre , att robust
dis "C_EST_RC=" _rc
dis "C_LASTRC_AFTER_EST=[${WOOLDID_LASTRC}]"
cap noi wooldid  
dis "C_TRAILBLANK_RC=" _rc
dis "C_LASTRC_AFTER_BLANK=[${WOOLDID_LASTRC}]"

* ---- 4. can level() override a hostile session level? ----
set level 99.9
cap noi wooldid y i t ttre , att robust level(95)
dis "C_OVERRIDE_RC=" _rc " e(level)=" e(level)
cap noi wooldid y i t ttre , att robust level(99.95)
dis "C_OVERRIDE2_RC=" _rc " e(level)=" e(level)
dis "C_CLEVEL_AFTER=" c(level)
set level 95
dis "=================== DONE pC ==================="
