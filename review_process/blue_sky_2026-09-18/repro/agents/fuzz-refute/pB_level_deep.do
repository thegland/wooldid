*=====================================================================
* pB: where does the level() failure really live, and how wide is it?
*     Also: FUZZ-06 redo with a legal permute() call; FUZZ-04 isolation.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- what does `local x = c(level)' render? ----
foreach L in 95 99.5 99.9 99.99 95.4 99.95 {
  set level `L'
  local w = c(level)
  dis "RENDER `L' -> c(level)=" c(level) "  local=[`w'] len=" strlen("`w'")
}
set level 95

* ---- is it wooldid, or is it reghdfe / regress / margins at c(level)=99.9? ----
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz-refute/gen_mypanel.do"
local mN 60
local mseed 808
include "`GEN'"
tempfile base
qui save "`base'", replace

set level 99.9
cap noi regress y x1
dis "B_REG_RC=" _rc
cap noi reghdfe y x1, absorb(i t)
dis "B_REGHDFE_RC=" _rc
qui regress y x1 i.t
cap noi margins, dydx(x1)
dis "B_MARGINS_RC=" _rc
cap noi wooldid y i t ttre , att robust
dis "B_WOOLDID_RC=" _rc
cap noi wooldid y i t ttre , att robust safety(off)
dis "B_WOOLDID_SAFOFF_RC=" _rc
cap noi wooldid y i t ttre , att unconditionalse
dis "B_WOOLDID_UNCOND_RC=" _rc
cap noi wooldid y i t ttre , att robust fastnose
dis "B_WOOLDID_FASTNOSE_RC=" _rc
set level 95
cap noi wooldid y i t ttre , att robust
dis "B_WOOLDID_AT95_RC=" _rc

* ---- map the whole two-decimal grid: which levels does syntax mangle? ----
cap prog drop lvlecho
prog def lvlecho
  syntax , [Level(real 0)]
  c_local echoed `level'
end
local nbad 0
local badlist
forvalues k = 1000/9999 {
  local L = `k'/100
  lvlecho , level(`L')
  if "`echoed'" != "`L'" {
    local nbad = `nbad' + 1
    if `nbad' <= 40 local badlist `badlist' `L'
  }
}
dis "B_NBAD_OF_9000=" `nbad'
dis "B_BADLIST_FIRST40= `badlist'"

* ---- FUZZ-06 redo: legal permute() call ----
use "`base'", clear
cap noi wooldid y i t ttre , att cluster(i) permute(reassign) permreps(24) rseed(9) dots(-5)
dis "B_PERM_NEG_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) permute(reassign) permreps(24) rseed(9) dots(0)
dis "B_PERM_ZERO_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) permute(reassign) permreps(24) rseed(9) dots(3)
dis "B_PERM_POS_RC=" _rc
cap noi wooldid y i t ttre , att cluster(i) bootstrap(standard) reps(24) rseed(9) dots(-5)
dis "B_BOOT_NEG_RC=" _rc

* ---- FUZZ-04 isolation: does a failed estimates replay poison the bare replay? ----
use "`base'", clear
cap noi wooldid y i t ttre , att robust
dis "B_F04_EST_RC=" _rc
cap noi wooldid
dis "B_F04_BARE_BEFORE_RC=" _rc
cap noi estimates store m1
cap noi wooldid
dis "B_F04_BARE_AFTERSTORE_RC=" _rc
cap noi estimates replay m1
dis "B_F04_REPLAY_RC=" _rc
cap noi wooldid
dis "B_F04_BARE_AFTERREPLAY_RC=" _rc
cap noi estimates restore m1
cap noi wooldid
dis "B_F04_BARE_AFTERRESTORE_RC=" _rc
dis "GLOBAL_LASTRC=[${WOOLDID_LASTRC}]"
* trailing blank claim
cap noi wooldid  
dis "B_F04_TRAILBLANK_RC=" _rc
dis "=================== DONE pB ==================="
