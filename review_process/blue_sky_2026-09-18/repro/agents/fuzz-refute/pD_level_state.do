*=====================================================================
* pD: the level() re-render, mapped, plus the session-state corruption
*     it causes on the restore path (ado:68-72 / 574-578).
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- A. which two-decimal levels does the 16-digit re-render break? ----
cap prog drop lvlecho
prog def lvlecho
  syntax , [Level(real 0)]
  c_local echoed `level'
end
local nbreak 0
local breaklist
forvalues k = 1000/9999 {
  local hi = floor(`k'/100)
  local lo = substr("`=mod(`k',100)+100'",2,2)
  local L "`hi'.`lo'"
  lvlecho , level(`L')
  cap set level `echoed'
  if _rc != 0 {
    local nbreak = `nbreak' + 1
    if `nbreak' <= 80 local breaklist `breaklist' `L'
  }
}
set level 95
dis "D_NBREAK_OF_9000=" `nbreak'
dis "D_BREAKLIST_FIRST80= `breaklist'"
* one-decimal grid
local n1 0
local b1
forvalues k = 100/999 {
  local hi = floor(`k'/10)
  local lo = mod(`k',10)
  local L "`hi'.`lo'"
  lvlecho , level(`L')
  cap set level `echoed'
  if _rc != 0 {
    local n1 = `n1' + 1
    local b1 `b1' `L'
  }
}
set level 95
dis "D_NBREAK1_OF_900=" `n1'
dis "D_BREAK1LIST= `b1'"
* c(level) round trip: which SESSION levels break `local w = c(level)'?
local n2 0
local b2
foreach L in 90 95 99 99.5 99.9 99.99 99.95 97.5 90.25 95.4 68.27 84.13 92.5 99.75 {
  set level `L'
  local w = c(level)
  cap set level `w'
  if _rc != 0 {
    local n2 = `n2' + 1
    local b2 `b2' `L'
  }
  set level 95
}
dis "D_SESSIONBREAK_N=" `n2' "  LIST= `b2'"

* ---- B. session-state corruption ----
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz-refute/gen_mypanel.do"
local mN 50
local mseed 31337
include "`GEN'"
tempfile base
qui save "`base'", replace

* baseline: a good session level is restored
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 90
cap noi wooldid y i t ttre , att robust level(99)
dis "D_GOOD_RC=" _rc "  e(level)=" e(level) "  c(level)_after=" c(level) "  LEVEL0=[${WOOLDID_LEVEL0}]"

* hostile session level, explicit legal level()
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 99.9
dis "D_BEFORE_c(level)=" c(level)
cap noi wooldid y i t ttre , att robust level(95)
dis "D_HOSTILE_RC=" _rc
dis "D_HOSTILE_ELEVEL=" e(level) "  e(cmd)=" e(cmd) "  e(N)=" e(N)
dis "D_HOSTILE_cAFTER=" c(level)
dis "D_HOSTILE_LEVEL0=[${WOOLDID_LEVEL0}]"
dis "D_HOSTILE_LASTRC=[${WOOLDID_LASTRC}]"
* the results are posted but rc != 0 -- can the user replay them?
cap noi wooldid
dis "D_HOSTILE_REPLAY_RC=" _rc
* and is c(level) now silently different from what the user set?
dis "D_LEAK: user set 99.9, session now at " c(level)
set level 95
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
dis "=================== DONE pD ==================="
