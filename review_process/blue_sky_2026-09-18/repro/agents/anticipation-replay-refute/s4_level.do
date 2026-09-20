clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop st
prog def st
  syntax anything(name=tag)
  dis "@@ `tag' rc=" _rc
end

use "dat.dta", clear

*==================== (a) level() option, session at 95 ====================
foreach L in 99.9 99.99 99.91 99.5 92.5 95 {
  set level 95
  macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
  cap noi wooldid y uid t ttre, att robust level(`L')
  local rc = _rc
  dis "@@A L=`L' rc=`rc' clevel_after=" c(level) " elevel=" e(level) " eN=" e(N) " ecmd=[`e(cmd)'] LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC] CLEAN=[$WOOLDID_CLEANUP]"
}

*==================== (b) NO level() option, session level varies ==========
foreach L in 99.9 99.91 99.99 32.22 95 99.5 {
  clear all
  use "dat.dta", clear
  macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
  set level `L'
  cap noi wooldid y uid t ttre, att robust
  local rc = _rc
  dis "@@B sess=`L' rc=`rc' clevel_after=" c(level) " elevel=" e(level) " eN=" e(N) " ecmd=[`e(cmd)'] LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"
  set level 95
}

*==================== (c) session 99.9 + level(95): worker SUCCEEDS ========
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 99.9
dis "@@C0 clevel_before=" c(level)
cap noi wooldid y uid t ttre, att robust level(95)
local rcC = _rc
dis "@@C rc=`rcC' clevel_after=" c(level) " elevel=" e(level) " eN=" e(N) " ecmd=[`e(cmd)'] LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC] CLEAN=[$WOOLDID_CLEANUP]"
mat lb = e(b)
dis "@@C b_att=" el(lb,1,1)

*=========== (c2) STICKY: the stale LEVEL0 survives; next calls =============
dis "@@C2 LEVEL0_still=[$WOOLDID_LEVEL0] clevel=" c(level)
cap noi wooldid y uid t ttre, att robust
dis "@@C2a no-level-call rc=" _rc " clevel=" c(level) " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"
cap noi wooldid y uid t ttre, att_i robust
dis "@@C2b second no-level-call rc=" _rc " clevel=" c(level) " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"
cap noi wooldid y uid t ttre, att robust level(90)
dis "@@C2c with-level-call rc=" _rc " clevel=" c(level) " elevel=" e(level) " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"
cap noi wooldid y uid t ttre, att robust
dis "@@C2d after-recovery no-level rc=" _rc " clevel=" c(level) " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"

*=========== (e) parse-stage failure while LEVEL0 is poisoned ===============
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 99.9
cap noi wooldid y uid t ttre, att robust level(95)
dis "@@E0 poisoned rc=" _rc " LEVEL0=[$WOOLDID_LEVEL0]"
cap noi wooldid y uid t ttre, att robust bogusoption
dis "@@E1 parsefail rc=" _rc " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"
cap noi wooldid nosuchvar uid t ttre, att robust
dis "@@E2 badvar rc=" _rc " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"

*=========== (d) staleness-guard gap: failed call then bare replay ==========
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
wooldid y uid t ttre, att robust
dis "@@D0 good rc=" _rc " eN=" e(N) " elevel=" e(level) " LASTRC=[$WOOLDID_LASTRC]"
set level 99.9
cap noi wooldid y uid t ttre, att_i robust
dis "@@D1 failed rc=" _rc " LASTRC=[$WOOLDID_LASTRC] LEVEL0=[$WOOLDID_LEVEL0] clevel=" c(level)
cap noi wooldid
dis "@@D2 bare-replay rc=" _rc " eN=" e(N) " elevel=" e(level) " ecmd=[`e(cmd)']"

*=========== (d-control) an ORDINARY failure does arm the guard =============
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
wooldid y uid t ttre, att robust
cap noi wooldid y uid t ttre, att robust bogusoption
dis "@@F1 ordinary-fail rc=" _rc " LASTRC=[$WOOLDID_LASTRC]"
cap noi wooldid
dis "@@F2 bare-replay rc=" _rc
set level 95
dis "S4_DONE"
