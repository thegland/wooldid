* probe11.do -- SELF-FOUND differential: the CELL-LEVEL surfaces (e(gteffects)) and the
* SUBGROUP surfaces under the two bootstrap legs.  Under cluster() the two legs draw the
* SAME resamples (probe1/cmp1), so every cell number they report must agree.  The reid leg
* harvests cells through a NEW-i -> OLD-i remap that keeps only the FIRST matching row
* (wooldid_bootcells, ado mata); the count leg matches cell keys directly.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 1111
local gN 45
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
save base11.dta, replace

local R 30
local SD 6060

cap prog drop dumpgt
prog def dumpgt
  syntax , tag(string)
  cap mat GT = e(gteffects)
  if _rc {
    dis "`tag'|NO_GTEFFECTS"
    exit
  }
  dis "`tag'|gt rows=" rowsof(GT) " cols=" colsof(GT)
  dis "`tag'|gt colnames=" "`: colfullnames GT'"
  mat list GT, format(%21x)
end

cap prog drop dumpse
prog def dumpse
  syntax , tag(string)
  tempname b v
  mat `b' = e(b)
  mat `v' = e(V)
  local K = colsof(`b')
  forvalues k = 1/`K' {
    local nm : word `k' of `: colfullnames `b''
    dis "`tag'|" "`nm'" "|b=" %21x `b'[1,`k'] "|se=" %21x sqrt(`v'[`k',`k'])
  }
end

dis "=== GT1: att gteffects cluster(unit)  (COUNT leg) ==="
use base11.dta, clear
cap noi wooldid ync cohort t ttre, att gteffects cluster(cohort) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "GT1_RC=" _rc
dumpgt , tag(GT1_COUNT)
dumpse , tag(GT1_COUNT)

dis "=== GT2: att att_i gteffects cluster(cohort) (REID leg) ==="
use base11.dta, clear
cap noi wooldid ync cohort t ttre, att att_i gteffects cluster(cohort) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "GT2_RC=" _rc
dumpgt , tag(GT2_REID)
dumpse , tag(GT2_REID)

dis "=== SG1: att subgroups(sg2) cluster(unit) (COUNT leg) ==="
use base11.dta, clear
cap noi wooldid ync unit t ttre, att subgroups(sg2) cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "SG1_RC=" _rc
dumpse , tag(SG1_COUNT)
cap mat PW = e(pwsgs)
if _rc == 0 mat list PW, format(%21x)

dis "=== SG2: att att_i subgroups(sg2) cluster(unit) (REID leg) ==="
use base11.dta, clear
cap noi wooldid ync unit t ttre, att att_i subgroups(sg2) cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "SG2_RC=" _rc
dumpse , tag(SG2_REID)
cap mat PW = e(pwsgs)
if _rc == 0 mat list PW, format(%21x)

dis "=== MW: multi-way cluster with bootstrap must be refused ==="
use base11.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit region) bootstrap(standard) reps(10) rseed(1) dots(0)
dis "MW_RC=" _rc

dis "PROBE11_DONE"
