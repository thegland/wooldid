* p04_offchanges.do  -- three questions in one run:
*  (4) does safety(off) move any NUMBER on a run that would have succeeded anyway?
*  (5) no never-treated control + safety(off): what number does the bypass unlock?
*  (6) the esfixedbaseperiod pre/post r459 (ado:908) bypassed by safety(off)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

local gN 80
local gT 10
local gt0 2000
local gonsets 2004 2006
local gdose unitconst
local gysd 0.30
include "`GEN'/gen_core.do"
save "p04_noisy.dta", replace

********** (4) numbers-unchanged battery **********
local c1  att cluster(unit)
local c2  att att_i att_it att_itime cluster(region) esprelength(2) espostlength(2)
local c3  att subgroups(sg2) cluster(unit) esprelength(1)
local c4  att controls(x1 x2) fe(region) cluster(unit) esprelength(2) espostlength(2) jointtests
local c5  att contreat(dose) contreatcontrols(basic) lattice(all 3) cluster(region)
local c6  att poisson cluster(unit)
local c7  att att_i cluster(unit) impvalues
local c8  att cluster(unit) esprelength(2) espostlength(2) makeplots adversarial esfixedbaseperiod

forval k = 1/8 {
  use "p04_noisy.dta", clear
  local ovar ync
  if `k' == 6 local ovar ycount
  cap drop _all_marker
  capture noisily wooldid `ovar' unit t ttre, `c`k''
  local rcA = _rc
  cap mat drop bA VA MA
  cap mat bA = e(b)
  cap mat VA = e(V)
  cap mat MA = e(wdidmainresults)
  capture noisily wooldid `ovar' unit t ttre, `c`k'' safety(off)
  local rcB = _rc
  cap mat drop bB VB MB
  cap mat bB = e(b)
  cap mat VB = e(V)
  cap mat MB = e(wdidmainresults)
  local db = 99
  local dv = 99
  local dm = 99
  cap local db = mreldif(bA, bB)
  cap local dv = mreldif(VA, VB)
  cap local dm = mreldif(MA, MB)
  dis "CFG`k'_RC_A=`rcA' RC_B=`rcB' MRELDIF_b=`db' MRELDIF_V=`dv' MRELDIF_M=`dm'"
}
dis "P04_BATTERY_DONE"

********** (5) no never-treated, safety(off) **********
use "p04_noisy.dta", clear
* exact DGP version so the truth is exact
clear
local gN 80
local gT 10
local gt0 2000
local gonsets 2003 2005
local gysd 0
local gaisd 0
include "`GEN'/gen_core.do"
qui replace ttre = 2007 if ttre >= .     // EVERY unit treated inside the window
qui replace post = (t >= ttre)
qui replace h = t - ttre
qui replace truetau = 0
qui replace truetau = 1.5 if post == 1 & ttre == 2003
qui replace truetau = 2.0 if post == 1 & ttre == 2005
qui replace truetau = 2.5 if post == 1 & ttre == 2007
qui replace ync = a_i + _gt + truetau
count if ttre >= .
dis "NEVER_TREATED_ROWS=" r(N)
save "p04_nonever.dta", replace

dis "=== NN default safety ==="
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "NN_DEF_RC=" _rc
dis "=== NN safety(off) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "NN_OFF_RC=" _rc
cap drop insmp
gen byte insmp = e(sample)
mat NN = e(wdidmainresults)
mat list NN, format(%21x)
dis "NN_N=" e(N) " NN_DROPPED=" e(N_dropped)
local TRUTHSMP insmp
local TRUTHI   unit
include "`GEN'/gen_truth.do"
dis "NN_TRUE_att=" `TRUE_att'
dis "NN_TRUE_att_hex=" %21x `TRUE_att'
dis "NN_END"

********** (6) esfixedbaseperiod pre/post r459 at ado:908 **********
use "p04_nonever.dta", clear
* every treated unit's onset is its FIRST observed period -> no pre-period at all
qui replace ttre = 2000
qui replace post = 1
qui replace truetau = 2
qui replace ync = a_i + _gt + truetau
dis "=== ESF default safety ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(1) espostlength(2)
dis "ESF_DEF_RC=" _rc
dis "=== ESF safety(off): the ado:908 refusal removed ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(1) espostlength(2) safety(off)
dis "ESF_OFF_RC=" _rc
cap mat list e(wdidmainresults), format(%21x)
dis "P04_DONE"
