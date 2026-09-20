* p05_bypass_numbers.do -- what NUMBER does each safety(off) bypass unlock?
*  (A) no never-treated control group   (ado:955, r459)
*  (B) esfixedbaseperiod pre/post       (ado:908, r459)
*  (C) dropped treated cell dummy       (wooldid_ccscreen, ado:1713, r322)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

*=================== (A) no never-treated ===================
local gN 80
local gT 10
local gt0 2000
local gonsets 2003 2005
local gysd 0
local gaisd 0
local gtrendsd 0
include "`GEN'/gen_core.do"
qui replace ync = ync - truetau            // strip the old effect, keep a_i + g_t
qui replace ttre = 2007 if ttre >= .       // EVERY unit treated inside the window
qui replace post = (t >= ttre)
qui replace h = t - ttre
qui replace truetau = 0
qui replace truetau = 1.5 if post == 1 & ttre == 2003
qui replace truetau = 2.0 if post == 1 & ttre == 2005
qui replace truetau = 2.5 if post == 1 & ttre == 2007
qui replace ync = ync + truetau
count if ttre >= .
dis "NEVER_TREATED_ROWS=" r(N)
save "p05_nonever.dta", replace

dis "=== A1 default safety ==="
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "A1_RC=" _rc
dis "=== A2 safety(off) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "A2_RC=" _rc
cap drop insmp
gen byte insmp = e(sample)
mat NN = e(wdidmainresults)
mat list NN, format(%10.7f)
dis "A2_N=" e(N) " A2_DROPPED=" e(N_dropped)
local TRUTHSMP insmp
local TRUTHI   unit
include "`GEN'/gen_truth.do"
dis "A2_TRUE_att=" `TRUE_att'
dis "A2_END"

*=================== (B) esfixedbaseperiod, no pre-period ===================
use "p05_nonever.dta", clear
qui replace ync = ync - truetau
qui replace ttre = 2000                    // onset = first observed period for everyone
qui replace post = 1
qui replace truetau = 2
qui replace ync = ync + truetau
dis "=== B1 default safety ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(1) espostlength(2)
dis "B1_RC=" _rc
dis "=== B2 safety(off) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(1) espostlength(2) safety(off)
dis "B2_RC=" _rc
cap mat list e(wdidmainresults), format(%10.7f)
dis "B2_END"

*=================== (C) collinear co-timed policy control ===================
clear
local gN 40
local gT 10
local gt0 2000
local gonsets 2004 2006
local gysd 0
local gaisd 0
local gtrendsd 0
include "`GEN'/gen_core.do"
* a SECOND policy that one treated state adopted in exactly the same year
qui sum unit if ttre < .
local U = r(min)
dis "COLLIN_UNIT=`U'"
gen byte policy2 = (unit == `U' & t >= ttre & ttre < .)
save "p05_collin.dta", replace

dis "=== C0 no control: the correct answer ==="
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "C0_RC=" _rc
cap drop insmp0
gen byte insmp0 = e(sample)
mat C0 = e(wdidmainresults)
mat list C0, format(%10.7f)
local TRUTHSMP insmp0
local TRUTHI   unit
include "`GEN'/gen_truth.do"
dis "C0_TRUE_att=" `TRUE_att'

dis "=== C1 controls(policy2), default safety ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) controls(policy2)
dis "C1_RC=" _rc

dis "=== C2 controls(policy2), safety(off) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) controls(policy2) safety(off) gteffects
dis "C2_RC=" _rc
cap drop insmp2
gen byte insmp2 = e(sample)
mat C2 = e(wdidmainresults)
mat list C2, format(%10.7f)
dis "C2_N=" e(N)
local TRUTHSMP insmp2
local TRUTHI   unit
include "`GEN'/gen_truth.do"
dis "C2_TRUE_att=" `TRUE_att'
cap mat GT = e(gteffects)
cap mat list GT, format(%9.5f)
dis "P05_DONE"
