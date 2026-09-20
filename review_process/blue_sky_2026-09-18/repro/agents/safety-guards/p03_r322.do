* p03_r322.do -- the r(322) dropped-treated-cell guard (wooldid_ccscreen, ado:1713).
* Realistic trip: region-by-year fixed effects where one region holds exactly one
* treated state, so that state's cell dummies are spanned by the region#t effects.
* QUESTION: what does safety(off) unlock -- a refusal removed, or a wrong ATT?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

local gN 30
local gT 10
local gt0 2000
local gonsets 2004 2006
local gysd 0
local gaisd 0
local gtrendsd 0
local gnever 0.40
include "`GEN'/gen_core.do"

* -------- the realistic pathology: region 6 contains exactly ONE treated state -----
qui sum unit if ttre < .
local lonestate = r(min)
dis "LONESTATE=`lonestate'"
gen long regx = region
qui replace regx = 6 if unit == `lonestate'
egen long regt  = group(regx t)      // region-by-year FE
egen long regt5 = group(region t)    // the clean control: every region has >1 state
save "p03.dta", replace

dis "=== CONTROL: clean region#t FEs, default safety. must equal the design truth ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) fe(regt5)
dis "K_RC=" _rc
scalar K_att = e(wdid_att_b)
gen byte insmpK = e(sample)
local TRUTHSMP insmpK
local TRUTHI   unit
include "`GEN'/gen_truth.do"
dis "K_TRUE=" %21x `TRUE_att'
mat KM = e(wdidmainresults)
mat list KM, format(%21x)
dis "K_END"

dis "=== TRIP: region 6 holds one treated state, default safety ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) fe(regt)
dis "T_RC=" _rc

dis "=== BYPASS: same run with safety(off) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) fe(regt) safety(off)
dis "S_RC=" _rc
mat SM = e(wdidmainresults)
mat list SM, format(%21x)
cap drop insmpS
gen byte insmpS = e(sample)
cap mat GT = e(gteffects)
dis "S_N=" e(N)
local TRUTHSMP insmpS
local TRUTHI   unit
include "`GEN'/gen_truth.do"
dis "S_TRUE_att=" %21x `TRUE_att'
dis "S_TRUE_att_printed=" `TRUE_att'
dis "S_END"

dis "=== how much weight does the unidentified state carry? ==="
qui count if insmpS & ttre<. & t>=ttre
local npost = r(N)
qui count if insmpS & ttre<. & t>=ttre & unit==`lonestate'
dis "POSTROWS_TOTAL=`npost'  POSTROWS_LONESTATE=" r(N)
qui sum truetau if insmpS & ttre<. & t>=ttre
dis "MEAN_TRUETAU_ALLPOST=" r(mean)
qui sum truetau if insmpS & ttre<. & t>=ttre & unit!=`lonestate'
dis "MEAN_TRUETAU_EXCL_LONE=" r(mean)

dis "=== does e(gteffects) carry structural zeros for the lone state? ==="
cap mat list GT
dis "P03_DONE"
