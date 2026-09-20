* p06_unlocks.do -- characterise five more safety(off) unlocks: is the output a
* correct number, a wrong number, a crash, or a number with no sign it is unreliable?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

local gN 60
local gT 10
local gt0 2000
local gonsets 2004 2006
local gysd 0.30
include "`GEN'/gen_core.do"
save "p06.dta", replace

dis "=== (a) unconditionalse + att_i : refusal, then bypass ==="
capture noisily wooldid ync unit t ttre, att_i robust unconditionalse
dis "UA_RC=" _rc
capture noisily wooldid ync unit t ttre, att_i robust unconditionalse safety(off)
dis "UB_RC=" _rc
mat U1 = e(wdidmainresults)
mat list U1, format(%10.7f)
capture noisily wooldid ync unit t ttre, att_i robust
dis "UC_RC=" _rc
mat U2 = e(wdidmainresults)
mat list U2, format(%10.7f)
capture noisily wooldid ync unit t ttre, att_i cluster(unit)
dis "UD_RC=" _rc
mat U3 = e(wdidmainresults)
mat list U3, format(%10.7f)
dis "UNCOND_END"

dis "=== (b) fractional onset (mid-year effective date on an annual panel) ==="
use "p06.dta", clear
gen double ttrefrac = ttre + 0.5
capture noisily wooldid ync unit t ttrefrac, att cluster(unit) esprelength(2) espostlength(2)
dis "FA_RC=" _rc
capture noisily wooldid ync unit t ttrefrac, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "FB_RC=" _rc
cap mat list e(wdidmainresults), format(%10.7f)
cap dis "FB_N=" e(N)
* the honest recodes, for comparison
capture noisily wooldid ync unit t ttre, att cluster(unit) esprelength(2) espostlength(2)
dis "FC_RC=" _rc
cap mat list e(wdidmainresults), format(%10.7f)
gen double ttreup = ttre + 1
capture noisily wooldid ync unit t ttreup, att cluster(unit) esprelength(2) espostlength(2)
dis "FD_RC=" _rc
cap mat list e(wdidmainresults), format(%10.7f)
dis "FRAC_END"

dis "=== (c) negative weights ==="
use "p06.dta", clear
gen double wneg = 1
qui replace wneg = -0.5 in 7/9
capture noisily wooldid ync unit t ttre [aw=wneg], att cluster(unit)
dis "WA_RC=" _rc
capture noisily wooldid ync unit t ttre [aw=wneg], att cluster(unit) safety(off)
dis "WB_RC=" _rc
cap dis "WB_N=" e(N) " WB_DROPPED=" e(N_dropped)
* negative weight on rows the run would have discarded anyway (missing outcome)
use "p06.dta", clear
gen double wneg2 = 1
gen double ymiss = ync
qui replace ymiss = . in 7/9
qui replace wneg2 = -0.5 in 7/9
capture noisily wooldid ymiss unit t ttre [aw=wneg2], att cluster(unit)
dis "WC_RC=" _rc
capture noisily wooldid ymiss unit t ttre [aw=wneg2], att cluster(unit) safety(off)
dis "WD_RC=" _rc
cap dis "WD_N=" e(N)
dis "WEIGHT_END"

dis "=== (d) bootstrap(bayes) + cell-weighted : does safety(off) print its caution? ==="
use "p06.dta", clear
capture noisily wooldid ync unit t ttre, att_i cluster(unit) bootstrap(bayes) reps(20)
dis "BA_RC=" _rc
capture noisily wooldid ync unit t ttre, att_i cluster(unit) bootstrap(bayes) reps(20) safety(off)
dis "BB_RC=" _rc
dis "BAYES_END"

dis "=== (e) subgroups splitting a unit over time : does safety(off) print its caution? ==="
use "p06.dta", clear
capture noisily wooldid ync unit t ttre, att subgroups(sgvary) cluster(unit)
dis "SA_RC=" _rc
capture noisily wooldid ync unit t ttre, att subgroups(sgvary) cluster(unit) safety(off)
dis "SB_RC=" _rc
dis "P06_DONE"
