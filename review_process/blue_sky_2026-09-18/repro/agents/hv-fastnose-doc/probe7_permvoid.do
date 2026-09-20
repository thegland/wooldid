* ===================================================================
* probe7_permvoid.do -- SELF-CONTAINED
* permute() without permtstat/permanalyticse sets `fastnose' internally
* (ado:5448).  The R6 dose-slope void screen (ado:11958) is gated on
* `fastnose'=="".  So does the OBSERVED permutation run also skip the
* screen and report ame rows of exactly 0 WITH permutation p-values?
* Control arm: permanalyticse (fastnose not set) -- screen should run.
* ===================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 24
gen long i = _n
gen int ttre = cond(mod(_n,2)==0, 5, .)
gen double dcoh = exp(0.2*mod(_n,5))        // CONSTANT within cohort
expand 5
bysort i: gen int unit = _n
expand 10
bysort i unit: gen int t = _n
gen double y0 = 0.1*i + 0.05*t + rnormal()
gen byte D = (ttre<. & t>=ttre)
gen double yV = y0 + cond(D==1, 0.5*dcoh, 0)
sort i unit t
compress
save "panel7.dta", replace

dis "=== W1 plain analytic (reference) ==="
use "panel7.dta", clear
ereturn clear
cap noi wooldid yV i t ttre, att robust contreat(dcoh) contreatcontrols(basic) esprelength(2) espostlength(2) esfixedbaseperiod
dis "W1_RC=" _rc
cap confirm matrix e(b_es_ame)
dis "W1_BES_ame_RC=" _rc
cap mat list e(wdidmainresults)

dis "=== W2 permute fast mode ==="
use "panel7.dta", clear
ereturn clear
cap noi wooldid yV i t ttre, att cluster(i) contreat(dcoh) contreatcontrols(basic) esprelength(2) espostlength(2) esfixedbaseperiod permute(shuffle) permreps(30) rseed(11)
dis "W2_RC=" _rc
cap confirm matrix e(b_es_ame)
dis "W2_BES_ame_RC=" _rc
cap mat list e(wdidmainresults)
cap mat list e(perm_p)

dis "=== W3 permute + permanalyticse (fastnose NOT set) ==="
use "panel7.dta", clear
ereturn clear
cap noi wooldid yV i t ttre, att cluster(i) contreat(dcoh) contreatcontrols(basic) esprelength(2) espostlength(2) esfixedbaseperiod permute(shuffle) permreps(30) rseed(11) permanalyticse
dis "W3_RC=" _rc
cap confirm matrix e(b_es_ame)
dis "W3_BES_ame_RC=" _rc
cap mat list e(wdidmainresults)

dis "=== PROBE7 DONE ==="
