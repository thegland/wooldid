*=============================================================================
* repro_FP2_cautions.do   SELF-CONTAINED.  Tests FP-2's asymmetry claim directly.
*   A  permute + permonsetcohorts + fe(coh)                 -> caution?
*   B  permute + permonsetcohorts + controls(cohbase)        -> caution?
*   C  bootstrap(standard), cluster(coh)==i, att_i, fe(coh)  -> the A10 caution?
*   D  bootstrap(standard), cluster(coh)==i, att_i, controls(cohbase)  -> caution?
*   E  bootstrap(standard), cluster(state) finer than i, att_i, fe(coh) -> caution?
* Grep the log for: "names the i variable"
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
clear
set obs 7
gen int state = _n
gen byte coh = 1
replace coh = 2 if inlist(state,4,5)
replace coh = 3 if inlist(state,6,7)
gen int onset = .
replace onset = 2005 if coh==1
replace onset = 2009 if coh==2
gen int t0 = 1998
replace t0 = 2000 if inlist(state,2,6)
replace t0 = 1999 if inlist(state,3,7)
replace t0 = 2001 if state==4
gen int t1 = 2015
replace t1 = 2014 if state==2
replace t1 = 2013 if state==4
replace t1 = 2012 if state==6
gen double a_s = .
replace a_s = 1.20 in 1
replace a_s = 0.40 in 2
replace a_s = 2.10 in 3
replace a_s = -0.60 in 4
replace a_s = 1.70 in 5
replace a_s = 0.10 in 6
replace a_s = -1.30 in 7
expand t1-t0+1
bys state: gen int t = t0 + _n - 1
drop t0 t1
gen double cohbase = cond(coh==1,0.80,cond(coh==2,1.45,2.10))
gen double post = (onset<. & t>=onset)
set seed 771131
gen double y = a_s + 0.35*cohbase + 0.45*rnormal() + 0.22*post
compress
tempfile PAN
save `PAN'
dis _n "@@@@ A"
use `PAN', clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts fe(coh)
dis "A_RC=" _rc
dis _n "@@@@ B"
use `PAN', clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts controls(cohbase)
dis "B_RC=" _rc
dis _n "@@@@ C"
use `PAN', clear
cap noi wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(20) fe(coh)
dis "C_RC=" _rc
dis _n "@@@@ D"
use `PAN', clear
cap noi wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(20) controls(cohbase)
dis "D_RC=" _rc
dis _n "@@@@ E"
use `PAN', clear
cap noi wooldid y coh t onset, att_i cluster(state) bootstrap(standard) reps(20) fe(coh)
dis "E_RC=" _rc
dis _n "REPRO DONE"
