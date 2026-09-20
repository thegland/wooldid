* probe8_bootvoid.do -- SELF-CONTAINED
* bootstrap() runs its replications (and its rep-0 comparator) under fastnose
* while the OBSERVED run is not a fastnose run.  On the R6-void geometry the
* observed run scrubs ame to missing but every replication returns ame = 0.
* What does the rewired table report?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
clear
set obs 30
gen long i = _n
gen int ttre = cond(mod(_n,2)==0, 6, .)
gen double dcoh = exp(0.2*mod(_n,5))
expand 6
bysort i: gen int unit = _n
expand 12
bysort i unit: gen int t = _n
gen double y0 = 0.1*i + 0.05*t + rnormal()
gen byte D = (ttre<. & t>=ttre)
gen double yV = y0 + cond(D==1, 0.5*dcoh, 0)
sort i unit t
compress
save "panel8.dta", replace
dis "TRUE_AME_RANGE:"
qui sum dcoh if D==1
dis "dcoh_min=" r(min) " dcoh_max=" r(max) " => true dy/ddose = 0.5 on treated cells"

dis "=== X1 bootstrap(standard) on the VOID geometry ==="
use "panel8.dta", clear
ereturn clear
cap noi wooldid yV i t ttre, att robust contreat(dcoh) contreatcontrols(basic) esprelength(2) espostlength(2) esfixedbaseperiod bootstrap(standard) reps(20) rseed(3)
dis "X1_RC=" _rc
cap confirm matrix e(b_es_ame)
dis "X1_BES_ame_RC=" _rc
cap mat list e(wdidmainresults)
dis "=== PROBE8 DONE ==="
