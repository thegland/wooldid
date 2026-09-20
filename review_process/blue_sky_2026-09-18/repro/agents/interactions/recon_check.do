clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* tiny staggered panel: 40 units, 5 clusters of 8 units, 10 periods
clear
set obs 40
gen unit = _n
gen clus = ceil(_n/8)
gen cohort = 1 + mod(_n,4)          // 4 cohorts, cross-cuts clus
gen onset  = cond(cohort==1, ., 3 + cohort)
expand 10
bysort unit: gen t = _n
gen y = 0.3*unit + 0.2*t + rnormal() + 1.0*(onset<. & t>=onset)
gen double w = 1 + mod(unit,3)

dis "SENTINEL A1: clustered bootstrap reps(10) -- expect resolution caution"
cap noi wooldid y unit t onset, att cluster(clus) bootstrap(standard) reps(10) rseed(7)
dis "A1 rc=" _rc

dis "SENTINEL A2: UNCLUSTERED bootstrap reps(10) -- resolution caution?"
cap noi wooldid y unit t onset, att bootstrap(standard) reps(10) rseed(7)
dis "A2 rc=" _rc

dis "SENTINEL B: att_it, bootstrap(standard), NO cluster() -- T2 lifted refusal"
cap noi wooldid y unit t onset, att_it bootstrap(standard) reps(20) rseed(7)
dis "B rc=" _rc
dis "B boot_failed=" e(boot_failed) " boot_completed=" e(boot_completed) " nclus=" e(boot_nclus)

dis "SENTINEL C: att_it, bootstrap(bayes), NO cluster() -- companion refusal still?"
cap noi wooldid y unit t onset, att_it bootstrap(bayes) reps(20) rseed(7)
dis "C rc=" _rc

dis "SENTINEL D: regtol under unconditionalse (silent no-op?)"
cap noi wooldid y unit t onset, att unconditionalse robust regtol(3)
dis "D rc=" _rc

dis "SENTINEL E: contreatpoly(3) with NO contreat() -- silent no-op?"
cap noi wooldid y unit t onset, att contreatpoly(3)
dis "E rc=" _rc

dis "SENTINEL F: robust alongside cluster() -- silently dropped?"
cap noi wooldid y unit t onset, att robust cluster(clus)
dis "F rc=" _rc
dis "F vcetype=" e(vcetype) " clustvar=" e(clustvar)

dis "SENTINEL ZEND"
