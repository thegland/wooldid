*! r04 — REFUTE FCW-04: is the reps() p-value-resolution caution really unreachable
*! without cluster()?  Six runs on one self-generated panel, reps(10) throughout.
*! Also tests bootdisplay(wald) and boottstat legs, which FCW-04 did not.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 5551
set sortseed 5551
set obs 12
gen long cty = _n
gen onset = .
replace onset = 4 if cty <= 4
replace onset = 6 if cty > 4 & cty <= 8
expand 4
sort cty
gen long unit = _n
gen double uf = rnormal(0,1)
expand 7
bysort unit (cty): gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double y = 2 + uf + .1*t + 1*D + rnormal(0,1)
keep y unit cty t onset

dis "===== RUN1: cell-weighted estimands, bootstrap(standard) reps(10), NO cluster()"
cap noi wooldid y unit t onset, att att_i att_it bootstrap(standard) reps(10) rseed(9) dots(0)
dis "RUN1_RC=" _rc

dis "===== RUN2: identical plus cluster(cty)"
cap noi wooldid y unit t onset, att att_i att_it cluster(cty) bootstrap(standard) reps(10) rseed(9) dots(0)
dis "RUN2_RC=" _rc

dis "===== RUN3: att only, NO cluster()"
cap noi wooldid y unit t onset, att bootstrap(standard) reps(10) rseed(9) dots(0)
dis "RUN3_RC=" _rc

dis "===== RUN4: att only, NO cluster(), bootdisplay(wald)"
cap noi wooldid y unit t onset, att bootstrap(standard) reps(10) rseed(9) bootdisplay(wald) dots(0)
dis "RUN4_RC=" _rc

dis "===== RUN5: att only, NO cluster(), boottstat"
cap noi wooldid y unit t onset, att bootstrap(standard) reps(10) rseed(9) boottstat dots(0)
dis "RUN5_RC=" _rc

dis "===== RUN6: bootstrap(bayes) reps(10) NO cluster()  (bayes needs cluster? record rc)"
cap noi wooldid y unit t onset, att bootstrap(bayes) reps(10) rseed(9) dots(0)
dis "RUN6_RC=" _rc
dis "===== DONE"
