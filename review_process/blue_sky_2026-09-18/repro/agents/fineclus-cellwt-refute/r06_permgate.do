*! r06 — REFUTE FCW-06: does T1's permute() relief reach "coarse i, finer cluster"?
*! FCW-06 tested ONE shape (i = state, several states per onset) and concluded the
*! relief never reaches the motivating case. Here BOTH shapes are tested on the same
*! panel: (X) i finer than the onset grouping, and (Y) i EXACTLY the onset grouping
*! with a strictly finer cluster() -- which is what the help tells the user to build.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 6061
set sortseed 6061

* 12 states x 4 counties, T=8. states 1-4 onset 4, 5-8 onset 6, 9-12 never.
set obs 12
gen long state = _n
gen onset = .
replace onset = 4 if state <= 4
replace onset = 6 if state > 4 & state <= 8
expand 4
sort state
gen long county = _n
gen double uf = rnormal(0,1)
expand 8
bysort county (state): gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double y = 2 + uf + .1*t + 1*D + rnormal(0,1)
* the onset grouping itself: never-treated is ONE level
gen long cohort = cond(missing(onset), 1, cond(onset==4, 2, 3))
keep y state county cohort t onset

dis "===== X1: i = state (finer than the onset grouping), cluster(county), permute(shuffle)"
cap noi wooldid y state t onset, att att_i cluster(county) permute(shuffle) permreps(30) rseed(1) dots(0)
dis "X1_RC=" _rc

dis "===== X2: same plus permonsetcohorts"
cap noi wooldid y state t onset, att att_i cluster(county) permute(shuffle) permreps(30) permonsetcohorts rseed(1) dots(0)
dis "X2_RC=" _rc

dis "===== Y1: i = cohort (IS the onset grouping), cluster(county) strictly finer, NO option"
cap noi wooldid y cohort t onset, att att_i cluster(county) permute(shuffle) permreps(30) rseed(1) dots(0)
dis "Y1_RC=" _rc

dis "===== Y2: i = cohort, cluster(county), permonsetcohorts  <- the T1 remedy"
cap noi wooldid y cohort t onset, att att_i cluster(county) permute(shuffle) permreps(30) permonsetcohorts rseed(1) dots(0)
dis "Y2_RC=" _rc
dis "Y2 perm_irebuilt = " e(perm_irebuilt)

dis "===== Z: i = state, cluster(county), bootstrap(standard) (the control FCW-06 ran)"
cap noi wooldid y state t onset, att att_i cluster(county) bootstrap(standard) reps(30) rseed(1) dots(0)
dis "Z_RC=" _rc
dis "===== DONE"
