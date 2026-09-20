* probe2a_caution.do -- two things the T2 relaxation newly exposes.
*  Q1 bootstrap(standard) with a cell-weighted estimand and NO cluster() used to be
*     refused at ado:6353 with error 198.  Does it now run?
*  Q2 does the reps() p-value-resolution caution still print on that run?  RF-2
*     (commit 813c414, 2026-09-13) moved that caution from an unconditional site into
*     the body of `if "`clus'" != "" {' (pinned build :6459-:6501, caution :6473-6478),
*     so an unclustered bootstrap should now print nothing at all.  The clustered run
*     below is the paired control: same reps(), same data, same estimands.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

clear
qui set obs 40
gen long unit = _n
gen int ttre = .
qui replace ttre = 3 if mod(unit,4)==1
qui replace ttre = 5 if mod(unit,4)==2
gen double tau = 0.3 + 0.05*unit if !missing(ttre)
qui replace tau = 0 if missing(tau)
qui expand 6
bysort unit: gen int t = _n
set seed 901
qui bysort unit (t): gen double mu = rnormal(0,0.5) if _n==1
qui by unit: replace mu = mu[1]
gen double eps = rnormal(0,0.8)
sort unit t
gen double y = tau*(t>=ttre & ttre<.) + 0.03*t + mu + eps

dis "@@@ RUN1: NO cluster(), cell-weighted estimands, reps(10)"
cap noi wooldid y unit t ttre, att att_i att_it bootstrap(standard) reps(10) rseed(11) dots(0)
dis "@@@ RUN1_RC=" _rc
dis "@@@ RUN2: SAME but cluster(unit) added"
cap noi wooldid y unit t ttre, att att_i att_it cluster(unit) bootstrap(standard) reps(10) rseed(11) dots(0)
dis "@@@ RUN2_RC=" _rc
dis "@@@ RUN3: NO cluster(), att only (no cell-weighted estimand), reps(10)"
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(11) dots(0)
dis "@@@ RUN3_RC=" _rc
