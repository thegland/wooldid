* probe07_cautions.do -- FP2's asymmetry claim, tested directly and NON-quietly.
* Does permute()+permonsetcohorts print anything when a covariate names i?
* Does the bootstrap's A10 caution actually print on the comparable specification?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
which wooldid

dis "@@@@ A: permute + permonsetcohorts + fe(coh)  (fe names i)"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts fe(coh)
dis "A_RC=" _rc

dis "@@@@ B: permute + permonsetcohorts + controls(cohbase) (constant within i, not named i)"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts controls(cohbase)
dis "B_RC=" _rc

dis "@@@@ C: bootstrap(standard) cluster(coh)==i, att_i, fe(coh)  -- the A10 nesting leg"
use panA.dta, clear
cap noi wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(20) fe(coh)
dis "C_RC=" _rc

dis "@@@@ D: bootstrap(standard) cluster(coh)==i, att_i, controls(cohbase) -- constant within i, NOT named i"
use panA.dta, clear
cap noi wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(20) controls(cohbase)
dis "D_RC=" _rc

dis "@@@@ E: bootstrap(standard) cluster(state) FINER than i, att_i, fe(coh)"
use panA.dta, clear
cap noi wooldid y coh t onset, att_i cluster(state) bootstrap(standard) reps(20) fe(coh)
dis "E_RC=" _rc
dis "PROBE07 DONE"
