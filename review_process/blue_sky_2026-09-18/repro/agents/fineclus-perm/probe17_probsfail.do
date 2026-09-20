* probe17: permute(probabilities) failure censoring -- is it rebuild-specific?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* (a) WITH the rebuild: cluster(state) finer than i
use "dataB.dta", clear
local pp (onset 4 6 .)
forval s = 1/6 {
  local pp `pp' (`s' 0.34 0.33 0.33)
}
cap noi wooldid y coh t onset, att cluster(state) permute(probabilities) permprobs(`pp') permreps(150) rseed(99) permonsetcohorts permsave(PF_reb) replace
dis "PF_REB rc=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att) " incomplete=" e(perm_incomplete)

* (b) WITHOUT the rebuild: cluster(coh) == i
use "dataB.dta", clear
local pp2 (onset 4 6 .)
foreach c in 1 2 3 {
  local pp2 `pp2' (`c' 0.34 0.33 0.33)
}
cap noi wooldid y coh t onset, att cluster(coh) permute(probabilities) permprobs(`pp2') permreps(150) rseed(99) permsave(PF_nor) replace
dis "PF_NOR rc=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att) " incomplete=" e(perm_incomplete)

* (c) what do the failing draws look like?
use "PF_reb.dta", clear
gen byte has4 = 0
gen byte has6 = 0
gen byte hasNT = 0
forval s = 1/6 {
  qui replace has4 = 1 if onset_`s' == 4
  qui replace has6 = 1 if onset_`s' == 6
  qui replace hasNT = 1 if missing(onset_`s')
}
gen byte failed = missing(main_att)
qui count if _rep > 0
dis "PF_TOTREP=" r(N)
qui count if _rep > 0 & failed
dis "PF_FAILED=" r(N)
qui count if _rep > 0 & hasNT == 0
dis "PF_NO_NEVERTREATED_DRAWN=" r(N)
qui count if _rep > 0 & has4 == 0 & has6 == 0
dis "PF_NO_TREATED_DRAWN=" r(N)
qui count if _rep > 0 & (has4 == 0 | has6 == 0)
dis "PF_LEVELSET_LOST=" r(N)
tab hasNT failed if _rep > 0
tab has4 failed if _rep > 0
