* probe10: contreat() under the rebuild; cross-cutting cluster; permprobs with rebuild
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* ---------- A. cross-cutting cluster under permonsetcohorts -> must be refused ----------
use "dataB.dta", clear
gen byte xclus = cond(inlist(state,1,3),1,cond(inlist(state,2,5),2,3))
cap noi wooldid y coh t onset, att cluster(xclus) permute(shuffle) permexact permonsetcohorts
dis "XC_RC=" _rc

* ---------- B. permprobs with the rebuild ----------
use "dataB.dta", clear
local pp (onset 4 6 .)
forval s = 1/6 {
  local pp `pp' (`s' 0.34 0.33 0.33)
}
cap noi wooldid y coh t onset, att cluster(state) permute(probabilities) permprobs(`pp') permreps(150) rseed(99) permonsetcohorts permsave(PR_B) replace
dis "PR_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att)
cap use "PR_B.dta", clear
cap count if missing(main_att)
dis "PR_NMISSING=" r(N)
* how many drawn assignments lose an onset level entirely?
cap {
  gen byte has4 = 0
  gen byte has6 = 0
  gen byte hasNT = 0
  forval s = 1/6 {
    qui replace has4 = 1 if onset_`s' == 4
    qui replace has6 = 1 if onset_`s' == 6
    qui replace hasNT = 1 if missing(onset_`s')
  }
  count if _rep > 0 & (has4 == 0 | has6 == 0)
  dis "PR_LEVELSETLOST=" r(N)
  count if _rep > 0 & hasNT == 0
  dis "PR_NONTdrawn=" r(N)
  count if _rep > 0 & missing(main_att)
  dis "PR_MISSINGB=" r(N)
  tab has4 has6 if _rep > 0
}

* ---------- C. contreat() under the rebuild ----------
use "dataB.dta", clear
set seed 4242
gen double dose = runiform()*2 if !missing(onset)
replace dose = 0 if missing(onset)
cap noi wooldid y coh t onset, att cluster(state) contreat(dose) permute(shuffle) permexact permonsetcohorts permsave(CT_A) replace
dis "CT_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att)
