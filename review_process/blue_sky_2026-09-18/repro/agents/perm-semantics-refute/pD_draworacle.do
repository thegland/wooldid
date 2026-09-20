*! SKEPTIC probe D -- are the permutation draws what they are documented to be?
*! External oracle: permsave() writes the onset each CLUSTER carried in each replication
*! (the onset_<id> columns).  Take replication r, relabel the onset variable by hand from
*! those columns, refit wooldid with no permute(), and compare the resulting coefficients
*! with the row permsave stored for r (and with e(permdist_b)).  If the stored draws are
*! the relabelled refits, they must agree to machine precision.
*! Also checks the _rep = 0 row against the observed fit (the rep-0 comparator).
*! Self-contained: rebuilds the panel and the permutation run from scratch.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pD_draworacle.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

clear
set obs 10
gen long dist = _n
gen double onset = cond(inlist(dist,4,9), 2017, .)
gen byte leaver = inlist(dist,2,6,10)
gen double a_i = rnormal()*2
expand 10
bysort dist: gen int year = 2010 + _n - 1
drop if leaver & year > 2014
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()
save pD_panel.dta, replace

cap noi wooldid y dist year onset, cluster(dist) att esprelength(2) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(100) ///
    rseed(99) permsave(pD_save) replace
dis "PD_RC=" _rc
dis "PD_reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
matrix _OBSB = e(perm_comparator_b)
matrix list _OBSB, format(%21x)

* ---- pick three SUCCESSFUL replications and one that failed ----
preserve
use pD_save.dta, clear
qui ds _rep onset_*, not
local coefvars `r(varlist)'
qui gen byte allmiss = 1
foreach v of local coefvars {
  qui replace allmiss = 0 if !missing(`v')
}
qui levelsof _rep if _rep >= 1 & allmiss == 0, local(okreps)
local pick1 : word 1 of `okreps'
local pick2 : word 17 of `okreps'
local pick3 : word 55 of `okreps'
dis "PD_PICKS = `pick1' `pick2' `pick3'"
qui levelsof _rep if _rep >= 1 & allmiss == 1, local(failreps)
local pf : word 1 of `failreps'
dis "PD_FAILPICK = `pf'"
global PICKS `pick1' `pick2' `pick3'
global FAILPICK `pf'
restore

foreach R of global PICKS {
  * pull the assignment for replication R
  preserve
  use pD_save.dta, clear
  qui keep if _rep == `R'
  forvalues d = 1/10 {
    local ons`d' = onset_`d'[1]
  }
  local saved_att = main_att[1]
  restore
  dis "PD_REP `R' assignment: `ons1' `ons2' `ons3' `ons4' `ons5' `ons6' `ons7' `ons8' `ons9' `ons10'"

  preserve
  use pD_panel.dta, clear
  gen double onset_perm = .
  forvalues d = 1/10 {
    qui replace onset_perm = `ons`d'' if dist == `d'
  }
  cap noi wooldid y dist year onset_perm, cluster(dist) att esprelength(2) espostlength(2) ///
      jointtests adversarial esfixedbaseperiod
  dis "PD_REFIT`R'_RC=" _rc
  if _rc == 0 {
    local refit_att = e(b)[1, colnumb(e(b), "main_att")]
    dis "PD_CMP rep=`R'  saved_main_att=" %21x `saved_att' "   refit_main_att=" %21x `refit_att'
    dis "PD_DIFF rep=`R'  reldiff=" %12.4e abs(`refit_att' - `saved_att')/max(abs(`saved_att'), 1e-30)
  }
  restore
}

* ---- the failed replication: does its assignment really fail on a hand refit? ----
preserve
use pD_save.dta, clear
qui keep if _rep == $FAILPICK
forvalues d = 1/10 {
  local fons`d' = onset_`d'[1]
}
restore
dis "PD_FAILREP $FAILPICK assignment: `fons1' `fons2' `fons3' `fons4' `fons5' `fons6' `fons7' `fons8' `fons9' `fons10'"
preserve
use pD_panel.dta, clear
gen double onset_perm = .
forvalues d = 1/10 {
  qui replace onset_perm = `fons`d'' if dist == `d'
}
cap noi wooldid y dist year onset_perm, cluster(dist) att esprelength(2) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod
dis "PD_FAILREFIT_RC=" _rc
restore

* ---- the _rep = 0 row against the observed fit ----
preserve
use pD_save.dta, clear
qui keep if _rep == 0
dis "PD_REP0_main_att = " %21x main_att[1]
restore
dis "PD_OBS_comparator_main_att = " %21x _OBSB[1, colnumb(_OBSB, "main_att")]

dis "PD_DONE"
