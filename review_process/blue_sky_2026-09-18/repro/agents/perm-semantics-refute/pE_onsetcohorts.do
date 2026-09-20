*! SKEPTIC probe E -- NEW ground: permonsetcohorts (T1, shipped 2026-09-16, whose own
*! assessment suite was deferred and has not run).  Under that option each replication
*! REBUILDS i from the onsets it drew.  Question: are the stored draws the refits of the
*! rebuilt design, or of something else?
*! Oracle: permsave() writes the onset each CLUSTER drew.  For a chosen replication,
*! rebuild the onset variable AND rebuild i by the documented map (drawn onset -> the
*! OBSERVED i code of the cohort whose observed onset equals it; never-treated -> the
*! never-treated cohort's own code, ado:7578-7590), then refit wooldid by hand and compare.
*! Also a negative control: with cluster() == i the rebuild must be a declared no-op
*! (e(perm_irebuilt) == 0).
*! Design: i IS the onset grouping (4 cohorts), cluster() is the finer unit id.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pE_onsetcohorts.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 24680
set sortseed 24680
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

* --------- 20 plants in 4 onset cohorts of 5, 2010-2019 ----------
clear
set obs 20
gen long plant = _n
gen byte coh = ceil(_n/5)          // 1..4 : i IS the onset grouping
gen double onset = .
replace onset = 2014 if coh == 1
replace onset = 2016 if coh == 2
replace onset = 2018 if coh == 3
* coh 4 never treated
gen double a_i = rnormal()*2
expand 10
bysort plant: gen int year = 2010 + _n - 1
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()
save pE_panel.dta, replace

* ---------------- the run under test: cluster finer than i -----------------
cap noi wooldid y coh year onset, cluster(plant) att esprelength(1) espostlength(1) ///
    permonsetcohorts permute(shuffle) permreps(60) rseed(808) permsave(pE_save) replace
dis "PE_RC=" _rc
dis "PE_irebuilt=" e(perm_irebuilt)
dis "PE_reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
matrix _OBSB = e(perm_comparator_b)
local obnames : colnames _OBSB
dis "PE_OBNAMES = `obnames'"

* ---------------- negative control: cluster() == i -----------------
preserve
use pE_panel.dta, clear
cap noi wooldid y coh year onset, cluster(coh) att esprelength(1) espostlength(1) ///
    permonsetcohorts permute(shuffle) permreps(20) rseed(808)
dis "PE_CTRL_RC=" _rc
dis "PE_CTRL_irebuilt=" e(perm_irebuilt)
restore

* ---------------- the oracle -----------------
preserve
use pE_save.dta, clear
qui ds _rep onset_*, not
local coefvars `r(varlist)'
dis "PE_COEFCOLS = `coefvars'"
qui gen byte allmiss = 1
foreach v of local coefvars {
  qui replace allmiss = 0 if !missing(`v')
}
qui levelsof _rep if _rep >= 1 & allmiss == 0, local(okreps)
local r1 : word 1 of `okreps'
local r2 : word 9 of `okreps'
local r3 : word 25 of `okreps'
global RPICKS `r1' `r2' `r3'
dis "PE_PICKS = `r1' `r2' `r3'"
restore

foreach R of global RPICKS {
  preserve
  use pE_save.dta, clear
  qui keep if _rep == `R'
  forvalues p = 1/20 {
    local o`p' = onset_`p'[1]
  }
  local sv_att  = main_att[1]
  restore
  dis "PE_REP `R' onsets: `o1' `o2' `o3' `o4' `o5' | `o6' `o7' `o8' `o9' `o10' | `o11' `o12' `o13' `o14' `o15' | `o16' `o17' `o18' `o19' `o20'"

  preserve
  use pE_panel.dta, clear
  gen double repons = .
  forvalues p = 1/20 {
    qui replace repons = `o`p'' if plant == `p'
  }
  * the documented rebuild: drawn onset -> observed i code of that onset's cohort
  gen double repi = .
  qui replace repi = 1 if repons == 2014
  qui replace repi = 2 if repons == 2016
  qui replace repi = 3 if repons == 2018
  qui replace repi = 4 if missing(repons)
  qui assert !missing(repi)
  cap noi wooldid y repi year repons, cluster(plant) att esprelength(1) espostlength(1)
  dis "PE_REFIT`R'_RC=" _rc
  if _rc == 0 {
    local rf_att = e(b)[1, colnumb(e(b), "main_att")]
    dis "PE_CMP rep=`R' saved=" %21x `sv_att' " refit=" %21x `rf_att'
    dis "PE_DIFF rep=`R' reldiff=" %12.4e abs(`rf_att' - `sv_att')/max(abs(`sv_att'),1e-30)
  }
  * control: what would a refit that did NOT rebuild i give?
  cap noi wooldid y coh year repons, cluster(plant) att esprelength(1) espostlength(1)
  dis "PE_NOREBUILD`R'_RC=" _rc
  if _rc == 0 {
    local nb_att = e(b)[1, colnumb(e(b), "main_att")]
    dis "PE_NOREBUILD rep=`R' saved=" %21x `sv_att' " norebuild=" %21x `nb_att'
    dis "PE_NBDIFF rep=`R' reldiff=" %12.4e abs(`nb_att' - `sv_att')/max(abs(`sv_att'),1e-30)
  }
  restore
}

dis "PE_DONE"
