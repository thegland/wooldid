* =================================================================
* probe_permuniv.do -- a probe NOBODY suggested: does permute() share
* the bootstrap's sample-universe defect?
* wooldid.ado:665-669 asserts that "the estimation sample, the permutation
* design and the bootstrap frame agree" because all copies of the user `if`
* share the exact-one rule. But ado:682 additionally zeroes `ifvar` for
* missing y/i/t/w, and permctl's `insmp` (ado:5957) / bootrun's `smp`
* (ado:8917) are built from the RAW if only. Test: units with an entirely
* missing outcome should then appear in the permutation ASSIGNMENT universe.
* permsave() writes one onset_<clustervalue> column per cluster, so the
* column count is a direct read-out of that universe.
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_permuniv.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 767676
set sortseed 767676

* 12 REAL units (4 treated 2004, 4 treated 2006, 4 never) + 6 PHANTOM units
* (3 treated 2004, 3 never) whose outcome is missing on every row.
set obs 18
gen long unit = _n
gen byte phantom = unit > 12
gen int ttre = .
replace ttre = 2004 if unit <= 4
replace ttre = 2006 if unit > 4 & unit <= 8
replace ttre = 2004 if unit > 12 & unit <= 15
gen double ufe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.6, 0))
gen double y  = 1 + ufe + 0.1*(year-1999) + te + rnormal()*0.6
replace y = . if phantom
save "sk_pu_A.dta", replace
drop if phantom
save "sk_pu_B.dta", replace

foreach arm in A B {
  use "sk_pu_`arm'.dta", clear
  cap erase "skpu`arm'.dta"
  dis "########## PERM ARM `arm' ##########"
  cap noi wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(40) rseed(8642) permsave(skpu`arm') dots(0)
  dis "PU_`arm'_RC=" _rc
  dis "PU_`arm'_N=" e(N) "  NCLUST=" e(N_clust)
  dis "PU_`arm'_ATT=" %21x el(e(b),1,1)
  cap confirm matrix e(permdist)
  preserve
  use "skpu`arm'.dta", clear
  qui ds onset_*
  local oc `r(varlist)'
  local no : word count `oc'
  dis "PU_`arm'_NONSETCOLS=" `no'
  cap confirm variable onset_14
  dis "PU_`arm'_HAS_ONSET_14_RC=" _rc
  cap qui sum onset_14 if _n > 1
  dis "PU_`arm'_ONSET14_NONMISS=" r(N)
  restore
}

dis "PROBE_PERMUNIV_DONE"
