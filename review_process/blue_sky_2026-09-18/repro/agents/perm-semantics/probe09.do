* probe09 -- permonsetcohorts (shipped 2026-09-16, assessment deferred).
* i IS the onset grouping and cluster() is finer than i, so each replication rebuilds i
* from the onsets it drew.  Two questions:
*  (1) do the reported p-values conform to the documented formulas on those draws?
*  (2) IS each draw the estimate wooldid would compute after rebuilding i by hand?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
mata:
void cmprep2(real scalar rr, string scalar rbname)
{
  real rowvector rb, dr
  real scalar j, d, worst
  string matrix nm, nd
  rb = st_matrix(rbname)
  dr = st_matrix("DUMPB")[rr+1, .]
  nm = st_matrixcolstripe(rbname)
  nd = st_matrixcolstripe("DUMPB")
  worst = 0
  for (j=1;j<=cols(rb);j++) {
    d = reldif(rb[1,j], dr[1,j])
    if (rb[1,j]>=. & dr[1,j]>=.) d = 0
    if (d > worst) worst = d
    printf("  rep %g col %g refit=%21x dump=%21x reldif=%g\n", rr, j, rb[1,j], dr[1,j], d)
  }
  printf("ONSETCOH REFIT rep=%g WORST reldif = %g\n", rr, worst)
}
end
* ---- data: i is EXACTLY the onset grouping; cluster is the unit, finer than i ----
local dN 24
local dT 8
local dons 2004 2006
local dnever 8
local dtau 1.2
local dsd 0.7
local dseed 313
include "dgp.do"
gen long icoh = .
replace icoh = 1 if ttre == 2004
replace icoh = 2 if ttre == 2006
replace icoh = 3 if missing(ttre)
tab icoh ttre, m
tempfile base
save `base'
use `base', clear
cap noi wooldid y icoh t ttre, cluster(unit) att esprelength(2) espostlength(2) permute(shuffle) permonsetcohorts permreps(99) rseed(17) permsave(oc_def) replace
dis "RC_OC=" _rc
dis "perm_irebuilt = " e(perm_irebuilt)
local htag ONSETCOH
local habsp 0
local hexact 0
include "harvest.do"
matrix list e(perm_n)
matrix list e(perm_p), format(%12.8f)

* ---- the refit check ----
use oc_def.dta, clear
qui ds onset_*
local ons `r(varlist)'
qui ds onset_* _rep, not
local bc `r(varlist)'
dis "dump rows=" _N "  onsetcols=`: word count `ons''  coefcols=`: word count `bc''"
mata: st_matrix("DUMPB", st_data(., "`bc'"))
mata: st_matrix("DUMPO", st_data(., "`ons'"))
local clusnames
foreach v of local ons {
  local clusnames `clusnames' `=subinstr("`v'", "onset_", "", 1)'
}
foreach rr in 0 1 3 9 42 {
  use `base', clear
  local k 0
  foreach cv of local clusnames {
    local ++k
    local newons = DUMPO[`=`rr'+1', `k']
    qui replace ttre = `newons' if unit == `cv'
  }
  * rebuild i from the drawn onsets exactly as the option asserts: one level of i per
  * distinct onset, the never-treated a single level of their own
  cap drop icoh2
  qui egen long icoh2 = group(ttre), missing
  cap noi wooldid y icoh2 t ttre, cluster(unit) att esprelength(2) espostlength(2) fastnose
  dis "ONSETCOH REFIT rep=`rr' rc=" _rc
  if _rc == 0 {
    matrix RB = e(b)
    mata: cmprep2(`rr', "RB")
  }
}
dis "PROBE09_DONE"
