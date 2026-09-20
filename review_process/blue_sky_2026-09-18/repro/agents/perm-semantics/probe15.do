* probe15 -- permexact with an element that SOME enumerated assignments cannot estimate.
* 8 clusters, 2 treated -> M = 28.  A subgroup level carried by only two clusters is
* absent from the 15 assignments that treat neither of them, so its coefficients have
* n = 13 < M.  The help: "an element that some enumerated assignments could not estimate
* keeps the exact share over the assignments that completed - its floor is 2/n on its own
* count n, and a caution names it.  Only an element whose reference set does not contain
* the observed statistic at all ... falls back to the add-one form, under its own caution."
* So: plain shares over the reduced n, e(perm_p_se) still missing, e(perm_incomplete)
* naming the element, e(perm_addone) empty.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
clear
set obs 8
gen long unit = _n
gen double ttre = cond(inlist(unit,1,2), 2004, .)
gen byte sg = cond(inlist(unit,1,5), 1, 0)
gen double ai = rnormal()
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ai + 0.25*(t-2000) + 2.0*post + 0.5*rnormal()
tab unit sg
tempfile base
save `base'
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(1) espostlength(1) subgroups(sg) permute(shuffle) permexact permsave(i_exact) replace
dis "RC_IEX=" _rc
dis "npossible=" e(perm_npossible) " completed=" e(perm_completed) " failed=" e(perm_failed)
local inc `e(perm_incomplete)'
local add `e(perm_addone)'
dis "perm_incomplete = [`inc']"
dis "perm_addone     = [`add']"
local htag IEXACT
local habsp 0
local hexact 1
include "harvest.do"
matrix list e(perm_n)
matrix list e(perm_p), format(%12.8f)
matrix list e(perm_p_se)
matrix list e(perm_p_ci)
mata:
  D = st_matrix("_PD")
  nm = st_matrixcolstripe("_PD")
  for (j=1;j<=cols(D);j++) printf("col %g %-26s non-missing draws=%g\n", j, nm[j,2], sum(D[.,j] :< .))
end
dis "PROBE15_DONE"
