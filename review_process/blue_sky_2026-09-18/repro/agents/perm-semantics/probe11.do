* probe11 -- ELEMENT-LEVEL missingness.  The help: "a subgroup level absent from some
* replications contributes fewer draws (e(perm_n))".  A subgroup level carried by only
* three of twenty-four clusters is absent from ~28% of shuffles, so its coefficients are
* missing in those draws.  Does every reported quantity use THAT element's own n?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
local dN 24
local dT 8
local dons 2004 2006
local dnever 8
local dtau 1.2
local dsd 0.7
local dseed 8181
include "dgp.do"
* a rare subgroup level: only three units carry sg==1
replace sg = 0
replace sg = 1 if inlist(unit, 3, 11, 19)
tab sg
tempfile base
save `base'
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(199) rseed(29) permstepdown permsave(e_sg) replace
dis "RC_ESG=" _rc
dis "reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
local htag ELEMSG
local habsp 0
local hexact 0
include "harvest.do"
matrix list e(perm_n)
matrix list e(perm_p), format(%12.8f)
matrix list e(perm_p_se), format(%12.8f)
mata:
  D = st_matrix("_PD")
  nm = st_matrixcolstripe("_PD")
  for (j=1;j<=cols(D);j++) printf("col %g %-26s non-missing draws=%g\n", j, nm[j,2], sum(D[.,j] :< .))
end
* the per-subgroup step-down matrices, by name
foreach nmx in perm_stepdown perm_stepdown_pre perm_stepdown_s0 perm_stepdown_s1 perm_stepdown_pre_s0 perm_stepdown_pre_s1 {
  cap confirm matrix e(`nmx')
  dis "`nmx' posted rc=" _rc
  if _rc == 0 {
    matrix TMP = e(`nmx')
    matrix list TMP, format(%12.8f)
  }
}
dis "PROBE11_DONE"
