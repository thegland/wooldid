* probe05 -- FAILING REPLICATIONS.  (a) whole-replication failures under a
* permute(probabilities) design that can draw zero treated clusters; (b) element-level
* missingness from horizons a reassigned onset makes inestimable.  In both cases:
* is the denominator the help's "number of COMPLETED replications"?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"

* ============ (a) probabilities design: P(no treated cluster) = .85^12 = .142 ==========
local dN 12
local dT 8
local dons 2004
local dnever 9
local dtau 1.5
local dsd 0.5
local dseed 777
include "dgp.do"
tab ttre, m
tempfile basea
save `basea'
local PP (ttre 2004 .)
forvalues u = 1/12 {
  local PP `PP' (`u' 0.15 0.85)
}
use `basea', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) permute(probabilities) permprobs(`PP') permreps(200) rseed(31) permsave(f_prob) replace
dis "RC_PROB=" _rc
dis "reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed) " npossible=" e(perm_npossible)
local htag FPROB
local habsp 0
local hexact 0
include "harvest.do"
matrix list e(perm_n)
mata:
  D = st_matrix("_PD")
  allmiss = 0
  anymiss  = 0
  for (r=1;r<=rows(D);r++) {
    if (sum(D[r,.] :>= .) == cols(D)) allmiss = allmiss + 1
    else if (sum(D[r,.] :>= .) > 0) anymiss = anymiss + 1
  }
  printf("FPROB dump rows=%g  all-missing rows=%g  partially-missing rows=%g\n", rows(D), allmiss, anymiss)
  printf("FPROB e(perm_completed)=%g  e(perm_failed)=%g\n", st_numscalar("e(perm_completed)"), st_numscalar("e(perm_failed)"))
end

* ============ (b) element-level: a long pre-window some draws cannot fill ==========
local dN 20
local dT 9
local dons 2003 2007
local dnever 8
local dtau 1.0
local dsd 0.6
local dseed 555
include "dgp.do"
tab ttre, m
tempfile baseb
save `baseb'
use `baseb', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(6) espostlength(5) permute(shuffle) permreps(199) rseed(41) permsave(f_elem) replace
dis "RC_ELEM=" _rc
dis "reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
local htag FELEM
local habsp 0
local hexact 0
include "harvest.do"
matrix list e(perm_n)
matrix list e(perm_p), format(%12.8f)
dis "PROBE05_DONE"
