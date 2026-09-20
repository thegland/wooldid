* probe14 -- ROW ALIGNMENT of the permutation draw matrices when replications FAIL.
* e(permdist_b) is compacted (failed rows dropped) but e(perm_jointW)/e(perm_advW) are
* written at full accumulator length.  If so, the two cannot be paired row by row, and
* neither indexes permsave's _rep.  The bootstrap side documents its compaction
* ("_rep indexes the draws exactly as e(bootdist) does after its missing rows are
* dropped"); the permutation side has no such sentence.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
local dN 12
local dT 8
local dons 2004
local dnever 9
local dtau 1.5
local dsd 0.5
local dseed 777
include "dgp.do"
local PP (ttre 2004 .)
forvalues u = 1/12 {
  local PP `PP' (`u' 0.15 0.85)
}
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) jointtests adversarial esfixedbaseperiod permute(probabilities) permprobs(`PP') permreps(120) rseed(31) permsave(a_prob) replace
dis "RC_A=" _rc
dis "reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
foreach m in permdist_b permdist perm_jointW perm_advW perm_sgW {
  cap confirm matrix e(`m')
  if _rc == 0 {
    matrix TMP = e(`m')
    dis "`m' : " rowsof(TMP) " x " colsof(TMP)
  }
  else dis "`m' : absent"
}
matrix JJ = e(perm_jointW)
matrix DD = e(permdist_b)
mata:
  J = st_matrix("JJ")
  D = st_matrix("DD")
  printf("jointW rows=%g non-missing=%g ; permdist_b rows=%g\n", rows(J), sum(J[.,1] :< .), rows(D))
  printf("ALIGNED row-for-row? %g\n", rows(J) == rows(D))
end
preserve
use a_prob.dta, clear
dis "permsave rows = " _N
restore
dis "PROBE14_DONE"
