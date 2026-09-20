* probe04 -- permstepdown: recompute the max-|studentized| step-down from the dump
* and check every invariant the help claims (monotone, FWER family, p_raw metric,
* separate pre/post families, listwise-complete reference set).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
mata:
// independent re-implementation of the DOCUMENTED step-down
void sdcheck(string scalar Dname, string scalar CMname, string scalar cols,
             string scalar Gname, real scalar exact)
{
  real matrix D, A, G, T, out
  real rowvector cm, idx, sd, tobs, ord
  real colvector ok, mx
  real scalar n, H, h, j, nex, pprev, pj, bad
  D = st_matrix(Dname)
  cm = st_matrix(CMname)
  idx = strtoreal(tokens(cols))
  A = D[., idx']
  cm = cm[1, idx']
  ok = rowmissing(A) :== 0
  n = sum(ok)
  A = select(A, ok)
  printf("  listwise-complete replications over the family: %g\n", n)
  if (exact) sd = sqrt(diagonal(variance(A)))' :* sqrt((rows(A)-1)/rows(A))
  else sd = sqrt(diagonal(variance(A \ cm)))'
  T = abs(A) :/ sd
  tobs = abs(cm) :/ sd
  H = cols(T)
  out = J(4, H, .)
  out[1,.] = tobs
  for (h=1;h<=H;h++) {
    nex = sum(T[.,h] :>= tobs[h])
    out[2,h] = (exact ? nex/n : (1+nex)/(n+1))
  }
  mx = rowmax(T)
  for (h=1;h<=H;h++) {
    nex = sum(mx :>= tobs[h])
    out[3,h] = (exact ? nex/n : (1+nex)/(n+1))
  }
  ord = order(-tobs', 1)'
  pprev = 0
  for (j=1;j<=H;j++) {
    mx = rowmax(T[., ord[j..H]])
    nex = sum(mx :>= tobs[ord[j]])
    pj = (exact ? nex/n : (1+nex)/(n+1))
    if (pj < pprev) pj = pprev
    out[4, ord[j]] = pj
    pprev = pj
  }
  G = st_matrix(Gname)
  bad = 0
  for (h=1;h<=H;h++) {
    for (j=1;j<=4;j++) {
      if (reldif(G[j,h], out[j,h]) > 1e-12) {
        printf("  SD-MISMATCH row %g col %g: got %21x want %21x\n", j, h, G[j,h], out[j,h])
        bad = 1
      }
    }
  }
  // invariants
  ord = order(-G[1,.]', 1)'
  pprev = -1
  for (j=1;j<=H;j++) {
    if (G[4,ord[j]] < pprev - 1e-15) {
      printf("  SD-NONMONOTONE at ordered position %g\n", j)
      bad = 1
    }
    pprev = G[4,ord[j]]
  }
  for (h=1;h<=H;h++) {
    if (G[4,h] < G[2,h] - 1e-15) {
      printf("  SD-BELOW-RAW col %g: stepdown %g < raw %g\n", h, G[4,h], G[2,h])
      bad = 1
    }
    if (G[4,h] > G[3,h] + 1e-15) {
      printf("  SD-ABOVE-SINGLESTEP col %g: stepdown %g > singlestep %g\n", h, G[4,h], G[3,h])
      bad = 1
    }
  }
  printf("  sdcheck: %g problem(s)\n", bad)
}
end
include "dgp.do"
tempfile base
save `base'
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(3) espostlength(3) jointtests permute(shuffle) permreps(199) rseed(11) permstepdown permsave(s_def) replace
dis "RC_SD=" _rc
matrix _PD = e(permdist_b)
matrix _CM = e(perm_comparator_b)
matrix SD  = e(perm_stepdown)
cap matrix SDPRE = e(perm_stepdown_pre)
cap confirm matrix SDPRE
dis "SDPRE exists rc=" _rc
matrix list SD, format(%12.8f)
cap matrix list SDPRE, format(%12.8f)
matrix list e(perm_n)
matrix list e(perm_p), format(%12.8f)
mata:
  nm = st_matrixcolstripe("_PD")
  for (j=1;j<=rows(nm);j++) printf("bcol %g = %s\n", j, nm[j,2])
  sn = st_matrixcolstripe("SD")
  for (j=1;j<=rows(sn);j++) printf("sdcol %g = %s\n", j, sn[j,2])
end
dis "PROBE04_DONE"
