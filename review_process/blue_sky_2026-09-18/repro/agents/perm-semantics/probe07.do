* probe07 -- the QUADRATIC-FORM tests and the PAIRWISE CONTRASTS.
* (i) every joint / equality / adversarial p must equal (1+#{W_r >= W_obs})/(n+1)
*     from the stored replication statistics and their observed twins;
* (ii) the joint pre-trend W must be the documented permutation-covariance quadratic form,
*     recomputed independently from e(permdist_b);
* (iii) every pairwise subgroup contrast kernel must be the permutation kernel of the
*     coefficient DIFFERENCE, recomputed from the dump.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
mata:
// (i) p from (W_r, W_obs)
void wpcheck(string scalar Wname, string scalar Wobsname, string scalar label)
{
  real matrix W, Wo
  real colvector w
  real scalar j, n, nex, p
  string matrix nm
  W = st_matrix(Wname)
  Wo = st_matrix(Wobsname)
  nm = st_matrixcolstripe(Wname)
  for (j=1;j<=cols(W);j++) {
    w = select(W[.,j], W[.,j] :< .)
    n = rows(w)
    nex = sum(w :>= Wo[1,j])
    p = (1+nex)/(n+1)
    printf("  %s col %g %-24s  n=%g  Wobs=%g  recomputed p=%21x  (%12.8f)\n", label, j, nm[j,2], n, Wo[1,j], p, p)
    st_numscalar("WP_" + strofreal(j), p)
  }
}
// (ii) the joint pre-trend quadratic form, recomputed from the draws
void jwcheck(string scalar cols, real scalar orow)
{
  real matrix D, A, CC, dev, Sig, Si, R
  real rowvector cm, idx, sc
  real colvector ok, Wr
  real scalar n, r, Wobs, memb, worst
  real matrix Wst
  D = st_matrix("_PD")
  cm = st_matrix("_CM")
  idx = strtoreal(tokens(cols))
  A = D[., idx']
  cm = cm[1, idx']
  ok = rowmissing(A) :== 0
  n = sum(ok)
  memb = 0
  CC = select(A, ok) \ cm
  dev = CC :- mean(CC)
  Sig = (dev' * dev)/(rows(CC) - 1)
  sc = sqrt(diagonal(Sig))'
  R = Sig :/ (sc' * sc)
  Si = invsym(R) :/ (sc' * sc)
  Wobs = cm * Si * cm'
  Wr = J(rows(A), 1, .)
  for (r=1;r<=rows(A);r++) {
    if (ok[r]) Wr[r] = A[r,.] * Si * A[r,.]'
  }
  Wst = st_matrix("JW")
  worst = 0
  for (r=1;r<=rows(A);r++) {
    if (Wr[r]<. & Wst[r,1]<.) {
      if (reldif(Wr[r], Wst[r,1]) > worst) worst = reldif(Wr[r], Wst[r,1])
    }
  }
  printf("  jointW recomputed: n=%g  Wobs=%21x  stored Wobs=%21x  reldif=%g\n",
         n, Wobs, st_matrix("JWO")[1,1], reldif(Wobs, st_matrix("JWO")[1,1]))
  printf("  jointW draws worst reldif over %g replications = %g\n", rows(A), worst)
  printf("  jointW recomputed p = %12.8f\n", (1+sum(select(Wr, Wr:<.) :>= Wobs))/(n+1))
}
// (iii) pairwise contrasts: identify the coefficient pair by its comparator difference
void pwcheck(real scalar lvl, real scalar absp)
{
  real matrix D, PK, A
  real rowvector cm, row
  real colvector d, ok, v
  real scalar j, ka, kb, a, b, sobs, found, n, m0, m1, ma, plo, phi, p, bad
  string matrix nm
  D = st_matrix("_PD")
  cm = st_matrix("_CM")
  PK = st_matrix("PWK")
  nm = st_matrixcolstripe("PWK")
  bad = 0
  for (j=1;j<=cols(PK);j++) {
    if (PK[13,j] != 1) continue
    sobs = PK[2,j]
    found = 0
    ka = 0
    kb = 0
    for (a=1;a<=cols(D);a++) {
      for (b=1;b<=cols(D);b++) {
        if (a != b & cm[1,a]<. & cm[1,b]<.) {
          if (cm[1,a] - cm[1,b] == sobs) {
            if (!found) {
              ka = a
              kb = b
            }
            found = found + 1
          }
        }
      }
    }
    if (ka == 0) {
      printf("  PW col %g %s: could not identify the coefficient pair\n", j, nm[j,2])
      continue
    }
    A = D[., (ka, kb)]
    ok = rowmissing(A) :== 0
    d = select(A[.,1] :- A[.,2], ok)
    n = rows(d)
    v = sort(d, 1)
    m0 = sum(v :<= sobs)
    m1 = sum(v :>= sobs)
    ma = sum(abs(v) :>= abs(sobs))
    plo = (1+m0)/(n+1)
    phi = (1+m1)/(n+1)
    if (absp == 1) p = (1+ma)/(n+1)
    else if (absp == 2) p = ma/n
    else if (absp == 3) p = plo
    else if (absp == 4) p = phi
    else p = min((1, 2*min((plo,phi))))
    printf("  PW col %g %-22s cols(%g,%g) matches=%g n: got %g want %g | p: got %21x want %21x\n",
           j, nm[j,2], ka, kb, found, PK[1,j], n, PK[3,j], p)
    if (PK[1,j] != n | PK[3,j] != p) {
      printf("    PW-MISMATCH\n")
      bad = bad + 1
    }
    if (reldif(PK[11,j], rq(v, (100-lvl)/200)) > 1e-12) {
      printf("    PW-MISMATCH bandlo got %21x want %21x\n", PK[11,j], rq(v,(100-lvl)/200))
      bad = bad + 1
    }
  }
  printf("  pwcheck: %g problem(s)\n", bad)
}
end
include "dgp.do"
tempfile base
save `base'
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(3) espostlength(3) subgroups(sg) jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(99) rseed(13) permsave(g_sg) replace
dis "RC_SG=" _rc
local htag SG
local habsp 0
local hexact 0
include "harvest.do"
mata:
  nm = st_matrixcolstripe("_PD")
  for (j=1;j<=rows(nm);j++) printf("bcol %g = %s\n", j, nm[j,2])
end
* the stored test statistics
cap confirm matrix e(perm_jointW)
dis "jointW rc=" _rc
matrix JW = e(perm_jointW)
matrix JWO = e(perm_jointWobs)
mata: wpcheck("JW", "JWO", "JOINT")
dis "e(joint_horiz)=" e(joint_horiz)
dis "e(joint_horiz_s0)=" e(joint_horiz_s0) "  e(joint_horiz_s1)=" e(joint_horiz_s1)
cap confirm matrix e(perm_sgW)
dis "sgW rc=" _rc
if _rc == 0 {
  matrix SW = e(perm_sgW)
  matrix SWO = e(perm_sgWobs)
  mata: wpcheck("SW", "SWO", "EQUALITY")
}
cap confirm matrix e(perm_advW)
dis "advW rc=" _rc
if _rc == 0 {
  matrix AW = e(perm_advW)
  matrix AWO = e(perm_advWobs)
  mata: wpcheck("AW", "AWO", "ADV")
}
ereturn list
* pairwise contrasts
cap confirm matrix e(perm_pwsgs_kernel)
dis "pwkernel rc=" _rc
if _rc == 0 {
  matrix PWK = e(perm_pwsgs_kernel)
  mata: pwcheck(95, 0)
}
dis "PROBE07_DONE"
