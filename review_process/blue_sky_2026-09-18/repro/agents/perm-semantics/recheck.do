*! recheck.do -- INDEPENDENT recomputation of every permute() reported quantity
*! from the saved null distribution, following the HELP FILE's definitions only.
*! Expects (copied out of e()): _PD (draws R x K), _CM (comparator 1 x K),
*! _OB (e(b) 1 x K), _PK (e(perm_kernel) 18 x K); optional _PDS / _CMS under permtstat.
mata:
// documented quantile, re-implemented independently: x = n*pr; if x is a positive
// integer take the mean of order statistics x and x+1, else order statistic ceil(x)
real scalar rq(real colvector sv, real scalar pr)
{
  real scalar n, x
  n = rows(sv)
  x = n*pr
  if (x == floor(x) & x >= 1) return((sv[x] + sv[min((x+1,n))])/2)
  return(sv[max((min((ceil(x),n)),1))])
}
void recheck(string scalar tag, real scalar lvl, real scalar absp, real scalar exact)
{
  real matrix D, K, S, cs
  real rowvector cm, ob
  real colvector v, s
  real scalar R, KK, k, n, sobs, m0, m1, ma, plo, phi, p, sd, qlo, qhi, qa
  real scalar al, zc, lb, ub, sesh, cc, mult, qq, pse, cpl, cph, oc, bad, nbad
  real scalar qU, qL, ts, r
  D = st_matrix("_PD")
  cm = st_matrix("_CM")
  ob = st_matrix("_OB")
  K  = st_matrix("_PK")
  S  = st_matrix("_PDS")
  cs = st_matrix("_CMS")
  ts = (rows(S) > 0)
  R = rows(D)
  KK = cols(D)
  al = (100 - lvl)/200
  zc = invnormal(1 - al)
  nbad = 0
  printf("\n==== RECHECK %s (lvl=%g absp=%g exact=%g R=%g K=%g ts=%g) ====\n", tag, lvl, absp, exact, R, KK, ts)
  for (k = 1; k <= KK; k++) {
    if (K[13,k] != 1) continue
    if (ts == 0) {
      s = D[.,k]
      sobs = cm[1,k]
    }
    else {
      s = J(R,1,.)
      for (r=1; r<=R; r++) {
        if (D[r,k]<. & S[r,k]<. & S[r,k]>0) s[r] = D[r,k]/S[r,k]
      }
      if (cs[1,k]>=. | cs[1,k]<=0 | cm[1,k]>=.) sobs = .
      else sobs = cm[1,k]/cs[1,k]
    }
    v = select(s, s :< .)
    n = rows(v)
    if (n < 2 | sobs >= .) continue
    v = sort(v,1)
    m0 = sum(v :<= sobs)
    m1 = sum(v :>= sobs)
    ma = sum(abs(v) :>= abs(sobs))
    if (exact) {
      plo = m0/n
      phi = m1/n
    }
    else {
      plo = (1+m0)/(n+1)
      phi = (1+m1)/(n+1)
    }
    mult = 1
    if (absp == 1) {
      p = (exact ? ma/n : (1+ma)/(n+1))
      cc = ma
    }
    else if (absp == 2) {
      p = ma/n
      cc = ma
    }
    else if (absp == 3) {
      p = plo
      cc = m0
    }
    else if (absp == 4) {
      p = phi
      cc = m1
    }
    else {
      p = min((1, 2*min((plo,phi))))
      cc = min((m0,m1))
      mult = 2
    }
    sd  = exact ? sqrt(variance(v)*(n-1)/n) : sqrt(variance(v))
    qlo = rq(v, al)
    qhi = rq(v, 1-al)
    qa  = rq(sort(abs(v),1), lvl/100)
    oc  = (qlo > 0 | qhi < 0)
    qU = (absp==1 | absp==2) ?  qa : qhi
    qL = (absp==1 | absp==2) ? -qa : qlo
    if (oc) {
      sesh = .
      lb = .
      ub = .
    }
    else {
      sesh = (qU - qL)/(2*zc)
      if (ts == 0) {
        lb = ob[1,k] - qU
        ub = ob[1,k] - qL
      }
      else {
        lb = ob[1,k] - qU*cs[1,k]
        ub = ob[1,k] - qL*cs[1,k]
        sesh = sesh*cs[1,k]
      }
    }
    qq = cc/n
    pse = .
    cpl = .
    cph = .
    if (exact == 0) {
      pse = mult*sqrt(qq*(1-qq)/n)
      cpl = (cc <= 0 ? 0 : invibeta(cc, n-cc+1, al))
      cph = (cc >= n ? 1 : invibeta(cc+1, n-cc, 1-al))
      cpl = min((1, mult*cpl))
      cph = min((1, mult*cph))
    }
    bad = 0
    if (K[1,k]  != n)   {
      printf("  MISMATCH col %g n: got %g want %g\n", k, K[1,k], n)
      bad=1
    }
    if (K[3,k]  != p)   {
      printf("  MISMATCH col %g p: got %21x want %21x\n", k, K[3,k], p)
      bad=1
    }
    if (K[14,k] != plo) {
      printf("  MISMATCH col %g p_lo: got %21x want %21x\n", k, K[14,k], plo)
      bad=1
    }
    if (K[15,k] != phi) {
      printf("  MISMATCH col %g p_hi: got %21x want %21x\n", k, K[15,k], phi)
      bad=1
    }
    if (reldif(K[6,k], sd)  > 1e-12) {
      printf("  MISMATCH col %g sd: got %21x want %21x\n", k, K[6,k], sd)
      bad=1
    }
    if (reldif(K[11,k], qlo) > 1e-12) {
      printf("  MISMATCH col %g bandlo: got %21x want %21x\n", k, K[11,k], qlo)
      bad=1
    }
    if (reldif(K[12,k], qhi) > 1e-12) {
      printf("  MISMATCH col %g bandhi: got %21x want %21x\n", k, K[12,k], qhi)
      bad=1
    }
    if ((K[10,k]==1) != (oc==1)) {
      printf("  MISMATCH col %g offcentre flag: got %g want %g\n", k, K[10,k], oc)
      bad=1
    }
    if (oc == 0) {
      if (reldif(K[7,k], lb) > 1e-12) {
        printf("  MISMATCH col %g ci_lb: got %21x want %21x\n", k, K[7,k], lb)
        bad=1
      }
      if (reldif(K[8,k], ub) > 1e-12) {
        printf("  MISMATCH col %g ci_ub: got %21x want %21x\n", k, K[8,k], ub)
        bad=1
      }
      if (reldif(K[9,k], sesh) > 1e-12) {
        printf("  MISMATCH col %g se_shift: got %21x want %21x\n", k, K[9,k], sesh)
        bad=1
      }
    }
    else {
      if (K[7,k]<. | K[8,k]<. | K[9,k]<.) {
        printf("  MISMATCH col %g offcentre but ci/se not blank\n", k)
        bad=1
      }
    }
    if (reldif(K[16,k], pse) > 1e-10 & !(K[16,k]>=. & pse>=.)) {
      printf("  MISMATCH col %g p_se: got %21x want %21x\n", k, K[16,k], pse)
      bad=1
    }
    if (reldif(K[17,k], cpl) > 1e-9 & !(K[17,k]>=. & cpl>=.)) {
      printf("  MISMATCH col %g cp_lo: got %21x want %21x\n", k, K[17,k], cpl)
      bad=1
    }
    if (reldif(K[18,k], cph) > 1e-9 & !(K[18,k]>=. & cph>=.)) {
      printf("  MISMATCH col %g cp_hi: got %21x want %21x\n", k, K[18,k], cph)
      bad=1
    }
    if (p == 0) printf("  ZEROP col %g p=0 absp=%g exact=%g n=%g\n", k, absp, exact, n)
    nbad = nbad + bad
  }
  printf("==== RECHECK %s : %g column(s) with a mismatch ====\n", tag, nbad)
}
end
