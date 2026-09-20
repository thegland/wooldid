mata:
real scalar wooldid_permq(real colvector v, real scalar pr)
{
  real scalar n, x
  n = rows(v)
  x = n*pr
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1, n))])/2)
  return(v[max((min((ceil(x), n)), 1))])
}
real scalar wooldid_permmemb(real matrix Am, real rowvector cm, real scalar M,
                             real scalar orow)
{
  if (M >= . | orow < 1 | orow > rows(Am)) return(0)
  if (hasmissing(cm) | hasmissing(Am[orow, .])) return(0)
  return(Am[orow, .] == cm)
}
real rowvector wooldid_permtail(real colvector draws, real scalar sobs, real scalar exact,
                                real scalar memb, real scalar al, real scalar zc,
                                real scalar absp)
{
  // absp is the CONVENTION CODE (F052/F064, 2026-09-06): 0 equal-tailed, 1 absolute add-one
  // (permabsp), 2 raw c/n |s| share with no add-one (permrawp), 3 left one-sided, 4 right
  // one-sided (permside()). Rows 14-18 (appended; rows 1-13 keep their meaning): 14 P_lo
  // and 15 P_hi (the one-sided shares under the add-one / plain-share rule), 16 se_p (the
  // Monte-Carlo standard error of the reported p's underlying share, doubled under the
  // equal-tailed 2*min rule), 17-18 the Clopper-Pearson interval of that share at level()
  // (doubled and capped at 1 under equal-tailed)
  real colvector v
  real scalar n, qlo, qhi, m0, m1, plo, phi, qa, ma, cc, qq, al2, cpl, cph, mult
  real rowvector out
  out = J(1, 18, .)
  out[13] = 0
  v = select(draws, draws :< .)
  n = rows(v)
  out[1] = n
  if (n < 2 | sobs >= .) return(out)
  v = sort(v, 1)
  out[2] = sobs
  out[6] = (exact ? sqrt(variance(v)*(n-1)/n) : sqrt(variance(v)))
  qlo = wooldid_permq(v, al)
  qhi = wooldid_permq(v, 1 - al)
  out[11] = qlo
  out[12] = qhi
  m0 = sum(v :<= sobs)
  m1 = sum(v :>= sobs)
  ma = sum(abs(v) :>= abs(sobs))
  // BS-23 (author ruling 2026-09-07): `exact' is COMPLETENESS and still governs the
  // population-SD convention below. `memb' is MEMBERSHIP - the observed statistic is one
  // of the n reference values - which is what the SHARE convention turns on. The add-one
  // form exists to inject an observed statistic the reference set does not contain; where
  // the set already contains it, adding it again counts it twice and doubles the floor.
  // Callers set memb only under an enumeration, and only after checking the observed row's
  // value equals sobs.
  if (exact | memb) {
    plo = m0/n
    phi = m1/n
  }
  else {
    plo = (1 + m0)/(n + 1)
    phi = (1 + m1)/(n + 1)
  }
  out[14] = plo
  out[15] = phi
  mult = 1
  if (absp == 1 | absp == 2) {
    // absolute-value conventions: the level/100 quantile of |s_r| on each side
    qa = wooldid_permq(sort(abs(v), 1), 1 - 2*al)
    if (absp == 1) out[3] = (exact | memb ? ma/n : (1 + ma)/(n + 1))
    else out[3] = ma/n
    out[4] = -qa
    out[5] = qa
    cc = ma
  }
  else {
    out[4] = qlo
    out[5] = qhi
    if (absp == 3) {
      out[3] = plo
      cc = m0
    }
    else if (absp == 4) {
      out[3] = phi
      cc = m1
    }
    else {
      out[3] = min((1, 2*min((plo, phi))))
      cc = min((m0, m1))
      mult = 2
    }
  }
  // Monte-Carlo precision of the reported p: the standard error and the level()
  // Clopper-Pearson interval of the underlying share cc/n (Stata's permute and ritest
  // report both); under the equal-tailed 2*min rule both are doubled, the upper end
  // capped at 1
  // under an enumeration the shares are exact, not Monte-Carlo estimates, so rows 16-18
  // keep the missing they were initialised to (bug sweep BS-22, extended to the
  // incomplete-but-member case by BS-23, author ruling 2026-09-07: an enumeration that
  // lost assignments is still deterministic, so its shares carry no sampling error either.
  // The only element under permexact that still gets a binomial SE is one whose observed
  // row is genuinely absent - the same element that takes the add-one form.)
  if (!exact & !memb) {
    qq = cc/n
    out[16] = mult*sqrt(qq*(1 - qq)/n)
    al2 = al
    cpl = (cc <= 0 ? 0 : invibeta(cc, n - cc + 1, al2))
    cph = (cc >= n ? 1 : invibeta(cc + 1, n - cc, 1 - al2))
    out[17] = min((1, mult*cpl))
    out[18] = min((1, mult*cph))
  }
  // the shift interval takes the estimate's sampling error to have the NULL permutation
  // distribution, which presumes that distribution is centred at zero. When its central
  // band excludes 0 the presumption fails (a quantity that is non-zero under reassigned
  // onsets - e.g. an ame row carrying a baseline dose slope - or a probabilities design
  // far from the observed assignment) and the interval would not even contain the
  // estimate: it is left blank and flagged (row 10); the rank p-value stands
  if (qlo > 0 | qhi < 0) out[10] = 1
  else {
    out[7] = sobs - out[5]
    out[8] = sobs - out[4]
    // row 9 is the standard error IMPLIED BY INVERTING the reported shift interval,
    // (q_hi - q_lo)/(2*z(level)), on the statistic's scale. It is NOT the dispersion of the
    // null draws (that is row 6, e(perm_sd)) and not an analytic SE. It is published under
    // its own name, e(perm_se_shift), and deliberately not posted into e(V) or e(se)
    // (author rulings 2026-09-03/04, revisited under PS-2 2026-09-10; the bootstrap side
    // publishes no interval-implied SE at all). It inherits the off-centre guard: when row
    // 10 fires, rows 7, 8 and 9 stay missing together, because an SE implied by an interval
    // that was refused is meaningless - the SE column vanishing there is the design, not a
    // bug
    if (zc < . & zc > 0) out[9] = (out[5] - out[4])/(2*zc)
  }
  out[13] = 1
  return(out)
}

end
