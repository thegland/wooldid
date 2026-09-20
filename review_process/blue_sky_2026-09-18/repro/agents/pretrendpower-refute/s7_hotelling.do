* s7_hotelling.do -- SELF-CONTAINED. No panel needed for the core demonstration.
*
* CLAIM: pretrendpower() applies KNOWN-covariance formulas to an ESTIMATED covariance.
*   ado:10516  local ptpprep = Ftail(keff, dfr, W/keff)        <- W = b'S^-1 b
*   ado:10559  local ptplam  = npnF(keff, dfr, invF(keff,dfr,.95), 1-power)
*   ado:10512  local ptpse   = 1/sqrt(m'S^-1 m)  ; CI uses invttail(dfr,.025)
* Those are exact only when S is the TRUE covariance.  When S is estimated with dfr
* degrees of freedom the exact reference for b'S^-1 b is Hotelling's T-squared:
*   W*(dfr-k+1)/(dfr*k) ~ F(k, dfr-k+1).
* PART A  computes the exact size of wooldid's own cutoff for a grid of (k, dfr).
* PART B  simulates the IDEAL model -- b ~ N(0,I), S ~ Wishart(I,dfr)/dfr, independent,
*         no panel, no clustering, no DiD -- and measures the size of wooldid's test,
*         the coverage of its slope CI, and the sampling distribution of its mdt column.
* PART C  runs wooldid once on a 16-state panel and shows dfr = 15, keff = 3, so the
*         PART A row k=3 dfr=15 is the row that design lands on.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

dis "================ PART A: exact size of wooldid's joint pre-trend cutoff ================"
dis "   k   dfr    wooldid W cutoff   nominal size   EXACT size   mdt_reported/mdt_needed(80%)"
foreach kk in 1 2 3 4 {
  foreach dd in 15 23 39 79 279 {
    local cut = `kk'*invF(`kk', `dd', .95)
    * exact Hotelling size of the rule  W > cut
    local sz = Ftail(`kk', `dd'-`kk'+1, `cut'*(`dd'-`kk'+1)/(`dd'*`kk'))
    * wooldid's noncentrality for 80% power, and the noncentrality that really delivers
    * 80% power against wooldid's OWN cutoff under the exact reference
    local lw = npnF(`kk', `dd', invF(`kk',`dd',.95), .20)
    local lt = npnF(`kk', `dd'-`kk'+1, `cut'*(`dd'-`kk'+1)/(`dd'*`kk'), .20)
    dis %4.0f `kk' %6.0f `dd' %18.4f `cut' "          0.0500" %13.4f `sz' %18.4f sqrt(`lw'/`lt')
  }
}

dis ""
dis "================ PART B: the IDEAL model, 100,000 draws ================"
mata:
mata clear
void idealrun(real scalar k, real scalar dfr, real scalar R)
{
  real matrix m, S, Si, L, Z
  real scalar r, W, c2, sl, se, cut, crit, nrej, ncov, i
  real colvector mdt, slv, sev
  m = (-3 \ -2 \ -1)
  if (k != 3) m = -range(k,1,1)
  cut  = k*invF(k, dfr, .95)
  crit = invttail(dfr, .025)
  mdt = J(R,1,.); slv = J(R,1,.); sev = J(R,1,.)
  nrej = 0; ncov = 0
  for (r = 1; r <= R; r++) {
    // S ~ Wishart(I, dfr)/dfr  via  S = Z'Z/dfr with Z (dfr x k) standard normal
    Z = rnormal(dfr, k, 0, 1)
    S = cross(Z, Z)/dfr
    Si = invsym(S)
    // b ~ N(0, I), independent of S
    L = rnormal(k, 1, 0, 1)
    W  = (L' * Si * L)
    c2 = (m' * Si * m)
    sl = (m' * Si * L)/c2
    se = 1/sqrt(c2)
    if (W > cut) nrej++
    if (abs(sl) <= crit*se) ncov++
    mdt[r] = 1/sqrt(c2)
    slv[r] = sl
    sev[r] = se
  }
  printf("k=%f dfr=%f:  size of wooldid's joint test = %9.4f   (nominal 0.0500)\n", k, dfr, nrej/R)
  printf("            GLS slope %s CI coverage       = %9.4f   (nominal 0.9500)\n", "95%", ncov/R)
  printf("            mean reported se / MC sd of slope = %7.4f\n", mean(sev)/sqrt(variance(slv)))
  // mdt_true: with S == I the truth is 1/sqrt(m'm^-1...) -- compute m'I^-1 m directly
  printf("            median reported (1/sqrt(m'S^-1 m)) / true (1/sqrt(m'm)) = %7.4f  (theory sqrt(chi2(dfr-k+1)/dfr) median)\n",
         (sort(mdt,1)[round(R/2)]) / (1/sqrt(m'*m)))
}
end
mata: idealrun(3, 15, 100000)
mata: idealrun(3, 79, 100000)
mata: idealrun(3, 279, 100000)
mata: idealrun(4, 15, 100000)

dis ""
dis "================ PART C: which row a 16-state DiD lands on ================"
clear
set obs 16
gen long st = _n
gen byte coh = 1 + floor((st-1)/4)
gen int ttre = .
replace ttre = 6  if coh == 1
replace ttre = 8  if coh == 2
replace ttre = 10 if coh == 3
gen double alpha_s = rnormal()*0.8
expand 12
bysort st: gen int t = _n
sort st t
gen double e = rnormal()
forvalues k = 2/12 {
  qui replace e = 0.5*e[_n-1] + sqrt(.75)*rnormal() if t == `k'
}
gen double y = 2 + alpha_s + 0.12*t - 0.02*t*t + 0.8*(ttre<. & t>=ttre) + e
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "PARTC_RC=" _rc
dis "PARTC e(df_r) = " e(df_r) "   e(roth_keff_att) = " e(roth_keff_att) "   e(N_clust) = " e(N_clust)
dis "PARTC -> the row k=3, dfr=15 of PART A: exact size " %6.4f Ftail(3,13,3*invF(3,15,.95)*13/45) " for a test the run labels 5%"
dis "DONE"
