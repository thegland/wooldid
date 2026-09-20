*==========================================================================
* REPRO: simultaneousbands produces a SIMULTANEOUS band that is NARROWER
*        than the POINTWISE interval printed beside it, at the same horizon
*        and the same level().
*
*  Run:  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_band_narrower.do
*
*  Setting: 20 states, 3 units each, 8 years, staggered onsets, a skewed
*  (log-normal) outcome error -- an ordinary shape for earnings, spending or
*  cost outcomes.  cluster(state), bootstrap(standard), boottstat,
*  simultaneousbands.
*
*  Mechanism: the pointwise percentile-t interval is EQUAL-TAILED,
*      [b - q_.975(t_r)*se , b - q_.025(t_r)*se],
*  while the simultaneous band is SYMMETRIC,
*      b +/- q_.95(max_h |t_rh|)*se     (wooldid_bootband, ado:16876).
*  For a symmetric t_r distribution q_.95(|t|) = q_.975(t), so the band
*  dominates the pointwise interval.  When the replication t-distribution is
*  SKEWED the two shapes cross, and the "simultaneous" band excludes values
*  that the pointwise interval at the same horizon includes.
*
*  Expected: a 95% simultaneous band must be at least as wide as the 95%
*            pointwise interval at every horizon of the family - that is what
*            makes it simultaneous.
*  Observed: at horizon pre3 the band's lower endpoint sits INSIDE the
*            pointwise lower endpoint by about 2.7 units.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

clear
set obs 20
gen int state = _n
gen int onset = .
forval s = 1/15 {
  local w = mod(`s'-1, 3) + 1
  local o : word `w' of 2003 2004 2005
  qui replace onset = `o' in `s'
}
gen double a_i = rnormal()*1.0
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand 8
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1.0 + 0.2*(t - onset)
qui replace tau = 0 if !post
gen double eps = (exp(rnormal()*1.2) - exp(1.2^2/2))*0.7
gen double y = a_i + u_i + tshock + tau*post + eps
gen int ttre = onset

wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) ///
        rseed(555) boottstat att esprelength(3) espostlength(3) ///
        simultaneousbands dots(0)

dis _n "esband_engine = " "`e(esband_engine)'"
dis     "critical value c = " e(esband_crit_att)

mat MR = e(wdidmainresults)
mat SB = e(esband_att)
mata:
  MR = st_matrix("MR"); rs = st_matrixrowstripe("MR")
  SB = st_matrix("SB"); cs = st_matrixcolstripe("SB")
  nv = 0
  printf("\n%-8s %14s %14s %14s %14s   %s\n", "horizon", "point_lb", "point_ub", "band_lb", "band_ub", "verdict")
  for (j = 1; j <= cols(SB); j++) {
    hit = 0
    for (r = 1; r <= rows(rs); r++) if (rs[r,1] == "ES_att" & rs[r,2] == cs[j,2]) hit = r
    if (hit == 0) continue
    plb = MR[hit,5]; pub = MR[hit,6]; blb = SB[1,j]; bub = SB[2,j]
    bad = (blb > plb + 1e-9) | (bub < pub - 1e-9)
    nv = nv + bad
    printf("%-8s %14.8f %14.8f %14.8f %14.8f   %s\n", cs[j,2], plb, pub, blb, bub,
           (bad ? "BAND NARROWER THAN POINTWISE" : "ok"))
  }
  printf("\nHORIZONS_WITH_BAND_INSIDE_POINTWISE = %f\n", nv)
  printf("FAIL = %f\n", (nv > 0))
end
