*=====================================================================
* SK-02 SELF-CONTAINED REPRO
* Claim: under boottstat, wooldid_bootband's studentized critical value
*   q_lev(max_h |b_r - b_obs| / se_r) has no guard against a replication whose
*   analytic standard error is near zero -- ado:16776 drops a replication only when
*   a horizon's se is MISSING or NON-POSITIVE.  At cluster() == i, where the
*   replication analytic SEs are themselves unstable, the critical value is routinely
*   an order of magnitude above the nominal ~2.5 and can reach 1e14, so the printed
*   simultaneous band is vacuous and nothing says so.
* Look for: SK02_CRIT (huge) beside SK02_NOMINAL (~2.1)
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , g(integer) t(integer) skew(integer)
  clear
  set obs `g'
  gen int state = _n
  gen int onset = .
  local c1 = max(1,round(0.25*`g'))
  local c2 = round(0.48*`g')
  local c3 = round(0.70*`g')
  qui replace onset = 2002 if state <= `c1'
  qui replace onset = 2003 if state >  `c1' & state <= `c2'
  qui replace onset = 2004 if state >  `c2' & state <= `c3'
  gen double a_i = rnormal()
  gen int nun = 2 + mod(state,3)
  expand nun
  bys state: gen int unit = _n + 100*state
  gen double u_i = rnormal()*0.6
  expand `t'
  bys unit: gen int tt = 1999 + _n
  rename tt t
  gen double tshock = 0
  qui bys t (unit): replace tshock = rnormal() if _n==1
  qui bys t (unit): replace tshock = tshock[1]
  gen byte post = (t >= onset) & !missing(onset)
  gen double tau = 1 + 0.2*(t-onset)
  qui replace tau = 0 if !post
  if `skew' gen double eps = (exp(rnormal()*1.3) - exp(1.3^2/2))*0.7
  else      gen double eps = rnormal()*0.7
  gen double y = a_i + u_i + tshock + tau*post + eps
  gen int ttre = onset
end

cap prog drop critrep
prog def critrep
  syntax , tag(string)
  tempname BD
  cap mat `BD' = e(esband_att)
  if _rc != 0 {
    dis "`tag'_NOBAND=1 why=`e(esband_why_att)'"
    exit
  }
  dis "`tag'_CRIT=" %16.6g e(esband_crit_att) "  nband=" e(esband_n_att) "  engine=`e(esband_engine)'"
  local m = colsof(`BD')
  local cn : colnames `BD'
  forval j = 1/`m' {
    local h : word `j' of `cn'
    dis "`tag'_BAND `h' [" %16.6g `BD'[1,`j'] "," %16.6g `BD'[2,`j'] "]  halfwidth=" %16.6g (`BD'[2,`j']-`BD'[1,`j'])/2
  }
end

dis "### A: G=15, log-normal errors, 5 horizons, boottstat + simultaneousbands"
set seed 12142
set sortseed 12142
mkpanel , g(15) t(9) skew(1)
qui count
dis "A_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(299) boottstat simultaneousbands dots(0)
dis "A_RC=" _rc
critrep , tag(A)
dis "SK02_NOMINAL=" %12.8f invnormal(0.975) "  (a well-behaved sup over 5 horizons should be ~2.5-3)"

dis "### B: the SAME shape with GAUSSIAN errors, G=20 -- is the inflation only skew?"
set seed 15713
set sortseed 15713
mkpanel , g(20) t(9) skew(0)
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(299) boottstat simultaneousbands dots(0)
dis "B_RC=" _rc
critrep , tag(B)
dis "SK02_DONE"
