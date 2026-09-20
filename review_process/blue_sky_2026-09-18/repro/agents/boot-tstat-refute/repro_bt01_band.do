*=====================================================================
* BT-01 / BT-03 VERDICT REPRO (skeptic's own)
* A: boottstat + simultaneousbands.  Reports, per horizon and SEPARATELY,
*    (i) whether a band endpoint lies INSIDE the pointwise percentile-t interval
*        (non-nesting -- the behaviour sthlp documents for percentile/bc/bca) and
*    (ii) whether the band is STRICTLY NARROWER in width (what sthlp:528's
*        boottstat exemption actually denies).
*    Mechanism: wooldid_boottp (ado:14332) makes an EQUAL-TAILED interval from the
*    signed t-quantiles; wooldid_bootband (ado:16776) makes a SYMMETRIC one from
*    q_lev(max_h |t|).  Different shapes, so they can cross.
* B: bootdisplay(wald) at G=10 -- the configuration BT-03 says must fail.  It does
*    not: the band's critical value comes out LARGER than invt(G-1,.975).
* Look for: A_SUM ... narrower=1   and   B_SUM ... inside_lo=0 inside_hi=0
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop bandcmp
prog def bandcmp
  syntax , tag(string)
  tempname W BD
  cap mat `BD' = e(esband_att)
  if _rc != 0 {
    dis "`tag'_NOBAND=1"
    exit
  }
  mat `W' = e(wdidmainresults)
  local rfn : rowfullnames `W'
  local cn : colnames `BD'
  local m = colsof(`BD')
  local nlo 0
  local nhi 0
  local nnar 0
  forval j = 1/`m' {
    local h : word `j' of `cn'
    local rr 0
    local k 0
    foreach r of local rfn {
      local ++k
      if "`r'" == "ES_att:`h'" local rr = `k'
    }
    if `rr' == 0 continue
    local plo = `W'[`rr',5]
    local phi = `W'[`rr',6]
    local blo = `BD'[1,`j']
    local bhi = `BD'[2,`j']
    local vlo = (`blo' > `plo' + 1e-12)
    local vhi = (`bhi' < `phi' - 1e-12)
    local nar = ((`bhi'-`blo') < (`phi'-`plo') - 1e-12)
    local nlo = `nlo' + `vlo'
    local nhi = `nhi' + `vhi'
    local nnar = `nnar' + `nar'
    dis "`tag'_H `h'  PW[" %13.7f `plo' "," %13.7f `phi' "] w=" %13.7f `phi'-`plo' ///
        "  BAND[" %13.7f `blo' "," %13.7f `bhi' "] w=" %13.7f `bhi'-`blo' ///
        "  inside_lo=" `vlo' " inside_hi=" `vhi' " narrower=" `nar'
  }
  dis "`tag'_SUM m=" `m' " inside_lo=" `nlo' " inside_hi=" `nhi' " narrower=" `nnar' ///
      "  crit=" %14.6g e(esband_crit_att) "  nband=" e(esband_n_att) "  bootdf=" e(boot_df) ///
      "  engine=`e(esband_engine)'"
end

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

dis "### A: boottstat + simultaneousbands, log-normal errors, G=20, 5 horizons"
set seed 8571
set sortseed 8571
mkpanel , g(20) t(9) skew(1)
qui count
dis "A_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(299) boottstat simultaneousbands dots(0)
dis "A_RC=" _rc
bandcmp , tag(A)

dis "### B: bootdisplay(wald) + simultaneousbands, G=10, 4 horizons (BT-03's claim)"
set seed 777333
set sortseed 777333
mkpanel , g(10) t(9) skew(0)
qui count
dis "B_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   bootstrap(standard) reps(999) bootdisplay(wald) simultaneousbands dots(0)
dis "B_RC=" _rc
dis "B_invt=" %14.9f invt(e(boot_df),0.975) "  B_bandcrit=" %14.9f e(esband_crit_att)
bandcmp , tag(B)
dis "### C: bootdisplay(wald) + simultaneousbands, G=8, only TWO horizons"
dis "###    (the family size, not the cluster count, is what makes BT-03 fire)"
set seed 30997
set sortseed 30997
mkpanel , g(8) t(9) skew(0)
qui count
dis "C_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(1) espostlength(1) ///
   bootstrap(standard) reps(999) bootdisplay(wald) simultaneousbands dots(0)
dis "C_RC=" _rc
dis "C_invt=" %14.9f invt(e(boot_df),0.975) "  C_bandcrit=" %14.9f e(esband_crit_att)
bandcmp , tag(C)
dis "BT01_DONE"
