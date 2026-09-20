* s3_mc_null.do -- SELF-CONTAINED. Own 16-state panel, own seeds.
* Null Monte Carlo of pretrendpower()'s joint pre-trend test, GLS slope CI and
* mdt_wald(80), PLUS the same 3-restriction cluster-robust joint test run by hand
* through reghdfe/testparm on the identical data (the "is this wooldid or is this
* Stata's standard F reference?" comparator).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local R 250

cap prog drop skpanel
prog def skpanel
  args delta rho
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
  gen double yeareff = 0.12*t - 0.02*t*t + 0.25*cos(t)
  gen byte D = (ttre < . & t >= ttre)
  gen double relyr = t - ttre
  gen double trendterm = cond(ttre < ., relyr, 0)
  gen double e = rnormal()
  forvalues k = 2/12 {
    qui replace e = `rho'*e[_n-1] + sqrt(1-`rho'*`rho')*rnormal() if t == `k'
  }
  gen double y = 2 + alpha_s + yeareff + 0.8*D + `delta'*trendterm + e
end

* the three pre-period coefficients sit at relyr -4,-3,-2 with reference -1, so
local mvec "(-3 \ -2 \ -1)"
* wooldid's own noncentrality for 80% power, k=3, df_r=15 (ado:10559)
local lam80 = npnF(3, 15, invF(3, 15, .95), 1 - 80/100)
dis "lam80 = " `lam80'
local cutW = 3*invF(3, 15, .95)
dis "wooldid's W cutoff = " `cutW'

* --- extract the pre-block statistics from a posted wooldid fit ---
cap prog drop pstat
prog def pstat, rclass
  args lam
  tempname B V S bp Si m w
  matrix `B'  = e(b_es_att)
  matrix `V'  = e(V_es_att)
  matrix `S'  = `V'[1..3, 1..3]
  matrix `S'  = (`S' + `S'')/2
  matrix `bp' = (`B'[1,1] \ `B'[1,2] \ `B'[1,3])
  matrix `Si' = invsym(`S')
  matrix `m'  = (-3 \ -2 \ -1)
  matrix `w'  = `bp''*`Si'*`bp'
  return scalar W  = `w'[1,1]
  matrix `w'  = `m''*`Si'*`m'
  return scalar c2 = `w'[1,1]
  matrix `w'  = `m''*`Si'*`bp'
  return scalar slope = `w'[1,1]/return(c2)
  return scalar se    = 1/sqrt(return(c2))
  return scalar mdt80 = sqrt(`lam'/return(c2))
  return scalar b1 = `bp'[1,1]
  return scalar b2 = `bp'[2,1]
  return scalar b3 = `bp'[3,1]
  return scalar v11 = `S'[1,1]
  return scalar v22 = `S'[2,2]
  return scalar v33 = `S'[3,3]
  return scalar v12 = `S'[1,2]
  return scalar v13 = `S'[1,3]
  return scalar v23 = `S'[2,3]
end

* ================= identity check on rep 1 =================
set seed 4410901
set sortseed 4410901
skpanel 0 0.5
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "ID_RC=" _rc
pstat `lam80'
matrix RP = e(roth_power_att)
dis "IDENT W     hand=" %21x r(W)     "  wooldid=" %21x e(roth_W_att)
dis "IDENT slope hand=" %21x r(slope) "  wooldid=" %21x e(roth_slope_att)
dis "IDENT se    hand=" %21x r(se)    "  wooldid=" %21x e(roth_slope_se_att)
dis "IDENT mdt80 hand=" %21x r(mdt80) "  wooldid=" %21x RP[1,2]
dis "IDENT relative errors: W " abs(r(W)/e(roth_W_att)-1) "  mdt " abs(r(mdt80)/RP[1,2]-1)
dis "IDENT dfr=" e(df_r) "  keff=" e(roth_keff_att)

* ================= timing of a second call =================
timer clear 2
timer on 2
cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
timer off 2
timer list 2

* ================= the loop =================
tempfile out
postfile PF int rep double(W pw c2 slope se mdt80 b1 b2 b3 v11 v22 v33 v12 v13 v23 pgf) using `out', replace
forvalues r = 1/`R' {
  set seed `= 900000 + 137*`r'' 
  set sortseed `= 900000 + 137*`r''
  skpanel 0 0.5
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc != 0 {
    dis "REP `r' wooldid rc=" _rc
    continue
  }
  pstat `lam80'
  local W = r(W)
  local pw = Ftail(3, 15, r(W)/3)
  local c2 = r(c2)
  local sl = r(slope)
  local se = r(se)
  local md = r(mdt80)
  local b1 = r(b1)
  local b2 = r(b2)
  local b3 = r(b3)
  local v11 = r(v11)
  local v22 = r(v22)
  local v33 = r(v33)
  local v12 = r(v12)
  local v13 = r(v13)
  local v23 = r(v23)
  * ---- the hand-rolled comparator: a plain pooled event study, cluster(st),
  * ---- 3 pre-period restrictions, Stata's own test/testparm F reference
  qui gen int rel = t - ttre
  qui gen byte nvt = (ttre == .)
  qui gen int relb = rel
  qui replace relb = -5 if rel <= -5 & !nvt
  qui replace relb =  3 if rel >=  3 & !nvt
  qui replace relb = -1 if nvt
  qui gen byte dm5 = (relb == -5 & !nvt)
  qui gen byte dm4 = (relb == -4 & !nvt)
  qui gen byte dm3 = (relb == -3 & !nvt)
  qui gen byte dm2 = (relb == -2 & !nvt)
  qui gen byte dp0 = (relb ==  0 & !nvt)
  qui gen byte dp1 = (relb ==  1 & !nvt)
  qui gen byte dp2 = (relb ==  2 & !nvt)
  qui gen byte dp3 = (relb ==  3 & !nvt)
  local pgf = .
  cap qui reghdfe y dm5 dm4 dm3 dm2 dp0 dp1 dp2 dp3, absorb(st t) cluster(st)
  if _rc == 0 {
    cap qui testparm dm4 dm3 dm2
    if _rc == 0 local pgf = r(p)
  }
  post PF (`r') (`W') (`pw') (`c2') (`sl') (`se') (`md') (`b1') (`b2') (`b3') (`v11') (`v22') (`v33') (`v12') (`v13') (`v23') (`pgf')
}
postclose PF
use `out', clear
export delimited using "s3_mc_null.csv", replace
count
dis "=== RESULTS, R = " _N " ==="
gen byte rej_w = (pw < .05)
gen byte rej_g = (pgf < .05)
qui sum rej_w
dis "wooldid joint pre-trend test size = " %8.4f r(mean) "   (n=" r(N) ")"
qui sum rej_g
dis "reghdfe testparm on 3 pre dummies size = " %8.4f r(mean) "   (n=" r(N) ")"
* GLS slope CI coverage at the reported t critical value
gen byte cov = (abs(slope) <= invttail(15,.025)*se)
qui sum cov
dis "GLS slope 95% CI coverage = " %8.4f r(mean)
qui sum slope
local sdsl = r(sd)
qui sum se
dis "mean reported slope se = " %8.4f r(mean) "   MC sd of slope = " %8.4f `sdsl' "   se/sd = " %8.4f r(mean)/`sdsl'
* empirical null distribution of W
_pctile W, p(90 95 99)
dis "empirical W quantiles: 90% " r(r1) "  95% " r(r2) "  99% " r(r3) "   (wooldid cutoff " `cutW' ")"
* mdt distribution
qui sum mdt80, detail
dis "mdt_wald(80): median " %8.5f r(p50) "  mean " %8.5f r(mean) "  sd " %8.5f r(sd) "  min " %8.5f r(min) "  max " %8.5f r(max)
dis "MDT80MEDIAN=" %12.8f r(p50)
* mdt_true from the MC covariance of the three pre coefficients
corr b1 b2 b3, cov
matrix SIGMC = r(C)
matrix MM = (-3 \ -2 \ -1)
matrix QQ = MM'*invsym(SIGMC)*MM
local c2true = QQ[1,1]
dis "c2_true (from MC covariance) = " `c2true'
local mdttrue = sqrt(`lam80'/`c2true')
dis "MDT80TRUE=" %12.8f `mdttrue'
qui sum mdt80, detail
dis "median mdt_hat / mdt_true = " %8.4f r(p50)/`mdttrue' "    (ideal-model prediction 0.9070)"
qui sum c2
dis "mean c2_hat / c2_true = " %8.4f r(mean)/`c2true' "    (ideal-model prediction 15/11 = 1.3636)"
dis "DONE"
