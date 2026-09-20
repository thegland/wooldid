* s5_ex7_ati.do -- SELF-CONTAINED.
* PART 1 (PTP-4): wooldid.sthlp example 7 verbatim, with i = the unit id (the natural
*   reading of "wooldid y i t ttre"), on an ordinary staggered panel.
* PART 2 (PTP-5): null Monte Carlo of the joint pre-trend test on the at-i geometry
*   (i = unit, cluster(unit)), against the ideal-model prediction for the same k and df_r.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop unitpanel
prog def unitpanel
  args nunit delta rho
  clear
  set obs `nunit'
  gen long i = _n
  gen byte g = 1 + mod(_n-1, 4)
  gen int ttre = .
  replace ttre = 6  if g == 1
  replace ttre = 8  if g == 2
  replace ttre = 10 if g == 3
  gen double alpha_i = rnormal()*0.8
  expand 12
  bysort i: gen int t = _n
  sort i t
  gen double yeareff = 0.12*t - 0.02*t*t + 0.25*cos(t)
  gen byte D = (ttre < . & t >= ttre)
  gen double relyr = t - ttre
  gen double trendterm = cond(ttre < ., relyr, 0)
  gen double e = rnormal()
  forvalues k = 2/12 {
    qui replace e = `rho'*e[_n-1] + sqrt(1-`rho'*`rho')*rnormal() if t == `k'
  }
  gen double y = 2 + alpha_i + yeareff + 0.8*D + `delta'*trendterm + e
end

* ================= PART 1: help example 7, i = unit =================
dis "=========== PART1 help example 7 (i = unit) ==========="
set seed 331177
set sortseed 331177
unitpanel 60 0 0.5
cap noi wooldid y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)
dis "P1_RC=" _rc
dis "P1 roth_ok  = " e(roth_ok)
dis "P1 roth_W   = " e(roth_W_att)
dis "P1 roth_prep= " e(roth_prep_att)
dis "P1 df_r     = " e(df_r) "   N_clust=" e(N_clust)
cap mat list e(roth_power_att), format(%12.6g)

* same run with a FINER cluster than i -- the caution's own recommended remedy
dis "=========== PART1b same design, cluster strictly finer than the cell ==========="
cap noi wooldid y i t ttre, att robust espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)
dis "P1b_RC=" _rc

* ================= PART 2: null MC on the at-i geometry =================
dis "=========== PART2 at-i null MC ==========="
local R 150
local lam80 = npnF(3, 15, invF(3, 15, .95), .20)
tempfile out
postfile PF int rep double(W pw dfr slope se) using `out', replace
forvalues r = 1/`R' {
  set seed `= 550000 + 211*`r''
  set sortseed `= 550000 + 211*`r''
  unitpanel 60 0 0.5
  cap qui wooldid y i t ttre, cluster(i) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc != 0 {
    dis "REP `r' rc=" _rc
    continue
  }
  local dfr = e(df_r)
  tempname B V S bp Si m w
  matrix `B'  = e(b_es_att)
  matrix `V'  = e(V_es_att)
  matrix `S'  = `V'[1..3, 1..3]
  matrix `S'  = (`S' + `S'')/2
  matrix `bp' = (`B'[1,1] \ `B'[1,2] \ `B'[1,3])
  matrix `Si' = invsym(`S')
  matrix `m'  = (-3 \ -2 \ -1)
  matrix `w'  = `bp''*`Si'*`bp'
  local W = `w'[1,1]
  matrix `w'  = `m''*`Si'*`m'
  local c2 = `w'[1,1]
  matrix `w'  = `m''*`Si'*`bp'
  local sl = `w'[1,1]/`c2'
  local se = 1/sqrt(`c2')
  local pw = Ftail(3, `dfr', `W'/3)
  post PF (`r') (`W') (`pw') (`dfr') (`sl') (`se')
}
postclose PF
use `out', clear
export delimited using "s5_ati_null.csv", replace
count
qui sum dfr
local dfr = r(mean)
dis "df_r (mean) = " `dfr'
gen byte rej = (pw < .05)
qui sum rej
dis "AT-I joint pre-trend test size = " %8.4f r(mean) "  (n=" r(N) ")"
gen byte cov = (abs(slope) <= invttail(`dfr',.025)*se)
qui sum cov
dis "AT-I GLS slope 95% CI coverage = " %8.4f r(mean)
qui sum slope
local sdsl = r(sd)
qui sum se
dis "AT-I mean reported se = " %8.4f r(mean) "  MC sd = " %8.4f `sdsl' "  se/sd = " %8.4f r(mean)/`sdsl'
dis "IDEAL-MODEL prediction for k=3, df_r=`dfr': size = " %8.4f Ftail(3, `dfr'-2, 3*invF(3,`dfr',.95)*(`dfr'-2)/(`dfr'*3))
dis "DONE"
