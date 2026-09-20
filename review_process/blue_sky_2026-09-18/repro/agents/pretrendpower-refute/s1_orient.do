* s1_orient.do -- SKEPTIC orientation. Own data generator, own geometry.
* Answers cheaply: PTP-2 (nissize formula), PTP-3 (poisson chi2 branch),
* PTP-4 (help example 7 caution), and times one wooldid run for the MC budget.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103
which wooldid

* ---------- my own 16-state policy panel: 4 onset cohorts of 4 states,
*            T=12, state random effects, common year path, AR(1) within state.
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

* =========== A. gaussian run, k=3 pre periods ===========
dis "=========== A: gaussian k=3 ==========="
skpanel 0 0.5
timer clear 1
timer on 1
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(50 80)
dis "A_RC=" _rc
timer off 1
timer list 1
dis "A df_r      = " e(df_r)
dis "A N_clust   = " e(N_clust)
dis "A roth_W    = " %21x e(roth_W_att)
dis "A roth_prep = " e(roth_prep_att)
dis "A keff      = " e(roth_keff_att)
dis "A nissize   = " %12.8f e(roth_nissize_att)
mat list e(roth_power_att), format(%12.6g)
* what the reported nissize would be under the 3-dim normal vs the t
dis "A 2*normal(-invttail(df,.025)) at k=1 equivalent = " %12.8f 2*normal(-invttail(e(df_r),.025))

* =========== B. gaussian run, k=1 pre period (PTP-2's exact case) ===========
dis "=========== B: gaussian k=1 ==========="
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(1) espostlength(2) pretrendpower(80)
dis "B_RC=" _rc
local bdf = e(df_r)
dis "B df_r    = " `bdf'
dis "B nissize = " %20.12f e(roth_nissize_att)
dis "B 2*normal(-invttail)= " %20.12f 2*normal(-invttail(`bdf',.025))
dis "B 2*ttail(df,invttail)= " %20.12f 2*ttail(`bdf',invttail(`bdf',.025))
dis "B DIFF = " %20.12f e(roth_nissize_att) - 2*normal(-invttail(`bdf',.025))
mat list e(roth_power_att), format(%12.6g)

* =========== C. poisson on the same design ===========
dis "=========== C: poisson ==========="
skpanel 0 0.5
gen double mu = exp(1 + 0.3*alpha_s + 0.05*t + 0.4*D)
gen long ycnt = rpoisson(mu)
cap noi wooldid ycnt coh t ttre, poisson fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "C_RC=" _rc
dis "C df_r    = " e(df_r)
dis "C N_clust = " e(N_clust)
dis "C e(V) dim/df check: e(df_m)=" e(df_m) "  e(N)=" e(N)
dis "C roth_W  = " e(roth_W_att)
dis "C prep    = " e(roth_prep_att)
mat list e(roth_power_att), format(%12.6g)
* what the F branch would have produced with df_r = N_clust-1
local kk = e(roth_keff_att)
local nc = e(N_clust)
local W  = e(roth_W_att)
local c2 = e(roth_slope_se_att)
dis "C keff=" `kk' "  slope_se=" `c2'
local lamF = npnF(`kk', `nc'-1, invF(`kk',`nc'-1,.95), 1-80/100)
local lamX = npnchi2(`kk', invchi2(`kk',.95), 1-80/100)
dis "C lambda chi2 = " `lamX' "   lambda F(df=" `nc'-1 ") = " `lamF'
dis "C mdt chi2 = " sqrt(`lamX')*`c2' "   mdt F = " sqrt(`lamF')*`c2'
dis "C ratio F/chi2 = " sqrt(`lamF'/`lamX')
* also: does gaussian on the SAME count outcome post df_r?
cap noi wooldid ycnt coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "C2_RC=" _rc
dis "C2 df_r = " e(df_r)

* =========== D. does poisson WITHOUT cluster post df_r? ===========
dis "=========== D: poisson robust ==========="
cap noi wooldid ycnt coh t ttre, poisson fe(st) robust esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "D_RC=" _rc
dis "D df_r = " e(df_r)
* and ppmlhdfe itself:
cap noi ppmlhdfe ycnt i.t D, absorb(st) cluster(st)
dis "D ppmlhdfe df_r = " e(df_r) "  N_clust=" e(N_clust)
cap noi poisson ycnt i.t D i.st, cluster(st)
dis "D poisson df_r = " e(df_r)

* =========== E. help example 7 verbatim, variables named as the help writes them ===========
dis "=========== E: help example 7 ==========="
skpanel 0 0.5
rename coh i
rename ttre ttre0
gen int ttre = ttre0
cap noi wooldid y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)
dis "E_RC=" _rc
dis "DONE"
