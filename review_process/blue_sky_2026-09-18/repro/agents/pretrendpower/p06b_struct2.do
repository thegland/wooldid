* ================================================================
* p06b_struct2.do -- reruns the p06 blocks that died because the
* never-treated cohort code was 0 (wooldid requires positive i),
* plus the poisson branch on a design that is actually estimable.
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

capture program drop mkstates
program define mkstates
  clear
  set obs 16
  gen long st = _n
  gen byte coh = 4
  replace coh = 1 if st <= 4
  replace coh = 2 if st > 4 & st <= 8
  replace coh = 3 if st > 8 & st <= 12
  gen int ttre = .
  replace ttre = 5 if coh == 1
  replace ttre = 7 if coh == 2
  replace ttre = 9 if coh == 3
  expand 10
  bysort st: gen int t = _n
  sort st t
  gen double f_t = 0.1*t + 0.3*sin(t)
  gen byte D = (ttre<. & t>=ttre)
  gen double e = rnormal()
  forvalues k = 2/10 {
    qui replace e = 0.6*e[_n-1] + sqrt(0.64)*rnormal() if t == `k'
  }
  gen double y = f_t + 1.0*D + e
end

* ---- B15: k = 1, small G.  The joint Wald pre-test and the
*      any-individually-significant pre-test are THE SAME TEST, so
*      "the actual size" of the latter is 0.05 by construction.
dis "TAG_B15"
mkstates
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(1) espostlength(2) pretrendpower(80)
dis "B15b_RC=" _rc
dis "B15b df_r = " e(df_r) "  npre = " e(roth_npre_att) "  nissize = " %12.8f e(roth_nissize_att)
dis "B15b 2*normal(-invttail(df_r,.025)) = " %12.8f 2*normal(-invttail(e(df_r),.025))
mat K1B = e(roth_power_att)
dis "B15b mdt_wald = " K1B[1,2] "   mdt_nis = " K1B[1,3] "   ratio = " K1B[1,3]/K1B[1,2]
dis "B15b invttail(15,.025) = " invttail(15,.025) "   invnormal(.975) = " invnormal(.975)

dis "TAG_B15C"
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "B15c_RC=" _rc
dis "B15c k=3 small-G: df_r = " e(df_r) "  nissize = " %12.8f e(roth_nissize_att)
mat K3B = e(roth_power_att)
dis "B15c mdt_wald = " K3B[1,2] "   mdt_nis = " K3B[1,3]
mat list e(V_es_att), format(%21x)

* ---- B12: which df does the Roth block use under bootstrap()+cluster()?
dis "TAG_B12"
mkstates
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80) bootstrap(standard) reps(199) rseed(9191)
dis "B12_RC=" _rc
dis "B12 bootstrap: df_r = " e(df_r) "  boot_df = " e(boot_df) "  keff = " e(roth_keff_att) "  prep = " e(roth_prep_att)
mat B12B = e(roth_power_att)
dis "B12 bootstrap mdt_wald(80) = " B12B[1,2] "  nissize = " e(roth_nissize_att)
dis "B12 boot_method = " e(boot_method) "  boot_stat = " e(boot_stat) "  roth_slope_se = " e(roth_slope_se_att)

* ---- B13: customvce(cluster) path
dis "TAG_B13"
mkstates
cap noi wooldid y coh t ttre, fe(st) customvce(cluster st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "B13_RC=" _rc
dis "B13 df_r = " e(df_r) "  prep = " e(roth_prep_att) "  nissize = " e(roth_nissize_att)

* ---- B10: poisson, on a design that is estimable (i = cohort, fe(unit))
dis "TAG_B10"
clear
set obs 120
gen long unit = _n
gen byte coh = 4
replace coh = 1 if mod(unit,4)==1
replace coh = 2 if mod(unit,4)==2
replace coh = 3 if mod(unit,4)==3
gen int ttre = .
replace ttre = 5 if coh==1
replace ttre = 7 if coh==2
replace ttre = 9 if coh==3
expand 10
bysort unit: gen int t = _n
gen byte D = (ttre<. & t>=ttre)
gen double mu = exp(2.5 + 0.03*t + 0.3*D)
gen double y = rpoisson(mu)
cap noi wooldid y coh t ttre, poisson fe(unit) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "B10_RC=" _rc
dis "B10 df_r = " e(df_r) "  keff = " e(roth_keff_att) "  prep = " e(roth_prep_att) "  nissize = " e(roth_nissize_att)
dis "TAG_BEND"
