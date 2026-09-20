* p13 — calibration: does Stata's OWN bootstrap prefix show the same k-dim Wald inflation?
clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916
clear
set obs 40
gen int state = _n
gen double sfe = rnormal()
expand 20
bysort state: gen int t = _n
gen double y = 1 + sfe + 0.05*t + rnormal()*0.8
* no xtset: the resampled frame has duplicate panel ids
dis "##### analytic cluster-robust #####"
qui regress y i.t, cluster(state)
dis "AN_DFR=" e(df_r)
qui test (2.t=0) (3.t=0) (4.t=0) (5.t=0) (6.t=0) (7.t=0)
dis "AN_TEST6 p=" r(p) " df=" r(df) " drop=" r(drop)
foreach R in 8 15 50 400 {
  dis "##### stata bootstrap prefix reps(`R') #####"
  cap noi bootstrap, reps(`R') cluster(state) seed(4242) nodots nowarn: regress y i.t
  dis "SBP`R'_RC=" _rc
  if _rc == 0 {
    matrix SV = e(V)
    mata: V=st_matrix("SV"); V=(V+V')/2; st_numscalar("rk",rank(V))
    dis "SBP`R'_RANK=" rk " of " colsof(SV) " NREPS=" e(N_reps) " DFR=" e(df_r)
    cap noi test (2.t=0) (3.t=0) (4.t=0) (5.t=0) (6.t=0) (7.t=0)
    dis "SBP`R'_TEST6 p=" r(p) " df=" r(df) " drop=" r(drop) " chi2/F=" r(chi2) r(F)
  }
}
dis "P13_DONE"
