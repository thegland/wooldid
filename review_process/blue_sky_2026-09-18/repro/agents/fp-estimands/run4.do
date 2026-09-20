* run4.do -- pin the mechanism behind the semielasticity engine split.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260919
set sortseed 20260919
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"

use "dataE.dta", clear
* rebuild wooldid's internal relative-time variables exactly (esprelength(3))
gen double relyr = t - ttre
gen double _t = relyr
replace _t = . if _t < -3
replace _t = . if t == 1
qui sum _t, meanonly
local mint = max(0, -r(min))
dis "MINT=" `mint'
replace _t = _t + `mint' + 1
replace _t = 0 if _t == .
gen double _i = i if _t > 0 & _t < .
replace _i = 0 if _i == .
gen byte tf = (_t > 0 & _t < .)

* the stage-1 model, both engines
qui reghdfe y i._i#i._t#i.tf, absorb(i.i i.t) tol(1e-12)
predict double xb_hdfe, xb
qui reg y i.i i.t i._i#i._t#i.tf
predict double xb_full, xb

* per-row treatment effect = xb(tf=1) - xb(tf=0); with the triple interaction this is
* just the cell coefficient, which we read off as the deviation from the base cell.
* Build it directly from the full-dummy fit by zeroing the cell.
gen double beta_row = .
qui levelsof _t if tf == 1, local(hl)
foreach h of local hl {
  qui levelsof _i if tf == 1 & _t == `h', local(gl)
  foreach g of local gl {
    qui replace beta_row = _b[`g'._i#`h'._t#1.tf] if _i == `g' & _t == `h' & tf == 1
  }
}
qui replace beta_row = 0 if beta_row == .

gen byte post = (t >= ttre) & ttre < . & tf == 1
qui count if post
dis "NPOST=" r(N)

* candidate aggregations of the semi-elasticity over the att subpop
foreach v in xb_hdfe xb_full {
  qui gen double se_`v' = beta_row/`v'
  qui sum se_`v' if post, meanonly
  dis "CAND_dydx_over_`v'=" %21x r(mean)
  qui gen double lr_`v' = ln(`v') - ln(`v' - beta_row)
  qui sum lr_`v' if post, meanonly
  dis "CAND_lndiff_`v'=" %21x r(mean)
  qui sum `v' if post, meanonly
  dis "MEAN_`v'=" %21x r(mean)
}
qui sum beta_row if post, meanonly
dis "MEAN_beta=" %21x r(mean)
dis "RATIO_meanbeta_over_mean_xbhdfe=" %21x 0
qui sum y if post, meanonly
dis "MEAN_y_post=" %21x r(mean)

* what wooldid reported (E1 reghdfe, E2 explicit dummies)
dis "WD_E1_main_att=3.666814271565411e-02"
dis "WD_E2_main_att=5.186377716829238e-02"

* margins directly, both engines, to see what margins itself returns
qui reghdfe y i._i#i._t#i.tf, absorb(i.i i.t) tol(1e-12)
cap noi margins , eydx(tf) subpop(if tf==1 & t >= ttre) noestimcheck nose
dis "MARG_HDFE_RC=" _rc
if _rc == 0 {
  matrix Mh = r(b)
  dis "MARG_HDFE_b=" %21x Mh[1, colsof(Mh)]
}
qui reg y i.i i.t i._i#i._t#i.tf
cap noi margins , eydx(tf) subpop(if tf==1 & t >= ttre) noestimcheck nose
dis "MARG_FULL_RC=" _rc
if _rc == 0 {
  matrix Mf = r(b)
  dis "MARG_FULL_b=" %21x Mf[1, colsof(Mf)]
}
dis "ALLDONE_RUN4"
