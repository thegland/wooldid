*==============================================================================
* SELF-CONTAINED REPRO -- FP-SEMIELAS-FE
*
* CLAIM: under the DEFAULT estimation engine (reghdfe, with i and t absorbed),
* `semielasticity` divides the treatment effect by a prediction that EXCLUDES the
* absorbed unit and time fixed effects.  The reported semi-elasticity is therefore
* the effect measured against the sample's grand mean, not against each
* observation's own predicted outcome, which is what the help promises: "the
* average percentage change ... between the two sets of predicted outcomes"
* (sthlp:243-245).  No caution is printed.
*
* The identical model estimated with EXPLICIT i.i / i.t dummies -- which is what
* customvce() and unconditionalse do -- returns the correct number.
*
* Run:  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_semielas_fe.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , spread(real)
  clear
  set obs 30
  gen int i = _n
  gen double ui = 120 - `spread'/2 + `spread'*runiform()
  gen int ttre = .
  replace ttre = 4 if i <= 10
  replace ttre = 6 if i > 10 & i <= 20
  expand 8
  bysort i: gen int t = _n
  gen double tau = 5*((t >= ttre) & ttre < .)   // constant level effect of 5
  gen double y = ui + tau + rnormal()*1.5       // NO calendar trend: isolates the unit FE
end

cap prog drop armreport
prog def armreport
  syntax , label(string)
  cap noi wooldid y i t ttre, att
  local lev_h = e(b)[1,1]
  cap noi wooldid y i t ttre, att customvce(robust)
  local lev_f = e(b)[1,1]
  cap noi wooldid y i t ttre, att semielasticity
  local se_h = e(b)[1,1]
  cap noi wooldid y i t ttre, att semielasticity customvce(robust)
  local se_f = e(b)[1,1]
  cap noi wooldid y i t ttre, att semielasticity unconditionalse robust
  local se_u = e(b)[1,1]
  dis "`label' SEMIEL  unconditionalse robust (third explicit-dummy arm) = " %12.9f `se_u'
  dis _n "`label' LEVEL   reghdfe=" %12.9f `lev_h' "  explicit=" %12.9f `lev_f' ///
         "  ratio=" %9.6f `lev_f'/`lev_h'
  dis    "`label' SEMIEL  reghdfe=" %12.9f `se_h'  "  explicit=" %12.9f `se_f'  ///
         "  ratio=" %9.6f `se_f'/`se_h'
end

dis _n "=========== ARM 1: unit levels spread over 120 +/- 110 ==========="
mkpanel , spread(220)
tempfile a1
save `a1'
armreport , label(ARM1)

dis _n "=========== ARM 2 (CONTROL): every unit at the same level ==========="
mkpanel , spread(0)
armreport , label(ARM2)

*----------------------------------------------------------------- mechanism
* Rebuild wooldid's internal relative-time variables for the DEFAULT window
* (esprelength(0)) and show by hand which prediction each engine divides by.
use `a1', clear
gen double _t = t - ttre
qui bysort i (t): gen byte first = (_n == 1)
replace _t = . if _t < 0
replace _t = . if first == 1
replace _t = _t + 1
replace _t = 0 if _t == .
gen double _i = i if _t > 0 & _t < .
replace _i = 0 if _i == .
gen byte tf = (_t > 0 & _t < .)

qui reghdfe y i._i#i._t#i.tf, absorb(i.i i.t) tol(1e-12)
predict double xb_hdfe, xb
qui reg y i.i i.t i._i#i._t#i.tf
predict double xb_full, xb
gen double beta_row = 0
qui levelsof _t if tf == 1, local(hl)
foreach h of local hl {
  qui levelsof _i if tf == 1 & _t == `h', local(gl)
  foreach g of local gl {
    qui replace beta_row = _b[`g'._i#`h'._t#1.tf] if _i == `g' & _t == `h' & tf == 1
  }
}
gen byte post = (t >= ttre) & ttre < . & tf == 1
qui corr xb_hdfe xb_full
dis _n "CORR(reghdfe xb , true fitted value) = " %9.6f r(rho)
qui sum xb_hdfe if post, meanonly
local m1 = r(sd)
qui sum xb_hdfe if post
dis "SD(reghdfe xb | post)      = " %10.5f r(sd)
qui sum xb_full if post
dis "SD(true fitted value|post) = " %10.5f r(sd)
qui gen double lh = ln(xb_hdfe) - ln(xb_hdfe - beta_row)
qui gen double lf = ln(xb_full) - ln(xb_full - beta_row)
qui sum lh if post, meanonly
dis "HAND mean[ln(xb_hdfe)-ln(xb_hdfe-b)] = " %12.9f r(mean) "  <- reproduces the DEFAULT answer"
qui sum lf if post, meanonly
dis "HAND mean[ln(yhat)   -ln(yhat   -b)] = " %12.9f r(mean) "  <- reproduces the EXPLICIT answer"

dis _n "EXPECTED: ARM2 SEMIEL ratio ~= 1 (nothing is absorbed), ARM1 ratio far from 1;"
dis "LEVEL ratio must be 1 in both arms (dydx cancels the fixed effects)."
dis "REPRO_DONE"
