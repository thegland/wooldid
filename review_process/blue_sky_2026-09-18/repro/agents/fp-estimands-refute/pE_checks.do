*==============================================================================
* SKEPTIC PROBE E -- two confirmatory checks
*  E1  poisson CONTROL arm: contreatcontrols(basic) on the same poisson data,
*      to show the xb-scale truth of 0.400 IS recoverable and the untreatedZero
*      poisson arm is the deviant.
*  E2  is the r(3598) wooldid-specific?  plain reghdfe + margins, eydx on the
*      X2 data, no wooldid.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31415
set sortseed 31415
which wooldid

* same DGP as repro_uz_ringer.do
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ui = 100 + 20*runiform()
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post
gen double y  = ui + 0.5*t + 3.0*d2 + rnormal()*0.3
gen double mu = exp(1.0 + 0.10*t + 0.40*d2)
gen int ycount = rpoisson(mu)
qui sum d2 if post==1
dis "E1 POISSON TRUTH (xb scale): ame = 0.400000   att = " %10.6f 0.40*r(mean)
tempfile P
save `P'

dis _n "##### E1a poisson + contreatcontrols(basic), default engine #####"
cap noi wooldid ycount i t ttre, att contreat(dose) contreatcontrols(basic) poisson
dis "E1a_RC=" _rc
if _rc==0 dis "E1a att=" %12.7f e(b)[1,1] " ame=" %12.7f e(b)[1,2]

use `P', clear
dis _n "##### E1b poisson + contreatcontrols(basic), customvce(robust) #####"
cap noi wooldid ycount i t ttre, att contreat(dose) contreatcontrols(basic) poisson customvce(robust)
dis "E1b_RC=" _rc
if _rc==0 dis "E1b att=" %12.7f e(b)[1,1] " ame=" %12.7f e(b)[1,2]

use `P', clear
dis _n "##### E1c poisson + contreatcontrols(untreatedZero), default engine (repeat) #####"
cap noi wooldid ycount i t ttre, att contreat(d2) contreatcontrols(untreatedZero) poisson
dis "E1c_RC=" _rc
if _rc==0 dis "E1c att=" %12.7f e(b)[1,1] " ame=" %12.7f e(b)[1,2]

*==============================================================================
* E2 -- is the r(3598) margins' own behaviour or wooldid's?
*==============================================================================
clear
set obs 40
gen int i = _n
gen int ttre = .
replace ttre = 6 if i <= 16
gen double alpha = cond(i<=16, 26 + 0.5*i, 6 + 0.1*(i-16))
expand 10
bysort i: gen int t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = alpha - 20*post + rnormal()*0.3
qui sum y
dis _n "E2 min(y) = " %10.5f r(min)

qui reghdfe y post, absorb(i.i i.t)
predict double xbh, xb
qui count if xbh <= 0
dis "E2 reghdfe prediction <= 0 on " r(N) " rows"
dis "--- plain margins, eydx after reghdfe (no wooldid) ---"
cap noi margins, eydx(post) noestimcheck
dis "E2_PLAIN_RC=" _rc
if _rc==0 {
  matrix Mp = r(b)
  dis "E2 plain margins eydx = " %14.9f Mp[1,1]
}
dis "--- and with a factor-variable spelling, as wooldid uses ---"
qui reghdfe y i.post, absorb(i.i i.t)
cap noi margins, eydx(post) noestimcheck
dis "E2_PLAINFAC_RC=" _rc

dis _n "PROBE_E_DONE"
