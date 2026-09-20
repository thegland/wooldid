*==============================================================================
* SKEPTIC PROBE C -- three things in one run
*  C1  FP-CTZERO-ENGINE mechanism: verbose on the EXPLICIT arm (which term(s) get
*      dropped) + the unconditionalse route, on CLEAN data.
*  C2  FP-SEMIELAS-FE escalation: a realistic design whose reghdfe prediction
*      _cons+Xb goes NEGATIVE while every true fitted value is positive.
*  C3  FP-NOPRE-LATEENTRY: a late-entering unit whose onset sits inside the
*      sample window but at its own first observed period.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do pC_mixed.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
which wooldid

*==============================================================================
* C1
*==============================================================================
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ui_unit = 100 + 20*runiform()
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post
gen double y = ui_unit + 0.5*t + 3.0*d2 + rnormal()*0.3
qui sum d2 if post==1
dis "C1 TRUTH: E[dose|post]=" %9.6f r(mean) "  TRUE att=" %9.6f 3*r(mean) "  TRUE ame=3.000000"
tempfile c1
save `c1'

dis _n "===== C1a EXPLICIT customvce(robust) VERBOSE ====="
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) customvce(robust) verbose
dis "C1a_RC=" _rc
if _rc==0 dis "C1a att=" %13.7f e(b)[1,1] " ame=" %13.7f e(b)[1,2]

use `c1', clear
dis _n "===== C1b unconditionalse robust ====="
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) unconditionalse robust
dis "C1b_RC=" _rc
if _rc==0 dis "C1b att=" %13.7f e(b)[1,1] " ame=" %13.7f e(b)[1,2]

use `c1', clear
dis _n "===== C1c DEFAULT reghdfe VERBOSE ====="
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) verbose
dis "C1c_RC=" _rc
if _rc==0 dis "C1c att=" %13.7f e(b)[1,1] " ame=" %13.7f e(b)[1,2]

*==============================================================================
* C2  semielasticity where reghdfe's _cons+Xb goes negative
*     every y > 0, every TRUE fitted value > 0, big negative treatment effect,
*     treated units high-baseline and controls low-baseline.
*==============================================================================
clear
set obs 40
gen int i = _n
gen int ttre = .
replace ttre = 6 if i <= 16
gen double alpha = cond(i<=16, 20 + 0.5*i, 2 + 0.1*(i-16))
expand 10
bysort i: gen int t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = alpha - 20*post + rnormal()*0.3
qui sum y
dis "C2: min(y)=" %9.5f r(min) "  mean(y)=" %9.5f r(mean)
qui count if y <= 0
dis "C2: n(y<=0) = " r(N)
tempfile c2
save `c2'

dis _n "===== C2a SEMIELAS default reghdfe ====="
cap noi wooldid y i t ttre, att semielasticity
dis "C2a_RC=" _rc
if _rc==0 dis "C2a semiel=" %14.9f e(b)[1,1]

use `c2', clear
dis _n "===== C2b SEMIELAS customvce(robust) ====="
cap noi wooldid y i t ttre, att semielasticity customvce(robust)
dis "C2b_RC=" _rc
if _rc==0 dis "C2b semiel=" %14.9f e(b)[1,1]

use `c2', clear
dis _n "===== C2c mechanism: reghdfe xb vs full fitted ====="
qui reghdfe y post, absorb(i.i i.t)
predict double xbh, xb
qui sum xbh
dis "C2c min(xb_hdfe)=" %10.5f r(min) "  max=" %10.5f r(max)
qui count if xbh <= 0
dis "C2c n(xb_hdfe <= 0) = " r(N)
qui reg y i.i i.t post
predict double xbf, xb
qui sum xbf
dis "C2c min(xb_full)=" %10.5f r(min) "  max=" %10.5f r(max)
qui count if xbf <= 0
dis "C2c n(xb_full <= 0) = " r(N)

*==============================================================================
* C3  late entry: unit first observed at t = its own onset
*==============================================================================
clear
set obs 30
gen int i = _n
gen int ttre = .
replace ttre = 5 if i <= 10
replace ttre = 4 if i > 10 & i <= 16       // these will ENTER LATE at t = 4
gen double alpha = 50 + i
expand 8
bysort i: gen int t = _n
* six units enter the panel late, at t = 4, which is exactly their onset
drop if i > 10 & i <= 16 & t < 4
gen byte post = (t >= ttre) & ttre < .
gen double eff = cond(i > 10 & i <= 16, 6, 2)
gen double y = alpha + 0.4*t + eff*post + rnormal()*0.5
qui count
dis "C3: total rows = " r(N)
tab ttre, m
tempfile c3
save `c3'

dis _n "===== C3a default run ====="
cap noi wooldid y i t ttre, att
dis "C3a_RC=" _rc
if _rc==0 dis "C3a att=" %13.7f e(b)[1,1] "  N=" e(N) "  Ndrop=" e(N_dropped)

use `c3', clear
dis _n "===== C3b same, but late entrants observed from t=3 (one pre period) ====="
* rebuild with the late entrants observed from t=3 so they survive
clear
set obs 30
gen int i = _n
gen int ttre = .
replace ttre = 5 if i <= 10
replace ttre = 4 if i > 10 & i <= 16
gen double alpha = 50 + i
expand 8
bysort i: gen int t = _n
drop if i > 10 & i <= 16 & t < 3
gen byte post = (t >= ttre) & ttre < .
gen double eff = cond(i > 10 & i <= 16, 6, 2)
gen double y = alpha + 0.4*t + eff*post + rnormal()*0.5
cap noi wooldid y i t ttre, att
dis "C3b_RC=" _rc
if _rc==0 dis "C3b att=" %13.7f e(b)[1,1] "  N=" e(N) "  Ndrop=" e(N_dropped)

dis _n "PROBE_C_DONE"
