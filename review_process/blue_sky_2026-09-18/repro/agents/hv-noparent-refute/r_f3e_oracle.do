* ============================================================================
* SKEPTIC FOLLOW-UP 4: the margins-eydx oracle DONE RIGHT.
* r_f3b's oracle failed to straddle zero (min yhat = .56) because u was
* per-observation noise that reg averaged out. Here u is GROUP-level, so the
* fitted values for x==1 genuinely span negative and positive.
* Question: on a straddling domain, is margins' eydx the mean of per-observation
* ratios dydx/yhat (a complete but meaningless number), the ratio of means, or
* does it refuse?  No wooldid involved.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f3e_oracle.do
* ============================================================================
clear all
set more off
set processors 2
set seed 8123
set sortseed 8123

clear
set obs 400
gen int g = mod(_n,20)
gen byte x = (_n>200)
gen double u = -3 + 6*(g/19)
gen double y = u + 1.0*x + rnormal()*0.05
qui reg y x i.g
predict double yh, xb
qui su yh if x==1
dis "STRADDLE_MINYHAT=" r(min)
dis "STRADDLE_MAXYHAT=" r(max)
scalar meanyh = r(mean)
count if x==1 & yh < 0
dis "STRADDLE_NNEG=" r(N)
gen double rat = _b[x]/yh
qui su rat if x==1
dis "STRADDLE_MEANOFRATIOS=" r(mean)
dis "STRADDLE_RATIOOFMEANS=" _b[x]/meanyh
capture noisily margins, eydx(x) subpop(if x==1) noestimcheck
dis "STRADDLE_MARGINS_RC=" _rc
capture noisily matrix list r(table)

* tighten the straddle so some yhat sit very close to zero
clear
set obs 400
gen int g = mod(_n,20)
gen byte x = (_n>200)
gen double u = -0.5 + 1.0*(g/19)
gen double y = u + 0.2*x + rnormal()*0.02
qui reg y x i.g
predict double yh, xb
qui su yh if x==1
dis "TIGHT_MINYHAT=" r(min)
dis "TIGHT_MAXYHAT=" r(max)
scalar meanyh2 = r(mean)
gen double rat = _b[x]/yh
qui su rat if x==1
dis "TIGHT_MEANOFRATIOS=" r(mean)
dis "TIGHT_RATIOOFMEANS=" _b[x]/meanyh2
capture noisily margins, eydx(x) subpop(if x==1) noestimcheck
dis "TIGHT_MARGINS_RC=" _rc
capture noisily matrix list r(table)

dis "DONE_ORACLE"
