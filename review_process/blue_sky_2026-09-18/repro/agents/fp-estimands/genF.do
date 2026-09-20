* genF.do -- Design F: a messy but plausible admin panel.
*  - t = 1..9, cohorts at 5 and 7 plus never-treated
*  - MID-WINDOW GAPS: several units miss one interior period (a reporting gap)
*  - ONSET-PERIOD GAP: one group of units misses its own relyr = 0 period
*  - one cohort has an onset AFTER the sample window (so oostreatedcontrols() bites)
clear
set obs 40
gen int i = _n
gen double ui = rnormal()*2
gen int ttre = .
replace ttre = 5  if i <= 10
replace ttre = 7  if i > 10 & i <= 18
replace ttre = 11 if i > 18 & i <= 22      // onset past the window (t max = 9)
* i in 23..40 never treated
expand 9
bysort i: gen int t = _n
* mid-window reporting gaps
drop if inlist(i,2,7,13,25,31) & t == 4
drop if inlist(i,4,15,28) & t == 6
* onset-period gap: units 6..10 do not report their own relyr = 0 year
drop if i >= 6 & i <= 10 & t == 5
gen double relyr = t - ttre
gen double tau = 0
replace tau = 1.2 + 0.3*relyr if relyr >= 0 & relyr < .
gen double y = ui + 0.2*t + tau + rnormal()*0.5
gen double w = 0.5 + runiform()*2
keep i t ttre y w
order i t ttre y w
sort i t
