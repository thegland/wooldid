* genC.do -- Design C: a REALISTIC admin panel (t = 1..8) in which one cohort enters the
* panel in the same period as its own treatment onset (a late-entering unit, or a unit whose
* record begins in the adoption year).  Cohort 3's first observed period IS its onset, so the
* pooled reference regime forces that row to be the reference (wooldid_reltime ado:12087).
* Effects differ by cohort so any mis-weighting across cohorts moves the printed number.
clear
set obs 36
gen int i = _n
gen double ui = rnormal()*2
gen int ttre = .
replace ttre = 4 if i <= 9              // ordinary cohort, pre-periods t=1..3
replace ttre = 6 if i > 9  & i <= 18    // ordinary cohort, pre-periods t=1..5
replace ttre = 3 if i > 18 & i <= 24    // enters the panel AT its own onset (see drop below)
* i in 25..36 never treated
expand 8
bysort i: gen int t = _n
* cohort 3 (i 19..24) is only observed from t = 3, which is its onset period
drop if i > 18 & i <= 24 & t < 3
gen double relyr = t - ttre
gen double tau = 0
replace tau = 1 if relyr >= 0 & relyr < . & i <= 18
replace tau = 3 if relyr >= 0 & relyr < . & i > 18 & i <= 24
gen double y = ui + 0.25*t + tau + rnormal()*0.4
gen double w = 0.5 + runiform()*2
gen double cwm = 0.4 + runiform()*2
keep i t ttre y w cwm
order i t ttre y w cwm
sort i t
