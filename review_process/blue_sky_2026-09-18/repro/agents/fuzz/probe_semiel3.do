* Decisive test: do the semi-elasticity forms work when the outcome (and every fitted
* value) is comfortably STRICTLY POSITIVE?  If they still fail, the feature is broken,
* not merely undefined on my data.
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"

set seed 20260917
set sortseed 20260917
clear
set obs 60
gen long unit   = _n
gen long cohort = 1 + mod(_n, 6)
gen byte region = 1 + mod(_n, 5)
gen double a_i  = 0.3*rnormal()
gen ttre = .
replace ttre = 2003 if inlist(cohort,1,2)
replace ttre = 2005 if inlist(cohort,3,4)
gen double dose_i = runiform()*3
replace  dose_i = 0 if runiform() < .4
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t  = 0.10*(t-2000)
gen byte post   = !missing(ttre) & t >= ttre
gen double dose = cond(post, dose_i, 0)
gen double w    = 0.5 + runiform()
* STRICTLY POSITIVE outcome, far from zero: min is about 40
gen double ypos = 50 + a_i + g_t + 0.8*dose + 0.2*rnormal()
summarize ypos, detail
di "MIN_Y = " r(min)

foreach o in "contreatelasticitytype(dydx)" "contreatelasticitytype(dyex)" ///
             "contreatelasticitytype(eydx)" "contreatelasticitytype(eyex)" ///
             "semielasticity" "semielasticity nofastcollapse" ///
             "contreatelasticitytype(eydx) nofastcollapse" {
  di _n "@@SBEGIN|`o'"
  capture noisily wooldid ypos cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) `o'
  di "@@SEND|`o'|rc=" _rc
}
* and with a positive outcome under untreatedZero (the 3598 signature)
di _n "@@SBEGIN|untreatedZero+semielasticity"
capture noisily wooldid ypos cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(untreatedZero) semielasticity
di "@@SEND|untreatedZero+semielasticity|rc=" _rc
* control: does plain margins-style semielasticity work on a NON-wooldid model here?
di _n "@@SBEGIN|manual margins eydx control"
qui reg ypos c.dose i.unit i.t
capture noisily margins, eydx(dose)
di "@@SEND|manual|rc=" _rc
di "SEMIEL3_DONE"
