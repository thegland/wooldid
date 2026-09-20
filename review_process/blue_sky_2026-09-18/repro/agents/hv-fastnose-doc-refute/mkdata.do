* ---------------------------------------------------------------
* Self-contained data generator for the hv-fastnose-doc refutation.
* Administrative shape: 40 states x 16 years x 20 individuals.
* dose = state-level policy generosity, CONSTANT WITHIN STATE (= within i).
* With i = state this is exactly the "dose constant within each cohort" class
* the help names at sthlp:902 (every treated cohort-period dose-slope term
* drops for multicollinearity).
* ---------------------------------------------------------------
clear all
set more off
set processors 2
set seed 20260917
set sortseed 20260917

local NS  40
local NY  16
local NI  20

set obs `=`NS'*`NY'*`NI''
gen long row = _n
gen int  st  = ceil(row/(`NY'*`NI'))
gen int  yr  = 2005 + mod(ceil(row/`NI')-1, `NY')
gen int  pid = mod(row-1, `NI') + 1

* onsets: 28 treated states in three waves, 12 never treated
gen int adopt = .
replace adopt = 2010 if inrange(st, 1, 10)
replace adopt = 2013 if inrange(st, 11, 20)
replace adopt = 2016 if inrange(st, 21, 28)

* state-level generosity, constant within state, plausible spread 0.5 - 3.0
bysort st: gen double u_st = runiform() if _n == 1
bysort st (u_st): replace u_st = u_st[1]
sort row
gen double dose = round(0.5 + 2.5*u_st, 0.01)

gen byte D = (yr >= adopt) & !missing(adopt)

bysort st: gen double afe = rnormal() if _n == 1
bysort st (afe): replace afe = afe[1]
sort row
gen double tfe = 0.03*(yr - 2005)

* structural dose effect is 0.40 -- NOT zero
gen double y = 2 + afe + tfe + 0.40*dose*D + rnormal(0, 1)

keep st yr pid adopt dose D y
label var dose "state policy generosity (constant within state)"
compress
save "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-fastnose-doc-refute/panel_states.dta", replace
di "MKDATA_OK N=" _N
