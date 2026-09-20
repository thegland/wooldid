* ---------------------------------------------------------------
* mkdata.do -- self-contained generator used by every probe here.
* Writes sk_base.dta: a plain staggered panel, 40 units x 10 years,
* cohorts 2004 (13 units), 2006 (13 units), never (14 units).
* No missingness, no weirdness. This is the "well-behaved" arm.
* ---------------------------------------------------------------
clear all
set more off
set seed 424242
set sortseed 424242
set obs 40
gen long unit = _n
gen int  ttre = .
replace  ttre = 2004 if unit <= 13
replace  ttre = 2006 if unit > 13 & unit <= 26
gen double ufe = rnormal()
gen long state = ceil(unit/4)          /* 10 states, 4 units each: a COARSER cluster */
expand 10
bysort unit: gen int year = 1999 + _n
gen double tfe = 0.12*(year-1999)
gen double te  = cond(missing(ttre), 0, cond(year >= ttre, 0.5 + 0.1*(year-ttre), 0))
gen double y   = 1 + ufe + tfe + te + rnormal()*0.7
gen double x1  = rnormal()
label var y "outcome"
compress
save "sk_base.dta", replace
