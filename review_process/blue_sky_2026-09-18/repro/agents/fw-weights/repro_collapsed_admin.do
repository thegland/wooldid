* SELF-CONTAINED REPRO — fw-weights/FW-02
* A researcher holds individual-level administrative microdata and collapses it to
* (unit, year) cells with a headcount, the standard way to make a large panel tractable.
* wooldid REFUSES [fw=n] and tells them to use aweights or pweights instead.  Neither
* reproduces the microdata inference: the point estimate is right, the standard error is
* not, and nothing in the output says so.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 8080
set sortseed 8080

* ---- individual-level microdata: 30 states x 12 years x ~25 people --------------
clear
set obs 30
gen state = _n
gen ttre = .
replace ttre = 2005 if state <=  8
replace ttre = 2008 if state >  8 & state <= 16
replace ttre = 2010 if state > 16 & state <= 22
gen double sfe = rnormal()
expand 12
bysort state: gen year = 2000 + _n
gen byte post = (year >= ttre) & !missing(ttre)
gen int npeople = 15 + int(runiform()*20)
expand npeople
gen double earn = 10 + sfe + 0.05*(year-2000) + 0.40*post + rnormal()
save "micro.dta", replace
count
dis "MICRO_N=" r(N)

* ---- the collapsed panel the researcher actually analyses ----------------------
use "micro.dta", clear
collapse (mean) earn (count) n = earn, by(state year ttre)
save "cells.dta", replace
count
dis "CELL_N=" r(N)

* ---- 1. the truth: run on the microdata --------------------------------------
use "micro.dta", clear
wooldid earn state year ttre, att cluster(state)
scalar b_micro = e(b)[1,1]
dis "MICRO_b=" b_micro "  se=" %12.8f sqrt(e(V)[1,1]) "  N=" e(N)
scalar se_micro = sqrt(e(V)[1,1])

* ---- 2. what the user asks for, and is refused -------------------------------
use "cells.dta", clear
cap noi wooldid earn state year ttre [fw=n], att cluster(state)
dis "FW_RC=" _rc

* ---- 3. the remedy the error message prescribes ------------------------------
use "cells.dta", clear
wooldid earn state year ttre [aw=n], att cluster(state)
scalar b_aw = e(b)[1,1]
dis "AW_b=" b_aw "  se=" %12.8f sqrt(e(V)[1,1]) "  N=" e(N)
scalar se_aw = sqrt(e(V)[1,1])

use "cells.dta", clear
wooldid earn state year ttre [pw=n], att cluster(state)
dis "PW_b=" e(b)[1,1] "  se=" %12.8f sqrt(e(V)[1,1]) "  N=" e(N)
scalar se_pw = sqrt(e(V)[1,1])

* ---- 3b. the same three runs with the DEFAULT iid vce ------------------------
use "micro.dta", clear
qui wooldid earn state year ttre, att
scalar se_micro_iid = sqrt(e(V)[1,1])
dis "MICRO_IID_se=" %12.8f se_micro_iid "  N=" e(N)
use "cells.dta", clear
qui wooldid earn state year ttre [aw=n], att
scalar se_aw_iid = sqrt(e(V)[1,1])
dis "AW_IID_se=" %12.8f se_aw_iid "  N=" e(N)
use "cells.dta", clear
qui wooldid earn state year ttre [pw=n], att
scalar se_pw_iid = sqrt(e(V)[1,1])
dis "PW_IID_se=" %12.8f se_pw_iid "  N=" e(N)

* ---- 4. verdict ---------------------------------------------------------------
dis "B_RELDIFF_aw_vs_micro = " %12.4e abs(b_aw-b_micro)/abs(b_micro)
dis "SE_RATIO aw/micro     = " %8.4f se_aw/se_micro
dis "SE_RATIO pw/micro     = " %8.4f se_pw/se_micro
dis "SE_RATIO aw/micro IID = " %8.4f se_aw_iid/se_micro_iid
dis "SE_RATIO pw/micro IID = " %8.4f se_pw_iid/se_micro_iid
dis "REPRO_DONE"
