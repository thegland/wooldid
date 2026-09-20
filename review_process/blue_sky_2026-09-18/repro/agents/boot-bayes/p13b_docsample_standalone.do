*=============================================================================
* p13b  SELF-CONTAINED REPRO (needs only the pinned build on the adopath).
*
* CLAIM: sthlp:557 says that under bootalpha(size) the bootsave() file's _rep = 0
* row "holds 1/N, N the number of observations in the estimation sample". It does
* not: ado:437 divides by the MAP-FRAME row count (the user-if sample intersected
* with nonmissing cluster()), which still contains every row the estimator drops.
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
which wooldid
set seed 20260916
set sortseed 20260916
clear
set obs 30
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2006 if st <= 7
replace onset = 2009 if st > 7 & st <= 14
replace onset = 2012 if st > 14 & st <= 20
expand 16
bysort st: gen int year = 1999 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.04*(year-2000) + 1.3*post + rnormal()*0.6
* ordinary administrative-panel item nonresponse on the outcome
gen double u = runiform()
replace y = . if u < 0.15
drop u
qui count
local NROWS = r(N)
qui count if !missing(y)
dis "rows in the data       = `NROWS'"
dis "rows with a non-missing outcome = " r(N)
qui wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(20) rseed(1234) bootsave(p13b_size) replace dots(0)
local eN = e(N)
dis "e(N) (estimation sample) = `eN'"
use p13b_size, clear
qui ds w_*
local w1 : word 1 of `r(varlist)'
qui su `w1' if _rep == 0
dis "bootsave _rep=0 value  = " %21x r(mean) "  = " r(mean)
dis "1 / e(N)               = " %21x 1/`eN'    "  = " 1/`eN'
dis "1 / rows in the data   = " %21x 1/`NROWS' "  = " 1/`NROWS'
dis "MATCHES e(N)?            " (reldif(r(mean), 1/`eN') < 1e-12)
dis "MATCHES data row count?  " (reldif(r(mean), 1/`NROWS') < 1e-12)
dis "PROBE_DONE"
