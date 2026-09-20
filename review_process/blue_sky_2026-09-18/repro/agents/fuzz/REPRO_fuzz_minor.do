*=======================================================================
* SELF-CONTAINED REPRO -- fuzz agent, the two low-severity items
*  (a) dots(-5) is refused under bootstrap() and silently accepted under permute()
*  (b) `update' on a full command line runs the installer and posts no e()
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_fuzz_minor.do
*=======================================================================
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 60
gen long unit   = _n
gen long cohort = 1 + mod(_n, 6)
gen byte region = 1 + mod(_n, 5)
gen double a_i  = rnormal()
gen ttre = .
replace ttre = 2003 if inlist(cohort,1,2)
replace ttre = 2005 if inlist(cohort,3,4)
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t = 0.3*(t-2000)
gen byte post  = !missing(ttre) & t >= ttre
gen double y   = a_i + g_t + 1.2*post + 0.3*rnormal()

di _n "=== (a) dots(-5) ==="
capture noisily qui wooldid y unit t ttre, att cluster(cohort) permute(shuffle) permreps(25) dots(-5)
di "   permute + dots(-5)   -> rc " _rc
capture noisily qui wooldid y unit t ttre, att cluster(cohort) bootstrap(standard) reps(25) rseed(3) dots(-5)
di "   bootstrap + dots(-5) -> rc " _rc
di "   EXPECTED: the same refusal on both."

di _n "=== (b) update on a full command line ==="
ereturn clear
capture noisily wooldid y unit t ttre, att robust update
di "   rc = " _rc
di "   e(cmd) = '" "`e(cmd)'" "'   (empty means nothing was estimated)"
di "   EXPECTED: either estimate and post e(), or refuse the combination."
