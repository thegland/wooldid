*=======================================================================
* SELF-CONTAINED REPRO -- fuzz agent
* CLAIM: wooldid refuses level(99.9) and level(99.99) with
*        "level() can have at most two digits after the decimal point",
*        although (a) the input has one and two decimals respectively,
*        (b) plain Stata's `set level 99.9' / `set level 99.99' both
*        succeed, and (c) wooldid's OWN validation message says
*        "level() must be between 10 and 99.99", i.e. it advertises
*        99.99 as legal.  level(97.5), level(90.25) and level(99.5)
*        -- the levels that are exactly representable in binary --
*        all work, which is the signature of the level being re-rendered
*        with full double precision before `set level'.
* RUN:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_fuzz_level999.do
*=======================================================================
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

di _n "=== plain Stata accepts every one of these ==="
foreach L in 99.9 99.99 97.5 90.25 99.5 {
  capture noisily set level `L'
  di "   set level `L' -> rc " _rc "   c(level) = " c(level)
}
set level 95

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
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t = 0.10*(t-2000)
gen byte post  = !missing(ttre) & t >= ttre
gen double y   = 50 + a_i + g_t + 1.5*post + 0.2*rnormal()

di _n "=== wooldid, level() ==="
foreach L in 95 99.5 97.5 90.25 99.9 99.99 {
  capture noisily qui wooldid y cohort t ttre, att cluster(cohort) level(`L')
  di "   level(`L') -> rc " _rc
}
di _n "EXPECTED: level(99.9) and level(99.99) accepted, like 99.5 / 97.5 / 90.25."
di "OBSERVED: r(198) 'level() can have at most two digits after the decimal point'."
