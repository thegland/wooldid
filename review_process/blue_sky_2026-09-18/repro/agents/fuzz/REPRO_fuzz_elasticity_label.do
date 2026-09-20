*=======================================================================
* SELF-CONTAINED REPRO -- fuzz agent
* CLAIM: with a BINARY treatment (no contreat()), the option
*        contreatelasticitytype(eyex) (or eydx) is silently INERT --
*        the estimates, standard errors, t, p and interval come back
*        bit-identical to the plain ATT -- but the results table is
*        relabelled "Main Results (Full Estimation Sample): Marginal
*        Treatment Effects Presented as Elasticities".
*        A reader takes a level ATT away as an elasticity.
*        wooldid refuses every other contreat companion in that
*        position (lattice(), contreatcontrols(), contreatwithin,
*        makecfxplots all error r(198) "requires contreat()"), so the
*        missing screen here is an inconsistency as well as a hazard.
*        Contrast semielasticity, which in the same position DOES
*        transform the number and labels it correctly.
* RUN:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_fuzz_elasticity_label.do
*=======================================================================
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

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
gen double g_t  = 0.10*(t-2000)
gen byte post   = !missing(ttre) & t >= ttre
gen double ypos = 50 + a_i + g_t + 1.5*post + 0.2*rnormal()

di _n "=== (1) plain ATT ==="
wooldid ypos cohort t ttre, att cluster(cohort)
matrix BASE = e(wdidmainresults)

di _n "=== (2) the SAME command plus contreatelasticitytype(eyex) ==="
di "    -- watch the header, then compare the numbers"
wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(eyex)
matrix EY = e(wdidmainresults)

di _n "=== comparison ==="
di "plain att              = " %21x BASE[1,1] "   (" BASE[1,1] ")"
di "'elasticity' att       = " %21x EY[1,1]   "   (" EY[1,1]   ")"
di "bit-identical?           " (BASE[1,1]==EY[1,1])
di "plain se               = " BASE[1,2]
di "'elasticity' se        = " EY[1,2]
di "bit-identical?           " (BASE[1,2]==EY[1,2])

di _n "=== (3) contrast: semielasticity DOES transform, and labels it correctly ==="
wooldid ypos cohort t ttre, att cluster(cohort) semielasticity

di _n "=== (4) contrast: every other contreat companion IS refused in this position ==="
foreach o in "lattice(all 3)" "contreatcontrols(basic)" "contreatwithin" "makecfxplots" {
  capture noisily wooldid ypos cohort t ttre, att cluster(cohort) `o'
  di "   `o' -> rc " _rc
}
capture noisily wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(eyex)
di "   contreatelasticitytype(eyex) -> rc " _rc  "   <-- NOT refused"

di _n "EXPECTED: either a refusal ('contreatelasticitytype() has no effect without contreat()'),"
di "          or a table header that does not claim elasticities."
di "OBSERVED: rc 0, numbers identical to the plain ATT, header says 'Presented as Elasticities'."
