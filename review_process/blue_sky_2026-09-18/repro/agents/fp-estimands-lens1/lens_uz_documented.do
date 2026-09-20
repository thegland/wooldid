*==============================================================================
* LENS-1 (DOCUMENTED-BEHAVIOUR) tiebreak probe for FP-CTZERO-ENGINE.
*
* Question this file answers:
*   (1) On a STRICTLY CONFORMING untreatedZero design (contreat exactly 0 outside
*       treated cohort-periods -- the premise sthlp:998 requires), does the default
*       engine recover the truth while customvce()/unconditionalse collapse?
*   (2) Is the cell collapse the axis, or the engine?  (DEFAULT vs DEFAULT
*       nofastcollapse: customvce/unconditionalse both turn the collapse off.)
*   (3) HETEROGENEOUS dose effects: is the explicit-arm answer "~0" (skeptic) or
*       "mean(true) - true_k for one dropped cell k" (my derivation)?  The second
*       is worse: a plausible nonzero wrong number rather than a suspicious null.
*   (4) Does the ringer reach the DEFAULT (ppmlhdfe) engine under poisson?
*       Structural: grep the verbose-printed command for the ringer term.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do lens_uz_documented.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

dis "=== BUILD ASSERTION ==="
which wooldid

*==============================================================================
* PART A -- homogeneous dose effect, conforming untreatedZero data
*==============================================================================
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ui = 100 + 20*runiform()          // UNIT-level intercept (before expand)
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post                     // EXACTLY 0 outside treated cohort-periods
gen double y  = ui + 0.5*t + 3.0*d2 + rnormal()*0.3
qui sum d2 if post==1
local ATT_A = 3*r(mean)
dis "TRUTH_A: ame = 3.000000   att = " %12.7f `ATT_A'
tempfile DA
save `DA'

* --- hand oracle, no wooldid: cell-specific dose slopes, NO pooled dose control
use `DA', clear
egen int cell = group(i t) if post==1
replace cell = 0 if post==0
qui reg y i.i i.t i.cell i.cell#c.d2
tempvar me
qui predict double `me', xb
* AME by hand = weighted mean of the cell dose slopes over treated-post rows
qui gen double slope = .
levelsof cell if cell>0, local(CL)
foreach c of local CL {
  qui replace slope = _b[`c'.cell#c.d2] if cell==`c'
}
qui sum slope if post==1
dis "ORACLE_A_no_pooled: ame = " %12.7f r(mean)
qui gen double attcontrib = slope*d2
qui sum attcontrib if post==1
dis "ORACLE_A_no_pooled: att = " %12.7f r(mean)

* --- the same design WITH a pooled dose control added (what the ringer does)
use `DA', clear
egen int cell = group(i t) if post==1
replace cell = 0 if post==0
qui cap noi reg y i.i i.t c.d2 i.cell i.cell#c.d2
dis "HAND_WITH_POOLED rc=" _rc
dis "HAND_WITH_POOLED pooled b(d2) = " %12.7f _b[d2]

* --- wooldid arms
foreach arm in "DEFAULT" "DEFAULT_NOFC" "CUSTOMVCE" "UNCONDSE" {
  use `DA', clear
  local o
  if "`arm'" == "DEFAULT_NOFC" local o nofastcollapse
  if "`arm'" == "CUSTOMVCE"    local o customvce(robust) verbose
  if "`arm'" == "UNCONDSE"     local o unconditionalse robust
  dis _n "##### A_`arm' #####"
  cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) `o'
  dis "A_RC_`arm'=" _rc
  if _rc==0 {
    mat b = e(b)
    dis "A_COLS_`arm': " "`: colnames b'"
    dis "A_RESULT_`arm': att=" %12.7f e(b)[1,1] "  ame=" %12.7f e(b)[1,2]
  }
}

* --- verbose on the DEFAULT arm too, to show what reghdfe is handed
use `DA', clear
dis _n "##### A_DEFAULT_VERBOSE #####"
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) verbose
dis "A_RC_DEFVERB=" _rc

*==============================================================================
* PART B -- HETEROGENEOUS dose effects (the discriminating test)
*   cohorts 1-3 slope 3.0, cohorts 4-6 slope 1.0
*==============================================================================
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ui = 100 + 20*runiform()
gen double beta = cond(i<=3, 3.0, 1.0)
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post
gen double y  = ui + 0.5*t + beta*d2 + rnormal()*0.3
qui sum beta if post==1
local AME_B = r(mean)
qui gen double attc = beta*d2
qui sum attc if post==1
local ATT_B = r(mean)
dis "TRUTH_B: ame = " %12.7f `AME_B' "   att = " %12.7f `ATT_B'
dis "TRUTH_B: if one cell of a beta=1 cohort is the dropped reference, predict ame ~ " %12.7f `AME_B'-1
dis "TRUTH_B: if one cell of a beta=3 cohort is the dropped reference, predict ame ~ " %12.7f `AME_B'-3
tempfile DB
save `DB'

foreach arm in "DEFAULT" "CUSTOMVCE" {
  use `DB', clear
  local o
  if "`arm'" == "CUSTOMVCE" local o customvce(robust) verbose
  dis _n "##### B_`arm' #####"
  cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) `o'
  dis "B_RC_`arm'=" _rc
  if _rc==0 dis "B_RESULT_`arm': att=" %12.7f e(b)[1,1] "  ame=" %12.7f e(b)[1,2]
}

*==============================================================================
* PART C -- poisson on the DEFAULT (ppmlhdfe) engine: does the ringer reach it?
*   xb-scale truth for the dose semi-elasticity is 0.30 exactly.
*==============================================================================
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ai = 1.0 + 0.4*runiform()
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post
gen double y  = exp(ai + 0.05*t + 0.30*d2 + rnormal()*0.05)
dis "TRUTH_C: xb-scale ame = 0.3000000"
tempfile DC
save `DC'

use `DC', clear
dis _n "##### C_POISSON_DEFAULT (verbose) #####"
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) poisson verbose
dis "C_RC_POIS=" _rc
if _rc==0 dis "C_RESULT_POIS: att=" %12.7f e(b)[1,1] "  ame=" %12.7f e(b)[1,2]

use `DC', clear
dis _n "##### C_LINEAR_DEFAULT (no poisson, reference that the design is fine) #####"
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero)
dis "C_RC_LIN=" _rc

dis _n "=== LENS1 PROBE COMPLETE ==="
