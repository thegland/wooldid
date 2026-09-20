*==============================================================================
* SELF-CONTAINED REPRO -- SK-SEMIELAS (independent confirmation of FP-SEMIELAS-FE)
* plus two extensions the seeded finding only inferred:
*   (X1) contreatelasticitytype(eyex) / (eydx) on the default engine, measured;
*   (X2) a strictly-positive-y design where reghdfe's _cons+Xb goes NEGATIVE and
*        the default engine dies r(3598) with a raw Mata traceback, while the
*        explicit-dummy engine returns the correct answer.
*
* CLAIM: wooldid's default engine is reghdfe with absorb(i.i i.t), so margins'
* prediction is _cons + Xb without the absorbed fixed effects.  eydx divides by
* that prediction, so the reported semi-elasticity measures the effect against the
* sample grand mean.  dydx differences the two predictions, so the fixed effects
* cancel and the level answers agree.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_semielas_own.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*==============================================================================
* PART 1 -- level agrees, semi-elasticity does not
*==============================================================================
capture program drop mkp
program define mkp
  syntax , HET(real)
  clear
  set obs 40
  gen int i = _n
  gen double alpha = 60 + `het'*(i-20)
  gen int ttre = .
  replace ttre = 4  if i <= 13
  replace ttre = 7  if i > 13 & i <= 26
  expand 10
  bysort i: gen int t = _n
  gen byte post = (t >= ttre) & ttre < .
  gen double y = alpha + `het'*0.5*t + 8*post + rnormal()*1
end

mkp , het(2)
tempfile H
save `H'

foreach q in LEVEL SEMI {
  foreach e in DEF EXP {
    use `H', clear
    local o1
    local o2
    if "`q'" == "SEMI" local o1 semielasticity
    if "`e'" == "EXP"  local o2 customvce(robust)
    dis _n "########## `q' / `e' ##########"
    cap noi wooldid y i t ttre, att `o1' `o2'
    dis "RC_`q'_`e'=" _rc
    if _rc==0 dis "RESULT `q'_`e' = " %21x e(b)[1,1] "  " %14.9f e(b)[1,1]
  }
}

* control arm: homogeneous unit and time levels -> the gap must collapse
mkp , het(0)
tempfile O
save `O'
foreach e in DEF EXP {
  use `O', clear
  local o2
  if "`e'" == "EXP" local o2 customvce(robust)
  dis _n "########## HOMOGENEOUS / SEMI / `e' ##########"
  cap noi wooldid y i t ttre, att semielasticity `o2'
  dis "RC_HOM_`e'=" _rc
  if _rc==0 dis "RESULT HOM_`e' = " %14.9f e(b)[1,1]
}

*==============================================================================
* PART 2 -- the mechanism, with no wooldid at all
*==============================================================================
use `H', clear
dis _n "########## MECHANISM ##########"
qui reghdfe y post, absorb(i.i i.t)
predict double xb_hdfe, xb
qui margins, eydx(post) noestimcheck
matrix Mh = r(b)
qui margins, dydx(post) noestimcheck
matrix Dh = r(b)
scalar bpost = _b[post]

qui reg y i.i i.t post
predict double xb_full, xb
qui margins, eydx(post) noestimcheck
matrix Me = r(b)
qui margins, dydx(post) noestimcheck
matrix De = r(b)

dis "reghdfe  dydx = " %14.9f Dh[1,1] "    explicit dydx = " %14.9f De[1,1] "    ratio = " %14.9f Dh[1,1]/De[1,1]
dis "reghdfe  eydx = " %14.9f Mh[1,1] "    explicit eydx = " %14.9f Me[1,1] "    ratio = " %14.9f Mh[1,1]/Me[1,1]
qui corr xb_hdfe xb_full
dis "corr(reghdfe prediction, true fitted value) = " %14.9f r(rho)
gen double rh = bpost/xb_hdfe
gen double rf = bpost/xb_full
qui sum rh, meanonly
dis "hand oracle mean(b/xb_hdfe) = " %14.9f r(mean) "   (must equal reghdfe eydx)"
qui sum rf, meanonly
dis "hand oracle mean(b/xb_full) = " %14.9f r(mean) "   (must equal explicit eydx)"
qui sum xb_hdfe
dis "SD(reghdfe prediction) = " %12.6f r(sd)
qui sum xb_full
dis "SD(true fitted value)  = " %12.6f r(sd)

*==============================================================================
* X1 -- contreatelasticitytype(eyex)/(eydx) on the default engine
*==============================================================================
clear
set obs 30
gen int i = _n
gen int ttre = .
replace ttre = 4 if i <= 10
replace ttre = 7 if i > 10 & i <= 20
gen double alpha = 40 + 3*i
expand 9
bysort i: gen int t = _n
expand 4
gen double dose = 1 + 3*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = alpha + 0.5*t + 2.0*dose*post + rnormal()*0.5
tempfile C
save `C'
foreach ct in dydx eydx eyex {
  foreach e in DEF EXP {
    use `C', clear
    local o2
    if "`e'" == "EXP" local o2 customvce(robust)
    dis _n "########## contreatelasticitytype(`ct') / `e' ##########"
    cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatelasticitytype(`ct') `o2'
    dis "RC_CT_`ct'_`e'=" _rc
    if _rc==0 dis "RESULT CT_`ct'_`e' : att=" %14.9f e(b)[1,1] "  ame=" %14.9f e(b)[1,2]
  }
}

*==============================================================================
* X2 -- strictly positive y, strictly positive true fitted values,
*       but reghdfe's _cons + Xb goes negative
*==============================================================================
clear
set obs 40
gen int i = _n
gen int ttre = .
replace ttre = 6 if i <= 16
gen double alpha = cond(i<=16, 26 + 0.5*i, 6 + 0.1*(i-16))
expand 10
bysort i: gen int t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = alpha - 20*post + rnormal()*0.3
qui sum y
dis _n "X2: min(y) = " %10.5f r(min) "   (strictly positive)"
tempfile X
save `X'
qui reghdfe y post, absorb(i.i i.t)
predict double xbh, xb
qui sum xbh
dis "X2: reghdfe prediction  min = " %10.5f r(min)
qui count if xbh <= 0
dis "X2: reghdfe prediction <= 0 on " r(N) " of 400 rows"
qui reg y i.i i.t post
predict double xbf, xb
qui sum xbf
dis "X2: true fitted value   min = " %10.5f r(min)
qui count if xbf <= 0
dis "X2: true fitted value  <= 0 on " r(N) " of 400 rows"

use `X', clear
dis _n "########## X2 semielasticity / default reghdfe ##########"
cap noi wooldid y i t ttre, att semielasticity
dis "RC_X2_DEF=" _rc
if _rc==0 dis "RESULT X2_DEF = " %14.9f e(b)[1,1]

use `X', clear
dis _n "########## X2 semielasticity / customvce(robust) ##########"
cap noi wooldid y i t ttre, att semielasticity customvce(robust)
dis "RC_X2_EXP=" _rc
if _rc==0 dis "RESULT X2_EXP = " %14.9f e(b)[1,1]
dis "X2 hand check: treated-post y ~ alpha-20, pre ~ alpha, alpha in 26.5..34"
dis "  so the correct semi-elasticity is mean ln((alpha-20)/alpha) ~ -1.1"

dis _n "REPRO_SEMIELAS_DONE"
