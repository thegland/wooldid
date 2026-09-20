*==============================================================================
* SKEPTIC PROBE A -- FP-SEMIELAS-FE, built from the claim, not from their script.
*
* CLAIM UNDER TEST: on the DEFAULT engine (reghdfe, absorb(i.i i.t)) the absorbed
* fixed effects are not in e(b), so margins' prediction is _cons + Xb and
* margins, eydx(targetflag) divides by a number that is NOT the fitted outcome.
* The same model with explicit i.i/i.t dummies (customvce/unconditionalse) divides
* by the true fitted value.  Falsifiable test:  level (dydx) answers must agree
* bit-for-bit across engines while eydx answers differ, and the gap must vanish
* when unit/time levels are made homogeneous.
*
* My own DGP: unit-level (NOT per-row) heterogeneous intercepts.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do pA_semielas.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103

which wooldid

*---------------------------------------------------------------- DGP builder
capture program drop mkpanel
program define mkpanel
  syntax , HET(real) [SHIFT(real 0)]
  clear
  set obs 40
  gen int i = _n
  gen double alpha = 60 + `het'*(i-20)          // unit-level level heterogeneity
  gen int ttre = .
  replace ttre = 4  if i <= 13
  replace ttre = 7  if i > 13 & i <= 26
  expand 10
  bysort i: gen int t = _n
  gen byte post = (t >= ttre) & ttre < .
  gen double y = alpha + `het'*0.5*t + 8*post + rnormal()*1 + `shift'
  * keep y strictly positive under both arms
end

*---------------------------------------------------------------- reporter
capture program drop rep1
program define rep1
  args tag
  local k = colsof(e(b))
  local nm : colnames e(b)
  dis "RESULT `tag' : b1=" %21x e(b)[1,1] "  b1dec=" %14.9f e(b)[1,1] "  N=" e(N) " cols=`k' names=`nm'"
end

*==============================================================================
* ARM 1: heterogeneous unit levels
*==============================================================================
mkpanel , het(2)
sum y
tempfile hetd
save `hetd'

dis _n "########## HET / LEVEL / default reghdfe ##########"
cap noi wooldid y i t ttre, att
dis "RC_H_LVL_DEF=" _rc
if _rc==0 rep1 H_LVL_DEF

use `hetd', clear
dis _n "########## HET / LEVEL / customvce(robust) ##########"
cap noi wooldid y i t ttre, att customvce(robust)
dis "RC_H_LVL_EXP=" _rc
if _rc==0 rep1 H_LVL_EXP

use `hetd', clear
dis _n "########## HET / SEMIELAS / default reghdfe ##########"
cap noi wooldid y i t ttre, att semielasticity
dis "RC_H_SEM_DEF=" _rc
if _rc==0 rep1 H_SEM_DEF

use `hetd', clear
dis _n "########## HET / SEMIELAS / customvce(robust) ##########"
cap noi wooldid y i t ttre, att semielasticity customvce(robust)
dis "RC_H_SEM_EXP=" _rc
if _rc==0 rep1 H_SEM_EXP

use `hetd', clear
dis _n "########## HET / SEMIELAS / unconditionalse ##########"
cap noi wooldid y i t ttre, att semielasticity unconditionalse
dis "RC_H_SEM_UNC=" _rc
if _rc==0 rep1 H_SEM_UNC

*==============================================================================
* ARM 2 (CONTROL): homogeneous unit and time levels -> gap should vanish
*==============================================================================
mkpanel , het(0)
tempfile homd
save `homd'

dis _n "########## HOM / SEMIELAS / default reghdfe ##########"
cap noi wooldid y i t ttre, att semielasticity
dis "RC_O_SEM_DEF=" _rc
if _rc==0 rep1 O_SEM_DEF

use `homd', clear
dis _n "########## HOM / SEMIELAS / customvce(robust) ##########"
cap noi wooldid y i t ttre, att semielasticity customvce(robust)
dis "RC_O_SEM_EXP=" _rc
if _rc==0 rep1 O_SEM_EXP

*==============================================================================
* MECHANISM DEMO -- NO wooldid.  Same claim, plain Stata.
*==============================================================================
use `hetd', clear
dis _n "########## MECHANISM: reghdfe vs explicit dummies ##########"
qui reghdfe y post, absorb(i.i i.t)
dis "reghdfe b(post) = " %21x _b[post] "   " %14.9f _b[post]
predict double xb_hdfe, xb
qui margins, eydx(post) noestimcheck
matrix Mh = r(b)
dis "reghdfe  margins eydx(post) = " %21x Mh[1,1] "   " %14.9f Mh[1,1]
qui margins, dydx(post) noestimcheck
matrix Dh = r(b)
dis "reghdfe  margins dydx(post) = " %14.9f Dh[1,1]

qui reg y i.i i.t post
dis "reg      b(post) = " %21x _b[post] "   " %14.9f _b[post]
predict double xb_full, xb
qui margins, eydx(post) noestimcheck
matrix Me = r(b)
dis "explicit margins eydx(post) = " %21x Me[1,1] "   " %14.9f Me[1,1]
qui margins, dydx(post) noestimcheck
matrix De = r(b)
dis "explicit margins dydx(post) = " %14.9f De[1,1]

dis "EYDX RATIO default/explicit = " %14.9f Mh[1,1]/Me[1,1]
dis "DYDX RATIO default/explicit = " %14.9f Dh[1,1]/De[1,1]

sum xb_hdfe xb_full y
corr xb_hdfe xb_full
dis "corr(xb_hdfe, xb_full) = " %14.9f r(rho)

* hand oracles: for a CONTINUOUS x, margins eydx = mean over sample of b/yhat
qui reghdfe y post, absorb(i.i i.t)
scalar bh = _b[post]
gen double r_hdfe = bh/xb_hdfe
gen double r_full = bh/xb_full
sum r_hdfe, meanonly
dis "HAND mean(b/xb_hdfe) = " %14.9f r(mean)
sum r_full, meanonly
dis "HAND mean(b/xb_full) = " %14.9f r(mean)

dis _n "PROBE_A_DONE"
