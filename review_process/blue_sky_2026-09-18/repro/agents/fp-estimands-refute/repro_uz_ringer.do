*==============================================================================
* SELF-CONTAINED REPRO -- SK-UZ-RINGER
*
* contreatcontrols(untreatedZero) combined with ANY of unconditionalse,
* customvce() or poisson silently installs a POOLED dose control that
* untreatedZero exists to exclude, and margins never switches it off.  A true
* dose effect of 3.0 is reported as ~0.
*
* ado:1233   the DF-1 "ringer" head is installed whenever
*              strpos(contreatcontrols,"basic") | unconditionalse | customvce() | poisson
*            -- the test does NOT exclude untreatedZero.
* ado:1300   contreatcontrolinserterFE = ccringerhead ccringertail
* ado:1313   the untreatedZero branch rewrites `contreatinj' only; it never clears
*            contreatcontrolinserterFE, which ado:1485 (reg/poisson) and, through
*            contreatcontrolinserterABS, ado:1462 (ppmlhdfe) then carry.
* The DF-1 comment justifies the ringer as "a numerical twin of columns the
* margins-visible block already carries through its ## chain (c.`contreat', ...)"
* -- true under basic, FALSE under untreatedZero, whose whole definition
* (sthlp:1033) is "by not including a g() function at all".
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_uz_ringer.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31415
set sortseed 31415
which wooldid

*-------------------------------------------------------------- clean DGP
* 9 cohorts (6 treated, 3 never treated), 7 periods, 14 records per cohort-period.
* UNIT-LEVEL intercept, sigma 0.3.  Dose is exactly 0 outside treated cohort-periods,
* which is untreatedZero's own premise.  True dose effect = 3.0.
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ui = 100 + 20*runiform()
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post
gen double y  = ui + 0.5*t + 3.0*d2 + rnormal()*0.3
qui sum d2 if post==1
local TATT = 3*r(mean)
dis "TRUTH:  ame = 3.000000     att = " %10.6f `TATT'
tempfile D
save `D'

*-------------------------------------------------------------- wooldid arms
foreach arm in "DEFAULT" "CUSTOMVCE" "UNCONDSE" {
  use `D', clear
  local o
  if "`arm'" == "CUSTOMVCE" local o customvce(robust)
  if "`arm'" == "UNCONDSE"  local o unconditionalse robust
  dis _n "########## untreatedZero  `arm' ##########"
  cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) `o'
  dis "RC_`arm'=" _rc
  if _rc==0 dis "RESULT `arm':  att=" %12.7f e(b)[1,1] "   ame=" %12.7f e(b)[1,2] "   (truth 4.25 / 3.00)"
}

* control: contreatcontrols(basic) on a dose that is defined everywhere -- the two
* engines agree to the last bit there, so the split is untreatedZero-specific
foreach arm in "DEFAULT" "CUSTOMVCE" {
  use `D', clear
  local o
  if "`arm'" == "CUSTOMVCE" local o customvce(robust)
  dis _n "########## basic (control arm)  `arm' ##########"
  cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) `o'
  dis "RC_BASIC_`arm'=" _rc
  if _rc==0 dis "RESULT basic `arm':  att=" %12.7f e(b)[1,1] "   ame=" %12.7f e(b)[1,2]
}

*-------------------------------------------------------------- poisson arm
* ado:1233 also fires on poisson, and there the ringer goes inside absorb(),
* so the DEFAULT engine is hit too.
use `D', clear
gen double mu = exp(1.0 + 0.10*t + 0.40*d2)
gen int ycount = rpoisson(mu)
qui sum d2 if post==1
dis _n "POISSON TRUTH (xb scale): ame = 0.400000   att = " %10.6f 0.40*r(mean)
tempfile P
save `P'
dis _n "########## untreatedZero  POISSON (default engine) ##########"
cap noi wooldid ycount i t ttre, att contreat(d2) contreatcontrols(untreatedZero) poisson
dis "RC_POIS=" _rc
if _rc==0 dis "RESULT POISSON:  att=" %12.7f e(b)[1,1] "   ame=" %12.7f e(b)[1,2]

*-------------------------------------------------------------- mechanism, wooldid-free
* Rebuild wooldid's untreatedZero design by hand and show that adding the pooled
* dose term is exactly what destroys the margins contrast.
use `D', clear
gen byte TF = (ttre < .)
gen int relyr = t - ttre
egen long cell = group(i relyr) if TF==1 & relyr != -1
replace cell = 0 if missing(cell)
tempfile M
save `M'

dis _n "########## MECHANISM (no wooldid) ##########"
qui reg y i.i i.t i.cell i.cell#c.d2
predict double p1, xb
preserve
  replace cell = 0
  predict double p0, xb
  gen double eff = p1 - p0
  qui sum eff if post==1
  dis "intended model  (no pooled dose term):  mean effect on treated-post = " %10.6f r(mean)
restore

use `M', clear
qui reg y i.i i.t c.d2 i.cell i.cell#c.d2
dis "pooled ringer coefficient on d2        = " %10.6f _b[d2]
predict double q1, xb
preserve
  replace cell = 0
  predict double q0, xb
  gen double eff = q1 - q0
  qui sum eff if post==1
  dis "shipped explicit model (pooled ringer): mean effect on treated-post = " %10.6f r(mean)
restore

dis _n "EXPECTED: all arms recover ~3.0 (ame) / ~4.25 (att)."
dis "OBSERVED: the reghdfe default arm does; customvce, unconditionalse and poisson do not."
dis "REPRO_UZ_DONE"
