*==============================================================================
* SELF-CONTAINED REPRO -- FP-CTZERO-ENGINE  (UNCERTAIN; reported as a lead)
*
* CLAIM: on a low-noise, well-identified continuous-treatment design with
* contreatcontrols(untreatedZero), the DEFAULT engine (reghdfe) and the
* explicit-dummy engine (customvce) return main_att and main_ame of OPPOSITE SIGN
* from the same data with the same R2 -- and only the explicit-dummy engine prints
* wooldid's multicollinearity caution.  A user on the default path gets a
* sign-flipped estimate at p < 0.001 with no warning at all.
*
* With contreatcontrols(basic) on the SAME data the two engines agree exactly
* (dydx ratio 1.00000), so this is specific to untreatedZero, not to the engines.
*
* I could not determine WHICH side is right; that needs someone who knows the
* untreatedZero model formulation (ado:1313).  Reported as a lead, not a verdict.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_contreat_engine_split.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260925
set sortseed 20260925
which wooldid

clear
set obs 9
gen int i = _n                       // i is a COHORT: contreat needs in-cell dose variation
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
expand 7
bysort i: gen int t = _n
expand 14                            // 14 records per cohort-period
gen double ui   = 100 + 20*runiform()
gen double dose = 0.2 + 2.5*runiform()
gen byte  post  = (t >= ttre) & ttre < .
gen double dz   = dose*post
gen double y    = ui + 0.5*t + 3.0*dz + rnormal()*0.3     // TRUE dose effect = 3.0
tempfile d
save `d'

foreach cc in untreatedZero basic {
  foreach eng in "" "customvce(robust)" {
    use `d', clear
    if "`cc'" == "untreatedZero" gen double d2 = dz
    else                         gen double d2 = dose
    local tag = cond("`eng'"=="", "DEFAULT(reghdfe)", "EXPLICIT(customvce)")
    dis _n "########## contreatcontrols(`cc')  engine=`tag' ##########"
    cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(`cc') `eng'
    dis "RC=" _rc
    if _rc == 0 {
      dis "`cc' `tag'  att=" %12.7f e(b)[1,1] "   ame=" %12.7f e(b)[1,2] ///
          "   R2=" %9.6f e(r2) "   N=" e(N)
    }
  }
}
dis _n "EXPECTED: under basic the two engines agree to the last bit."
dis "OBSERVED under untreatedZero: opposite signs, identical R2, and the"
dis "multicollinearity caution appears on the EXPLICIT arm only."
dis "REPRO_DONE"
