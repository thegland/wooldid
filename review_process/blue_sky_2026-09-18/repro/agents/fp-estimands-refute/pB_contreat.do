*==============================================================================
* SKEPTIC PROBE B -- FP-CTZERO-ENGINE, rebuilt from the claim on CLEAN data.
*
* The seeded repro calls its design "low-noise, well-identified".  It is not:
* its unit intercept `ui' is generated AFTER expand, so it is per-ROW noise with
* SD ~5.8, swamping the stated rnormal()*0.3.  Here `ui' is genuinely unit-level.
*
* Also: the DGP truth is computable.  dose ~ U(0.2,2.7) so E[dose]=1.45 and with a
* true dose effect of 3.0 and no pure onset shift,
*     TRUE ame = 3.00      TRUE att = 3.0 * E[dose | treated-post] ~= 4.35
* so we can ask WHICH engine is deviant rather than only that they differ.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do pB_contreat.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4471
set sortseed 4471
which wooldid

*---------------------------------------------------------------------------
* B1: CLEAN design -- unit-level intercept, genuine sigma 0.3
*---------------------------------------------------------------------------
capture program drop mkct
program define mkct
  syntax , ROWNOISE(real) [NPER(int 7) NREC(int 14)]
  clear
  set obs 9
  gen int i = _n
  gen int ttre = .
  replace ttre = 3 if i <= 3
  replace ttre = 5 if i > 3 & i <= 6
  gen double ui_unit = 100 + 20*runiform()          // UNIT-LEVEL intercept
  expand `nper'
  bysort i: gen int t = _n
  expand `nrec'
  gen double ui = ui_unit + `rownoise'*rnormal()
  gen double dose = 0.2 + 2.5*runiform()
  gen byte post = (t >= ttre) & ttre < .
  gen double dz = dose*post
  gen double y = ui + 0.5*t + 3.0*dz + rnormal()*0.3
end

capture program drop rct
program define rct
  args tag
  dis "RES `tag'  att=" %13.7f e(b)[1,1] "  ame=" %13.7f e(b)[1,2] ///
      "  R2=" %9.6f e(r2) "  N=" e(N)
end

foreach rn in 0 {
  mkct , rownoise(`rn')
  qui sum dz if post==1
  dis "TRUTH CHECK: E[dose|post] = " %9.6f r(mean) "  => TRUE att = " %9.6f 3*r(mean) "  TRUE ame = 3.000000"
  tempfile dclean
  save `dclean'

  foreach cc in untreatedZero basic {
    foreach eng in "" "customvce(robust)" {
      use `dclean', clear
      if "`cc'" == "untreatedZero" gen double d2 = dz
      else                         gen double d2 = dose
      local tag = cond("`eng'"=="", "DEFAULT", "EXPLICIT")
      dis _n "##### CLEAN cc=`cc' eng=`tag' #####"
      cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(`cc') `eng'
      dis "RC_CLEAN_`cc'_`tag'=" _rc
      if _rc==0 rct "CLEAN `cc' `tag'"
    }
  }
}

*---------------------------------------------------------------------------
* B2: same but with the seeded harness's per-row noise, to show the artifact
*---------------------------------------------------------------------------
mkct , rownoise(5.77)
tempfile dnoisy
save `dnoisy'
foreach cc in untreatedZero {
  foreach eng in "" "customvce(robust)" {
    use `dnoisy', clear
    gen double d2 = dz
    local tag = cond("`eng'"=="", "DEFAULT", "EXPLICIT")
    dis _n "##### NOISY cc=`cc' eng=`tag' #####"
    cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(`cc') `eng'
    dis "RC_NOISY_`cc'_`tag'=" _rc
    if _rc==0 rct "NOISY `cc' `tag'"
  }
}

*---------------------------------------------------------------------------
* B3: wooldid-free -- does reghdfe keep an omitted (structurally all-zero)
*     interaction slope in e(b) colnames with an "o." flag, the way reg does?
*     wooldid_ccscreen (ado:11818) scans `: colnames e(b)' for "o.".  If reghdfe
*     strips such names, the multicollinearity caution is DEAD on the default path.
*---------------------------------------------------------------------------
dis _n "##### B3: e(b) colname survey #####"
clear
set obs 800
gen int g = mod(_n,8)+1
gen int tt = mod(floor(_n/8),5)+1
gen double x = runiform()
gen double z = x*(g>4)                    // structurally 0 for g = 1..4
gen double yy = 2*g + tt + 1.5*z + rnormal()

qui reg yy i.g i.tt i.g#c.z
local bn : colnames e(b)
dis "REG colnames: `bn'"
dis "REG has o. :  " (strpos("`bn'","o.")!=0)

qui reghdfe yy i.g#c.z, absorb(i.g i.tt)
local bh : colnames e(b)
dis "RGHDFE colnames: `bh'"
dis "RGHDFE has o. : " (strpos("`bh'","o.")!=0)

* and a genuinely collinear (not structurally zero) duplicate regressor
gen double zdup = z
qui reg yy i.g i.tt z zdup
local b2 : colnames e(b)
dis "REG dup colnames: `b2'"
qui reghdfe yy z zdup, absorb(i.g i.tt)
local b3 : colnames e(b)
dis "RGHDFE dup colnames: `b3'"

dis _n "PROBE_B_DONE"
