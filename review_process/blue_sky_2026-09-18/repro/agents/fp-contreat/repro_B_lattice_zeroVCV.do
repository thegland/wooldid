* ===========================================================================
* FINDING B (self-contained repro; run:  stata-mp -e do <thisfile>)
*
* CLAIM: with contreat() + lattice() on a realistic panel, wooldid's DEFAULT
* (reghdfe) estimation path returns an all-zero covariance matrix and loses
* every standard error, while the IDENTICAL model fitted by plain regress on
* the IDENTICAL data (customvce()/unconditionalse) returns finite ones.
* Same N, same df_r, same rss, same r2, same point estimates to 8 decimals.
* So the zero is an artifact of the absorbing engine, not of the design; it is
* neither the documented cell-containment geometry (sthlp:1076, which needs
* cluster() at or above i - this fires under bare `robust`) nor a
* lattice-construction defect (a user-supplied lattice variable does it too).
*
* The fixture is the realistic one: 12 states, 5-12 counties each, 20 years,
* unit-level lognormal dose with 10% time noise, staggered onsets, true dose
* effect exactly +0.4 on every treated post observation.
* ===========================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 121212
set sortseed 121212

clear
set obs 12
gen int st = _n
gen int nc = 5 + mod(_n*7,8)
gen double ttre = .
qui replace ttre =  8 if inlist(st,1,2,3)
qui replace ttre = 12 if inlist(st,4,5,6)
qui replace ttre = 15 if inlist(st,7,8)
expand nc
gen long unit = _n
expand 20
bysort unit: gen int year = _n
bysort unit (year): gen double ubase = exp(rnormal(0,0.6)) if _n==1
bysort unit (year): replace ubase = ubase[1]
gen double dose_it = ubase*(1 + 0.10*rnormal())
bysort unit (year): gen double ufe = rnormal(0,1) if _n==1
bysort unit (year): replace ufe = ufe[1]
gen double y_dose = 3 + ufe + 0.2*year + 0.5*dose_it + rnormal(0,1)
qui replace y_dose = y_dose + 0.4*dose_it + 0.3 if ttre<. & year>=ttre   /* TRUE ame = 0.4 */
qui gquantiles latu = dose_it , xtile nquantiles(3)   /* a user lattice with the same bins */
tempfile F
qui save `F'

foreach spec in ///
   "lattice(all 3) robust" ///
   "lattice(all 3) robust regtol(14)" ///
   "lattice(all 3) ccnoabsorb robust" ///
   "lattice(latu) robust" ///
   "lattice(all 3) customvce(robust)" ///
   "lattice(all 3) unconditionalse robust" ///
   "robust" {
  use `F', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) ///
      `spec' esprelength(2) espostlength(2)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "SPEC = `spec'"
    dis as result "   N=" e(N) "  df_r=" e(df_r) "  rss=" %12.6g e(rss) "  r2=" %8.6f e(r2)
    dis as result "   att      b=" %14.8f M[1,1] "  se=" %14.8g M[1,2]
    dis as result "   att_pre  b=" %14.8f M[2,1] "  se=" %14.8g M[2,2]
    dis as result "   ame      b=" %14.8f M[3,1] "  se=" %14.8g M[3,2] "   (TRUE ame = 0.4)"
    dis as result "   ame_pre  b=" %14.8f M[4,1] "  se=" %14.8g M[4,2]
  }
  else dis as error "SPEC = `spec'   RC=" _rc
}
dis "DONE"
