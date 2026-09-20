* ============================================================================
* R6_B : SKEPTIC's OWN repro of seeded findings fp-contreat-B1 / B2 / C
*   B1 claim: lattice()+contreat() on the default reghdfe path returns an
*             all-zero e(V) (so every SE is lost) while plain regress on the
*             SAME fit returns finite SEs -> "an artifact of the engine".
*   B2 claim: on such a fixture the ame is ~45x the truth with no usable warning
*   C  claim: the all-zero-VCV cautions blame clustered SEs that were never
*             requested and point at a caution that never printed.
*
* My fixture is built independently: 12 states x (5-12 counties) x 10 years,
* county-level dose, staggered onsets, i = state.  A knob `tvar' controls how
* much WITHIN-county (over time) dose variation there is: the advisor's
* discriminating variable.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R6_B_lattice.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242

capture program drop mkdata
program define mkdata
  args tvar
  clear
  set obs 12
  gen int st = _n
  gen int ncty = 5 + mod(_n*7, 8)          /* 5..12 counties per state */
  expand ncty
  gen long cty = _n
  expand 10
  bysort cty: gen int year = _n
  gen double ttre = .
  qui replace ttre = 4 if st <= 3
  qui replace ttre = 6 if st > 3 & st <= 6
  qui replace ttre = 8 if st > 6 & st <= 8
  * st 9..12 never treated
  gen byte post = (ttre < . & year >= ttre)
  * county-level dose, plus `tvar' of within-county time variation
  bysort cty (year): gen double dbase = 2*runiform() if _n == 1
  bysort cty (year): replace dbase = dbase[1]
  gen double dose = dbase + `tvar'*rnormal(0, 0.5)
  gen double aunit = 0
  bysort cty (year): replace aunit = 3*runiform() if _n == 1
  bysort cty (year): replace aunit = aunit[1]
  * TRUE dose effect is exactly +0.4 on every treated post observation, linear,
  * with no bin structure at all
  gen double y = 5 + aunit + 0.25*year + 0.7*dose + 0.4*dose*post + 0.5*post + rnormal(0,1)
end

foreach tv in 0 1 3 {
  mkdata `tv'
  tempfile F
  qui save `F'
  qui count
  local NN = r(N)
  dis as result "===== tvar=`tv'  N=`NN'  (TRUE dose effect on treated post rows = 0.400) ====="
  * how much within-(county) dose variation is there?
  qui reghdfe dose , absorb(i.cty)
  dis as result "  within-county dose rmse = " %10.6f e(rmse)

  * no lattice baseline
  use `F', clear
  cap noi qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) robust
  dis "RC_nolat_`tv'=" _rc
  if _rc == 0 {
    mat M = e(wdidmainresults)
    mat V = e(V)
    mata: st_local("mv", strofreal(max(abs(st_matrix("V")))))
    dis as result "  NOLAT robust     ame=" %16.8f M[2,1] " se=" %14.8f M[2,2] " max|e(V)|=`mv'"
  }

  foreach spec in "robust" "customvce(robust)" "unconditionalse robust" "cluster(cty)" {
    use `F', clear
    cap noi qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) `spec'
    dis "RC_lat_`tv'_`spec'=" _rc
    if _rc == 0 {
      mat M = e(wdidmainresults)
      mat V = e(V)
      mata: st_local("mv", strofreal(max(abs(st_matrix("V")))))
      dis as result "  LAT3 `spec' | ame=" %16.8f M[2,1] " se=" %14.8f M[2,2] ///
          " att=" %14.8f M[1,1] " attse=" %12.8f M[1,2] " max|e(V)|=`mv'"
    }
  }
}

* ---------- the caution text under BARE robust (finding C) ----------------
dis as result "===== NOISY leg: lattice(all 3) + robust, no cluster() anywhere ====="
mkdata 0
cap noi wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust
dis "RC_NOISY=" _rc

* ---------- outside wooldid: is the zero the engine or the design? --------
dis as result "===== outside wooldid, tvar=0 ====="
mkdata 0
qui egen cell = group(st year) if post == 1
qui replace cell = 0 if cell >= .
* dose-quantile lattice, 3 bins, built the plain way
qui xtile lat = dose , nq(3)
gen double ringer = dose
gen byte nullv = 1
qui reghdfe y ib0.cell##i.lat##c.dose , absorb(i.st i.year i.nullv##c.ringer i.lat##c.ringer) vce(robust)
mat VH = e(V)
mata: st_local("mvh", strofreal(max(abs(st_matrix("VH")))))
dis as result "  reghdfe  max|e(V)|=`mvh'  df_r=" e(df_r) "  rss=" %14.4f e(rss) "  rank=" e(rank) " N=" e(N)
qui reg y i.st i.year i.nullv##c.ringer i.lat##c.ringer ib0.cell##i.lat##c.dose , vce(robust)
mat VR = e(V)
mata: st_local("mvr", strofreal(max(abs(st_matrix("VR")))))
dis as result "  regress  max|e(V)|=`mvr'  df_r=" e(df_r) "  rss=" %14.4f e(rss) "  rank=" e(rank) " N=" e(N)

dis as result "===== outside wooldid, tvar=3 ====="
mkdata 3
qui egen cell = group(st year) if post == 1
qui replace cell = 0 if cell >= .
qui xtile lat = dose , nq(3)
gen double ringer = dose
gen byte nullv = 1
qui reghdfe y ib0.cell##i.lat##c.dose , absorb(i.st i.year i.nullv##c.ringer i.lat##c.ringer) vce(robust)
mat VH = e(V)
mata: st_local("mvh", strofreal(max(abs(st_matrix("VH")))))
dis as result "  reghdfe  max|e(V)|=`mvh'  df_r=" e(df_r) "  rss=" %14.4f e(rss) "  rank=" e(rank) " N=" e(N)
qui reg y i.st i.year i.nullv##c.ringer i.lat##c.ringer ib0.cell##i.lat##c.dose , vce(robust)
mat VR = e(V)
mata: st_local("mvr", strofreal(max(abs(st_matrix("VR")))))
dis as result "  regress  max|e(V)|=`mvr'  df_r=" e(df_r) "  rss=" %14.4f e(rss) "  rank=" e(rank) " N=" e(N)

dis "DONE_R6"
