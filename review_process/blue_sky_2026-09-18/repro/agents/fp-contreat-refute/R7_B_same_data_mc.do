* ============================================================================
* R7 : (a) the lattice all-zero-V, reghdfe vs regress ON THE SAME DATASET, with
*          an estimability check on the ame;
*      (b) a 12-draw Monte Carlo that separates VARIANCE from BIAS in the
*          "45x the truth" claim (fp-contreat-B2), and compares the reported
*          conditional robust SE against the realized sampling SD of the ame.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R7_B_same_data_mc.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 13579
set sortseed 13579

capture program drop mkdata
program define mkdata
  args tvar
  clear
  set obs 12
  gen int st = _n
  gen int ncty = 5 + mod(_n*7, 8)
  expand ncty
  gen long cty = _n
  expand 10
  bysort cty: gen int year = _n
  gen double ttre = .
  qui replace ttre = 4 if st <= 3
  qui replace ttre = 6 if st > 3 & st <= 6
  qui replace ttre = 8 if st > 6 & st <= 8
  gen byte post = (ttre < . & year >= ttre)
  bysort cty (year): gen double dbase = 2*runiform() if _n == 1
  bysort cty (year): replace dbase = dbase[1]
  gen double dose = dbase + `tvar'*rnormal(0, 0.5)
  gen double aunit = 0
  bysort cty (year): replace aunit = 3*runiform() if _n == 1
  bysort cty (year): replace aunit = aunit[1]
  gen double y = 5 + aunit + 0.25*year + 0.7*dose + 0.4*dose*post + 0.5*post + rnormal(0,1)
end

* ---------------- (a) SAME dataset, both engines --------------------------
forval rep = 1/6 {
  mkdata 1
  tempfile F
  qui save `F'
  local SED = 0
  use `F', clear
  cap noi qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust
  local rcA = _rc
  local ameA = .
  local seA  = .
  if `rcA' == 0 {
    mat M = e(wdidmainresults)
    local ameA = M[2,1]
    local seA  = M[2,2]
  }
  use `F', clear
  cap noi qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) customvce(robust)
  local rcB = _rc
  local ameB = .
  local seB  = .
  if `rcB' == 0 {
    mat M = e(wdidmainresults)
    local ameB = M[2,1]
    local seB  = M[2,2]
  }
  * outside wooldid, same data
  use `F', clear
  qui egen cell = group(st year) if post == 1
  qui replace cell = 0 if cell >= .
  qui xtile lat = dose , nq(3)
  qui gen double ringer = dose
  qui gen byte nullv = 1
  qui reghdfe y ib0.cell##i.lat##c.dose , absorb(i.st i.year i.nullv##c.ringer i.lat##c.ringer) vce(robust)
  mat VH = e(V)
  mata: st_local("mvh", strofreal(max(abs(st_matrix("VH")))))
  local rkh = e(rank)
  local rssh = e(rss)
  qui margins , dydx(dose) subpop(if post == 1) noestimcheck nose
  mat r = r(b)
  local amH = r[1,1]
  qui reg y i.st i.year i.nullv##c.ringer i.lat##c.ringer ib0.cell##i.lat##c.dose , vce(robust)
  mat VR = e(V)
  mata: st_local("mvr", strofreal(max(abs(st_matrix("VR")))))
  local rkr = e(rank)
  local rssr = e(rss)
  qui margins , dydx(dose) subpop(if post == 1) noestimcheck nose
  mat r = r(b)
  local amR = r[1,1]
  dis as result "rep`rep' | WD def ame=" %12.6f `ameA' " se=" %10.6f `seA' ///
     " | WD cvce ame=" %12.6f `ameB' " se=" %10.6f `seB'
  dis as result "      out: reghdfe max|V|=`mvh' rank=`rkh' rss=" %12.4f `rssh' " margins-ame=" %12.6f `amH'
  dis as result "      out: regress max|V|=`mvr' rank=`rkr' rss=" %12.4f `rssr' " margins-ame=" %12.6f `amR'
}

* ---------------- (b) Monte Carlo: variance vs bias -----------------------
tempname MC
postfile `MC' double ame double se double amenolat double senolat using mcres , replace
forval rep = 1/12 {
  mkdata 1
  tempfile G
  qui save `G'
  local a = .
  local s = .
  local a2 = .
  local s2 = .
  use `G', clear
  cap qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) customvce(robust)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    local a = M[2,1]
    local s = M[2,2]
  }
  use `G', clear
  cap qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) customvce(robust)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    local a2 = M[2,1]
    local s2 = M[2,2]
  }
  dis as result "MC rep`rep': lat ame=" %12.6f `a' " se=" %10.6f `s' "   nolat ame=" %12.6f `a2' " se=" %10.6f `s2'
  post `MC' (`a') (`s') (`a2') (`s2')
}
postclose `MC'
use mcres , clear
qui sum ame
dis as result "MC LATTICE ame: mean=" %12.6f r(mean) " SD=" %12.6f r(sd) " min=" %12.6f r(min) " max=" %12.6f r(max) "  (TRUE 0.400)"
qui sum se
dis as result "MC LATTICE reported se: mean=" %12.6f r(mean) " min=" %12.6f r(min) " max=" %12.6f r(max)
qui sum amenolat
dis as result "MC NOLATTICE ame: mean=" %12.6f r(mean) " SD=" %12.6f r(sd)
qui sum senolat
dis as result "MC NOLATTICE reported se: mean=" %12.6f r(mean)

dis "DONE_R7"
