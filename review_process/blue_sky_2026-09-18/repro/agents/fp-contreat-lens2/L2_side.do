* ============================================================================
* LENS 2 side checks: what does the user SEE, and does the zero leak into a
* silently wrong NUMBER anywhere downstream?
*  S1 what wooldid prints (caution text) on the default zero-V path
*  S2 what e(V) actually contains (zeros -> downstream t=inf/p=0 would be a
*     silent wrong number; missings -> loud)
*  S3 a downstream test/lincom on the posted e(b)/e(V)
*  S4 boottstat on the zero-V spec: studentized t = (b_r - b)/0 -- number or refusal?
*  S5 incidence: does the zero survive cluster(cty) (FINER than i -- the spec a
*     real applied economist would actually write), and lattice(all 2)?
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do L2_side.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242

which wooldid

clear
set obs 12
gen int st = _n
gen int ncty = 5 + mod(_n*7, 8)
expand ncty
gen long cty = _n
gen byte reg3 = 1 + mod(_n,3)
expand 10
bysort cty: gen int year = _n
gen double ttre = .
qui replace ttre = 4 if st <= 3
qui replace ttre = 6 if st > 3 & st <= 6
qui replace ttre = 8 if st > 6 & st <= 8
gen byte post = (ttre < . & year >= ttre)
bysort cty (year): gen double dbase = 2*runiform() if _n == 1
bysort cty (year): replace dbase = dbase[1]
gen double dose = dbase + rnormal(0, 0.5)
gen double aunit = 0
bysort cty (year): replace aunit = 3*runiform() if _n == 1
bysort cty (year): replace aunit = aunit[1]
gen double y = 5 + aunit + 0.25*year + 0.7*dose + 0.4*dose*post + 0.5*post + rnormal(0,1)
tempfile F
qui save `F'

dis as result "=============== S1: what the user SEES on the default path ==============="
cap noi wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust
dis "S1_RC=" _rc

dis as result "=============== S2: contents of e(V) ==============="
cap {
  mat V = e(V)
  local nr = rowsof(V)
  dis "S2 e(V) dim = `nr'"
  mata: V = st_matrix("V")
  mata: st_local("nmiss", strofreal(sum(V:>=.)))
  mata: st_local("mx", strofreal(max(abs(V:*(V:<.)))))
  dis "S2 e(V) n_missing_cells = `nmiss'   max|finite cell| = `mx'"
  mat b = e(b)
  dis "S2 e(b) colnames:"
  mat li b
  dis "S2 e(V) diag:"
  mat d = vecdiag(V)
  mat li d
}
dis "S2_RC=" _rc

dis as result "=============== S3: downstream test / lincom on posted e(b),e(V) ==============="
cap noi test [#1]_b[#1] = 0
dis "S3a_RC=" _rc
cap noi lincom [#1]_b[#1]
dis "S3b_RC=" _rc
qui cap {
  mat b = e(b)
  local c1 : word 1 of `: colnames b'
}
cap noi test `c1' = 0
dis "S3c_RC=" _rc

dis as result "=============== S4: boottstat on the zero-V spec ==============="
use `F', clear
cap noi wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust bootstrap(standard) reps(30) boottstat rseed(99)
dis "S4_RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "S4 boottstat att  b=" M[1,1] " se=" M[1,2] " p=" M[1,4]
  dis "S4 boottstat ame  b=" M[2,1] " se=" M[2,2] " p=" M[2,4]
  mat li M
}

dis as result "=============== S5: incidence across realistic specs ==============="
foreach spec in "lattice(all 3) robust" ///
                "lattice(all 3) cluster(cty)" ///
                "lattice(all 2) robust" ///
                "lattice(all 2) cluster(cty)" ///
                "lattice(reg3) robust" ///
                "lattice(reg3) cluster(cty)" ///
                "lattice(all 5) robust" ///
                "robust" ///
                "cluster(cty)" {
  use `F', clear
  cap qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) `spec'
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "S5 [`spec']  att b=" %12.7f M[1,1] " se=" %12.7g M[1,2] "  ame b=" %12.6f M[2,1] " se=" %12.7g M[2,2] "  N=" e(N) " df_r=" e(df_r)
  }
  else dis as error "S5 [`spec'] RC=" _rc
}
dis "DONE_L2"

dis as result "=============== S6: is i-above-the-unit FORCED by contreat()? ==============="
use `F', clear
cap noi qui wooldid y cty year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust
dis "S6a i=cty lattice RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "S6a att b=" M[1,1] " se=" M[1,2] "  ame b=" M[2,1] " se=" M[2,2]
}
use `F', clear
cap noi qui wooldid y cty year ttre , att contreat(dose) contreatcontrols(basic) robust
dis "S6b i=cty nolattice RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "S6b att b=" M[1,1] " se=" M[1,2] "  ame b=" M[2,1] " se=" M[2,2]
}
dis "DONE_L2B"
