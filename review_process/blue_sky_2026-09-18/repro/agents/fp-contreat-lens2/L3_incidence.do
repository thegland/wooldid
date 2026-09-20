* LENS 2 part 3: incidence of the zero-V across the specs a real applied
* economist would actually write (S5), and whether i-above-the-unit is FORCED
* by contreat() (S6). The boottstat leg is in L4 -- it terminated the L2 run.
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
    dis as result "S5 [`spec']  attb=" %11.6f M[1,1] " attse=" %11.6g M[1,2] " ameb=" %11.5f M[2,1] " amese=" %11.6g M[2,2] " N=" e(N) " dfr=" e(df_r)
  }
  else dis as error "S5 [`spec'] RC=" _rc
}

dis as result "=============== S6: is i-above-the-unit forced by contreat()? ==============="
use `F', clear
cap noi qui wooldid y cty year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust
dis "S6a i=cty lattice RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "S6a attb=" M[1,1] " attse=" M[1,2] " ameb=" M[2,1] " amese=" M[2,2]
}
use `F', clear
cap noi qui wooldid y cty year ttre , att contreat(dose) contreatcontrols(basic) robust
dis "S6b i=cty nolattice RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "S6b attb=" M[1,1] " attse=" M[1,2] " ameb=" M[2,1] " amese=" M[2,2]
}
dis "DONE_L3"
