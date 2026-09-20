clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31415
set sortseed 31415

* ================= (a) contreatpoly(2) =====================================
cap prog drop mkpoly
prog def mkpoly
  clear
  set obs 320                       /* 4 groups x 10 units x 8 periods */
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/10)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  gen double b1 = 0
  gen double b2 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.05*(t-ttre) if ttre<. & t>=ttre
  qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre<. & t>=ttre
  qui replace b2 = -0.15 + 0.05*i + 0.02*(t-ttre) if ttre<. & t>=ttre
  gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose - 0.2*dose^2 + k0 + b1*dose + b2*dose^2
end
foreach path in "robust" "customvce(robust)" {
  qui mkpoly
  gen double dm = b1 + 2*b2*dose
  qui sum dm if ttre<. & t>=ttre
  local tru = r(mean)
  qui sum b1 if ttre<. & t>=ttre
  cap qui wooldid y i t ttre , att contreat(dose) contreatpoly(2) contreatcontrols(basic) `path' regtol(14)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "POLY|`path'|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|r2=" %14.12f e(r2)
  }
  else dis "POLY|`path'|RC=" _rc
}

* ================= (b) user lattice, bin-specific slopes ===================
cap prog drop mklat
prog def mklat
  clear
  set obs 320
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/10)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  gen byte lat = 1 + (dose > 2)
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.4*lat        if ttre<. & t>=ttre
  qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) + 0.6*(lat==2) if ttre<. & t>=ttre
  * lattice-specific baseline dose slope, which contreatcontrols(basic)+lattice must absorb
  gen double y = 2 + 0.6*i + 0.35*t + (0.5 + 0.7*(lat==2))*dose + k0 + b1*dose
end
foreach path in "robust" "customvce(robust)" {
  qui mklat
  qui sum b1 if ttre<. & t>=ttre
  local tru = r(mean)
  cap qui wooldid y i t ttre , att contreat(dose) lattice(lat) contreatcontrols(basic) `path' regtol(14)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "LATU|`path'|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|r2=" %14.12f e(r2) "|se=" %12.6g M[2,2]
  }
  else dis "LATU|`path'|RC=" _rc
}

* ================= (c) untreatedZero =======================================
cap prog drop mkuz
prog def mkuz
  clear
  set obs 320
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/10)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 0
  qui replace dose = 1 + 2*runiform() if ttre<. & t>=ttre     /* zero when untreated */
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre<. & t>=ttre
  qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre<. & t>=ttre
  gen double y = 2 + 0.6*i + 0.35*t + k0 + b1*dose
end
foreach path in "robust" "customvce(robust)" {
  qui mkuz
  qui sum b1 if ttre<. & t>=ttre
  local tru = r(mean)
  cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) `path' regtol(14)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "UZ|`path'|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|r2=" %14.12f e(r2)
  }
  else dis "UZ|`path'|RC=" _rc
}
dis "DONE"
