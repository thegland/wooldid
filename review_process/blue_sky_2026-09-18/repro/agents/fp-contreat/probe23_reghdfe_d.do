* Mechanism + remedy for FINDING A: reghdfe is called WITHOUT the `d` option
* (ado:1476), so margins' prediction after it is xb WITHOUT the absorbed fixed
* effects. ppmlhdfe on the poisson path IS called with `d` (ado:1462/1468).
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 60601
set sortseed 60601

clear
set obs 1000
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/5)
gen double ttre = .
qui replace ttre = 4 if inrange(i, 1, 5)
qui replace ttre = 6 if inrange(i, 6,10)
qui replace ttre = 8 if inrange(i,11,15)
gen double sfe = 40 + 3*i
gen double dose = 1 + 2*runiform()
gen double kappatrue = 0
gen double betatrue  = 0
qui replace kappatrue = 1.0 + 0.20*(t-ttre) if ttre < . & t >= ttre
qui replace betatrue  = 2.0 - 0.30*(t-ttre) if ttre < . & t >= ttre
gen double y = sfe + 0.5*t + 0.9*dose + kappatrue + betatrue*dose
gen double relyr = t - ttre
gen byte tf = (ttre<. & relyr>=0)
egen long cell = group(i relyr) if tf
qui replace cell = 0 if cell==.
gen double ringer = dose
gen byte nullv = 1

gen double q = betatrue/y
qui sum q if tf
dis "TRUTH eydx = " %18.12f r(mean)

dis "===== reghdfe WITHOUT d ====="
qui reghdfe y ib0.cell##c.dose , absorb(i.i i.t i.nullv##c.ringer) vce(robust) tol(1e-14)
cap drop xbh
qui predict double xbh , xb
qui gen double dev = y - xbh
qui sum dev
dis "NO_D: mean(y-xb)=" %12.6f r(mean) " sd=" %12.6f r(sd)
qui margins , eydx(dose) subpop(if tf) noestimcheck nose
mat r1 = r(b)
dis "NO_D eydx margins = " %18.12f r1[1,1]

dis "===== reghdfe WITH d ====="
qui reghdfe y ib0.cell##c.dose , absorb(i.i i.t i.nullv##c.ringer) vce(robust) tol(1e-14) d
cap drop xbh2
qui predict double xbh2
qui gen double dev2 = y - xbh2
qui sum dev2
dis "WITH_D: mean(y-pred)=" %12.6f r(mean) " sd=" %12.6f r(sd)
cap noisily qui margins , eydx(dose) subpop(if tf) noestimcheck nose
dis "WITH_D_RC=" _rc
if _rc==0 {
  mat r2 = r(b)
  dis "WITH_D eydx margins = " %18.12f r2[1,1]
}

* ===== poisson remap check (BS-6): under poisson without poisexpresults the help
* says eydx is COMPUTED AS dydx and eyex AS dyex - "the option changes the table
* header, not the estimate". Verify the four runs collapse to two numbers. =====
clear
set obs 384
gen long unit = ceil(_n/8)
gen int  t    = _n - (unit-1)*8
gen byte i    = ceil(unit/12)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
gen double dose = 1 + 2*runiform()
gen double xbi = 1 + 0.15*i + 0.05*t + 0.10*dose
qui replace xbi = xbi + 0.10 + 0.20*dose if ttre<. & t>=ttre
gen double ypo = exp(xbi + rnormal(0,0.2))
foreach et in dydx dyex eyex eydx {
  cap qui wooldid ypo i t ttre , att contreat(dose) contreatcontrols(basic) poisson ///
      contreatelasticitytype(`et') robust regtol(9)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "PREMAP|`et'|ame=" %21x M[2,1]
  }
  else dis "PREMAP|`et'|RC=" _rc
}
dis "DONE23"
