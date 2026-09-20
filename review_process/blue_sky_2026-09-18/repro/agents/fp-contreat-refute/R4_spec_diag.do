* R4: what exactly does wooldid fit under contreat()+basic?  verbose echo, and
* three candidate outside-wooldid spellings compared on dydx.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

clear
set obs 720
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/12)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
qui replace ttre = 8 if i == 3
gen byte post = (ttre < . & t >= ttre)
gen double dose = 1 + 2*runiform()
gen double bit = 0
qui replace bit = 0.60 + 0.30*i - 0.09*(t-ttre) if post == 1
gen double kit = 0
qui replace kit = 0.40 + 0.15*i + 0.05*(t-ttre) if post == 1
bysort unit (t): gen double aunit = 20 + 6*runiform() if _n == 1
bysort unit (t): replace aunit = aunit[1]
gen double y = aunit + 0.35*t + 0.9*dose + kit + bit*dose + rnormal(0, 0.5)
gen double ringer = dose
gen byte nullv = 1
egen cell = group(i t) if post == 1
qui replace cell = 0 if cell >= .
tempfile F
qui save `F'

dis "=== VERBOSE wooldid, default path ==="
use `F', clear
cap noi wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) robust regtol(14) verbose
dis "RCV=" _rc

dis "=== SPELLING 1: c.ringer + ib0.cell + ib0.cell#c.dose, absorb(i.i i.t) ==="
use `F', clear
reghdfe y c.ringer ib0.cell ib0.cell#c.dose , absorb(i.i i.t) vce(robust)
qui margins , dydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "SP1 dydx = " %18.12f r[1,1]

dis "=== SPELLING 2: ib0.cell##c.dose, absorb(i.i i.t i.nullv##c.ringer) ==="
use `F', clear
reghdfe y ib0.cell##c.dose , absorb(i.i i.t i.nullv##c.ringer) vce(robust)
qui margins , dydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "SP2 dydx = " %18.12f r[1,1]
qui margins , eydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "SP2 eydx = " %18.12f r[1,1]
qui predict double xb2 , xb
qui sum xb2
dis as result "SP2 xb mean " %12.6f r(mean) " sd " %12.6f r(sd)

dis "=== SPELLING 3: reg y i.i i.t i.nullv##c.ringer ib0.cell##c.dose ==="
use `F', clear
reg y i.i i.t i.nullv##c.ringer ib0.cell##c.dose , vce(robust)
qui margins , dydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "SP3 dydx = " %18.12f r[1,1]
qui margins , eydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "SP3 eydx = " %18.12f r[1,1]
qui predict double xb3 , xb
qui sum xb3
dis as result "SP3 xb mean " %12.6f r(mean) " sd " %12.6f r(sd)
qui sum y
dis as result "  y  mean " %12.6f r(mean) " sd " %12.6f r(sd)

dis "DONE_R4"
