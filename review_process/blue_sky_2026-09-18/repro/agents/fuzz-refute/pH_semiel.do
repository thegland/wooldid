*=====================================================================
* pH: FUZZ-02 -- semielasticity when the fitted values straddle zero.
* Is the r(459) wooldid's, or plain margins' pass-through?
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* REALISTIC: outcome is an annual growth rate, which straddles zero
clear
set seed 606
set sortseed 606
qui set obs 60
gen long unit = _n
gen double a_i = rnormal(0,0.5)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.4 & u <= 0.6
qui replace ttre = 2005 if u > 0.6 & u <= 0.8
qui replace ttre = 2007 if u > 0.8
qui expand 10
bysort unit: gen int t = 1999 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double g = a_i + 0.3*rnormal() + 0.6*D + 0.2*x1
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i
qui sum g
dis "PH_G_MIN=" r(min) "  PH_G_MAX=" r(max)
tempfile gd
qui save "`gd'", replace

cap noi qui wooldid g i t ttre , att cluster(i)
dis "PH_CTRL_RC=" _rc
cap noi wooldid g i t ttre , att cluster(i) semielasticity
dis "PH_SEMI_i_RC=" _rc
cap noi wooldid g icoh t ttre , att cluster(icoh) semielasticity
dis "PH_SEMI_coh_RC=" _rc
cap noi wooldid g icoh t ttre , att cluster(icoh) semielasticity nofastcollapse
dis "PH_SEMI_coh_nofc_RC=" _rc
cap noi wooldid g i t ttre , att cluster(i) semielasticity safety(off)
dis "PH_SEMI_i_safoff_RC=" _rc
cap noi wooldid g i t ttre , att robust semielasticity safety(off)
dis "PH_SEMI_i_rob_safoff_RC=" _rc

* --- is plain margins the same? fit a comparable model by hand ---
use "`gd'", clear
qui reghdfe g x1 i.icoh#i.t#i.D, absorb(i t) vce(cluster i)
cap noi margins, dydx(x1)
dis "PH_MARG_DYDX_RC=" _rc
cap noi margins, eydx(x1)
dis "PH_MARG_EYDX_RC=" _rc
qui predict double xbh, xbd
qui sum xbh
dis "PH_XB_MIN=" r(min) " PH_XB_MAX=" r(max)

* --- decisive control: strictly positive outcome, same design ---
use "`gd'", clear
qui replace g = g + 50
cap noi qui wooldid g i t ttre , att cluster(i) semielasticity
dis "PH_POS_SEMI_RC=" _rc
qui reghdfe g x1 i.icoh#i.t#i.D, absorb(i t) vce(cluster i)
cap noi margins, eydx(x1)
dis "PH_POS_MARG_EYDX_RC=" _rc
dis "=================== DONE pH ==================="
