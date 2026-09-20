clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- self-contained staggered-DiD panel (60 units, T=10, 3 onsets) ----
clear
set seed 31337
set sortseed 31337
qui set obs 60
gen long unit = _n
gen double a_i = rnormal(0,1)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.35 & u <= 0.5667
qui replace ttre = 2005 if u > 0.5667 & u <= 0.7833
qui replace ttre = 2007 if u > 0.7833
qui expand 10
bysort unit: gen int t = 1999 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double y = a_i + 0.4*(t-2000) + 2*D + 0.2*x1 + 0.3*rnormal()
gen double dose = D*(0.5 + runiform())
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i

* ==== NEW: level() restore corrupts c(level) and reports success as failure ====
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
* control: an ordinary session level is restored correctly
set level 90
cap noi qui wooldid y i t ttre , att robust safety(off) level(99)
dis "CONTROL rc=" _rc " e(level)=" e(level) " c(level)_after=" c(level) " (should be 90)"
* the defect: session level 99.9, a perfectly legal level(95) request
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 99.9
dis "BEFORE c(level)=" c(level)
cap noi wooldid y i t ttre , att robust safety(off) level(95)
dis "AFTER rc=" _rc
dis "   e(cmd)=" e(cmd) "  e(N)=" e(N) "  e(level)=" e(level) "   <- the estimation SUCCEEDED and posted"
dis "   c(level)=" c(level) "   <- the user's 99.9 is GONE, silently replaced"
dis "   WOOLDID_LEVEL0 leaked = [${WOOLDID_LEVEL0}]"
set level 95
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
