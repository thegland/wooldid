clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- self-contained staggered-DiD panel (60 units, T=10, 3 onsets) ----
clear
set seed 5150
set sortseed 5150
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

* ============ FUZZ-06: dots(-5) validated on one engine only ============
cap noi qui wooldid y i t ttre , att cluster(i) permute(shuffle) permreps(24) rseed(9) dots(-5)
dis "PERMUTE dots(-5) rc=" _rc "   (0 = silently accepted)"
cap noi qui wooldid y i t ttre , att cluster(i) permute(shuffle) permreps(24) rseed(9) dots(0)
dis "PERMUTE dots(0)  rc=" _rc
cap noi wooldid y i t ttre , att cluster(i) bootstrap(standard) reps(24) rseed(9) dots(-5)
dis "BOOTSTRAP dots(-5) rc=" _rc "   (198 = refused)"
