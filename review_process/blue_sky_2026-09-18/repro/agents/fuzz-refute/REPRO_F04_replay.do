clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- self-contained staggered-DiD panel (60 units, T=10, 3 onsets) ----
clear
set seed 808
set sortseed 808
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

* ============ FUZZ-04: estimates replay, and the LASTRC poisoning ============
macro drop WOOLDID_LASTRC
cap noi qui wooldid y i t ttre , att robust safety(off)
dis "EST rc=" _rc
cap noi qui wooldid
dis "BARE BEFORE rc=" _rc "   (bare replay works)"
* REFUTATION of the trailing-blank sub-claim: two trailing blanks still replay
cap noi qui wooldid  
dis "BARE + TWO TRAILING BLANKS rc=" _rc "   (claim said this defeats the gate)"
cap noi estimates store m1
cap noi qui wooldid
dis "BARE AFTER STORE rc=" _rc
cap noi estimates replay m1
dis "ESTIMATES REPLAY rc=" _rc
dis "   WOOLDID_LASTRC is now [${WOOLDID_LASTRC}]"
cap noi wooldid
dis "BARE AFTER FAILED REPLAY rc=" _rc "   <- the WORKING route is now refused"
cap noi estimates restore m1
cap noi wooldid
dis "BARE AFTER RESTORE rc=" _rc "   <- still refused; only a refit clears it"
cap noi qui wooldid y i t ttre , att robust safety(off)
cap noi qui wooldid
dis "BARE AFTER REFIT rc=" _rc
* controls: other post-estimation commands are unaffected
cap noi qui estimates table m1
dis "ESTIMATES TABLE rc=" _rc
