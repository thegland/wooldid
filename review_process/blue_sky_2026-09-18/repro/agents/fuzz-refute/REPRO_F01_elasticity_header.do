clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- self-contained staggered-DiD panel (60 units, T=10, 3 onsets) ----
clear
set seed 4242
set sortseed 4242
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

* ============ FUZZ-01: contreatelasticitytype() in BINARY mode ============
* sthlp:908 promises: "without contreat(varname) the option is otherwise inert,
* its value read by nothing and no message printed."
cap noi wooldid y i t ttre , att robust safety(off)
dis "BASE_RC=" _rc
scalar b0 = e(b)[1,1]
matrix V0 = e(V)
scalar s0 = sqrt(V0[1,1])
dis "BASE_ATT=" %21x b0 " SE=" %21x s0
foreach ct in dydx dyex eydx eyex {
  cap noi wooldid y i t ttre , att robust safety(off) contreatelasticitytype(`ct')
  scalar bx = e(b)[1,1]
  matrix Vx = e(V)
  scalar sx = sqrt(Vx[1,1])
  dis "F01 `ct' rc=" _rc "  bitidentical_b=" (bx==b0) "  bitidentical_se=" (sx==s0) "  e(wdid_cetype)=[" e(wdid_cetype) "]"
}
* the stored header is replayed too:
cap noi wooldid y i t ttre , att robust safety(off) contreatelasticitytype(eyex)
dis "---- bare replay below repeats the same header ----"
cap noi wooldid
* what the OTHER inert companions do in the same position (claim: all refused - false)
foreach o in "lattice(all 3)" "contreatcontrols(basic)" "contreatwithin" "makecfxplots" "cfxlattice(dose)" "latticeignoreweights" "ccnoabsorb" "contreatpoly(2)" "cfxplottypes(0)" "makecfxplotsbysg" "cfxplotnose" {
  cap noi qui wooldid y i t ttre , att robust safety(off) `o'
  dis "COMPANION [`o'] rc=" _rc
}
