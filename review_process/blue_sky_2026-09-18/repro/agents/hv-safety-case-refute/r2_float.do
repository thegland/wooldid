* r2_float.do -- SELF-CONTAINED. Skeptic repro of hv-ctw-float.
* Does contreatwithin + contreat() suppress the BS-34 float-above-2^24 caution?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* ---- an administrative panel: 8-digit provider codes held in FLOAT (gen's default) ----
clear
set obs 60
gen float id = 20000000 + _n          // float cannot hold consecutive ints here
egen double g = group(id)             // g counts the DISTINCT (rounded) codes
qui sum g
dis "DISTINCT_FLOAT_IDS=" r(max) " out of 60"
gen double idd = 20000000 + _n        // the double twin of the same codes
* onset assigned FROM the rounded code so onset stays constant within each merged id
gen double ttre = .
replace ttre = 4 if mod(g,3) == 1
replace ttre = 6 if mod(g,3) == 2
gen double ue = rnormal()
expand 8
bys id (ttre): gen double t = _n
replace t = mod(_n-1,8)+1
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
tempfile P
save "`P'"

prog def shownum
  args tag
  tempname b
  mat `b' = e(b)
  dis "`tag'_N=" e(N)
  dis "`tag'_ATT=" %21x `b'[1,1]
  dis "`tag'_SE=" %21x sqrt(e(V)[1,1])
end

* LEG1: contreat, NO contreatwithin -> caution expected
use "`P'", clear
dis "LEGSTART_L1"
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort) cluster(id)
dis "LEG_L1_RC=" _rc
cap shownum L1

* LEG2: + contreatwithin -> the claim says the caution vanishes
use "`P'", clear
dis "LEGSTART_L2"
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin cluster(id)
dis "LEG_L2_RC=" _rc
cap shownum L2

* LEG3: + contreatwithin + the user EXPLICITLY asking for safety(on)
use "`P'", clear
dis "LEGSTART_L3"
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin cluster(id) safety(on)
dis "LEG_L3_RC=" _rc
cap shownum L3

* LEG4: the panel the user MEANT (double ids), same options as LEG2
use "`P'", clear
dis "LEGSTART_L4"
cap noi wooldid y idd t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin cluster(idd)
dis "LEG_L4_RC=" _rc
cap shownum L4

* LEG5: binary design, contreatwithin absent, float id -> caution control
use "`P'", clear
dis "LEGSTART_L5"
cap noi wooldid y id t ttre, att cluster(id)
dis "LEG_L5_RC=" _rc
cap shownum L5

* LEG6: is ANYTHING returned that records the internal relaxation?
use "`P'", clear
dis "LEGSTART_L6"
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin cluster(id)
dis "LEG_L6_RC=" _rc
dis "L6_safety_e=" e(safety)
dis "L6_wdid_safety_e=" e(wdid_safety)
dis "L6_cmdline=" e(cmdline)
qui ereturn list
dis "ALLDONE"
