* p5_conseq.do -- SELF-CONTAINED.  Consequence surface of the two guards that
* contreatwithin + contreat() silently disables.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case"
which wooldid
dis "BUILDCHECK_DONE"

*---- ONE panel with an 8-digit float id ----
clear
set obs 60
gen long unit = _n
gen float ibig = 20000000 + unit
gen double ttre = .
replace ttre = 5 if ibig <= 20000020
replace ttre = 8 if ibig > 20000020 & ibig <= 20000040
gen double ufe = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort unit: gen double t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
gen double ibig_d = 20000000 + unit
save "`WD'/panelC.dta", replace
qui tab ibig
dis "DISTINCT_FLOAT=" r(r)
qui tab ibig_d
dis "DISTINCT_DOUBLE=" r(r)

dis "=== C1 float id + cluster(i), contreatwithin (no caution) ==="
use "`WD'/panelC.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin cluster(i)
dis "C1_RC=" _rc
cap noi dis "C1_ATT=" %21x el(e(b),1,1)
cap noi dis "C1_SE=" %21x sqrt(el(e(V),1,1))
cap noi dis "C1_SEd=" sqrt(el(e(V),1,1))
cap noi dis "C1_NCLUST=" e(nclusters)

dis "=== C2 double id + cluster(i), contreatwithin (the panel the user meant) ==="
use "`WD'/panelC.dta", clear
gen double i = ibig_d
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin cluster(i)
dis "C2_RC=" _rc
cap noi dis "C2_ATT=" %21x el(e(b),1,1)
cap noi dis "C2_SE=" %21x sqrt(el(e(V),1,1))
cap noi dis "C2_SEd=" sqrt(el(e(V),1,1))
cap noi dis "C2_NCLUST=" e(nclusters)

dis "=== C3 float id + cluster(i), NO contreatwithin (caution prints) ==="
use "`WD'/panelC.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i)
dis "C3_RC=" _rc

*---- negative / non-integer identifiers: informative refusal vs opaque crash ----
dis "=== D1 non-integer i, contreatwithin ==="
use "`WD'/panelC.dta", clear
gen double i = unit
qui replace i = i + 0.5 if unit == 3
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "D1_RC=" _rc

dis "=== D2 non-integer i, NO contreatwithin ==="
use "`WD'/panelC.dta", clear
gen double i = unit
qui replace i = i + 0.5 if unit == 3
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "D2_RC=" _rc

dis "=== D3 negative i, contreatwithin ==="
use "`WD'/panelC.dta", clear
gen double i = unit
qui replace i = -9 if unit == 3
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "D3_RC=" _rc

dis "=== D4 negative i, NO contreatwithin ==="
use "`WD'/panelC.dta", clear
gen double i = unit
qui replace i = -9 if unit == 3
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "D4_RC=" _rc

dis "ALLDONE_P5"
