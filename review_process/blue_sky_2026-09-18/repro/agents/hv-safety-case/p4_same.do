* p4_same.do -- SELF-CONTAINED.  Same-data version of the contreatwithin probe.
* The panel is generated ONCE and saved; every leg reloads it, so the only thing
* that differs between legs is the identifier or the option.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case"
which wooldid
dis "BUILDCHECK_DONE"

*---- ONE panel, generated once ----
clear
set obs 60
gen long unit = _n
gen double ttre = .
replace ttre = 5 if unit <= 15
replace ttre = 8 if unit > 15 & unit <= 30
gen double ufe = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort unit: gen double t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
save "`WD'/panel.dta", replace

*===================================================================
* B: i == 0 aliases the never-treated code inside wooldid (_i)
*===================================================================
dis "=== B0 ids 1..60, NO contreatwithin (reference run, guard passes) ==="
use "`WD'/panel.dta", clear
gen double i = unit
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "B0_RC=" _rc
cap noi dis "B0_ATT_HEX=" %21x el(e(b),1,1)
cap noi dis "B0_N=" e(N)

dis "=== B1 ids 0..59, NO contreatwithin (expect refusal 411) ==="
use "`WD'/panel.dta", clear
gen double i = unit - 1
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "B1_RC=" _rc

dis "=== B2 ids 0..59, WITH contreatwithin ==="
use "`WD'/panel.dta", clear
gen double i = unit - 1
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B2_RC=" _rc
cap noi dis "B2_ATT_HEX=" %21x el(e(b),1,1)
cap noi dis "B2_N=" e(N)

dis "=== B3 ids 1..60 (SAME DATA, ids shifted by one), WITH contreatwithin ==="
use "`WD'/panel.dta", clear
gen double i = unit
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B3_RC=" _rc
cap noi dis "B3_ATT_HEX=" %21x el(e(b),1,1)
cap noi dis "B3_N=" e(N)

dis "=== B4 ids 101..160 (SAME DATA, ids shifted by 100), WITH contreatwithin ==="
use "`WD'/panel.dta", clear
gen double i = unit + 100
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B4_RC=" _rc
cap noi dis "B4_ATT_HEX=" %21x el(e(b),1,1)
cap noi dis "B4_N=" e(N)

*===================================================================
* A: float 8-digit id, onset assigned from the ROUNDED id so the
*    onset-constancy screen does not fire.  Same data both legs.
*===================================================================
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
qui tab ibig
dis "A_DISTINCT_FLOAT=" r(r)
qui tab ibig_d
dis "A_DISTINCT_DOUBLE=" r(r)
save "`WD'/panelA.dta", replace

dis "=== A1 float id, contreat, NO contreatwithin (caution expected) ==="
use "`WD'/panelA.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "A1_RC=" _rc
cap noi dis "A1_ATT_HEX=" %21x el(e(b),1,1)

dis "=== A2 float id, contreat + contreatwithin (caution?) ==="
use "`WD'/panelA.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "A2_RC=" _rc
cap noi dis "A2_ATT_HEX=" %21x el(e(b),1,1)

dis "=== A3 double id (what the caution asks for), contreat + contreatwithin ==="
use "`WD'/panelA.dta", clear
gen double i = ibig_d
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "A3_RC=" _rc
cap noi dis "A3_ATT_HEX=" %21x el(e(b),1,1)

dis "=== A4 float id + contreatwithin + safety(on) explicitly ==="
use "`WD'/panelA.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin safety(on)
dis "A4_RC=" _rc

dis "ALLDONE_P4"
