* r7_bare198.do -- SELF-CONTAINED.
* A panel whose treated units have no usable pre-period. Under default safety the
* run is refused with a named r(459). contreatwithin + contreat() rewrites the
* working safety to off (ado:4829) above that screen (ado:908), and the run then
* dies with a BARE r(198) carrying no message at all.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

clear
set obs 40
gen double id = _n
gen double ttre = .
replace ttre = 1 if _n <= 30
gen double ue = rnormal()
expand 8
bys id: gen double t = _n
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
tempfile Q
save "`Q'"

use "`Q'", clear
dis "LEGSTART_1_default"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort)
dis "LEG_1_default_RC=" _rc

use "`Q'", clear
dis "LEGSTART_2_ctw"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort) contreatwithin
dis "LEG_2_ctw_RC=" _rc

use "`Q'", clear
dis "LEGSTART_3_safetyoff_noctw"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort) safety(off)
dis "LEG_3_safetyoff_noctw_RC=" _rc

use "`Q'", clear
dis "LEGSTART_4_binary_safetyoff"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) safety(off)
dis "LEG_4_binary_safetyoff_RC=" _rc

use "`Q'", clear
dis "LEGSTART_5_ctw_verbose"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort) contreatwithin verbose
dis "LEG_5_ctw_verbose_RC=" _rc

dis "ALLDONE"
