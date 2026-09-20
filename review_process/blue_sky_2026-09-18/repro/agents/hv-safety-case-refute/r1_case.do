* r1_case.do -- SELF-CONTAINED. Skeptic repro of hv-safety-case-1.
* Does an unrecognised CASING of safety() ever get silently accepted, and
* which message does the user see first?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* ---- my own panel: 40 units, 10 periods, staggered onset, 10 never-treated ----
clear
set obs 40
gen unit = _n
gen ttre = .
replace ttre = 4 if unit <= 10
replace ttre = 6 if unit > 10 & unit <= 20
replace ttre = 8 if unit > 20 & unit <= 30
gen ue = rnormal()
expand 10
bys unit: gen t = _n
gen dose = runiform() + 0.5
gen post = (t >= ttre) & ttre < .
gen y = 1 + 0.3*t + ue + 0.8*post + rnormal()*0.5
gen wt = 1 + runiform()
tempfile P
save "`P'"

foreach s in "off" "OFF" "Off" "oFF" "on" "ON" "On" {
  use "`P'", clear
  dis "LEGSTART_plain_`s'"
  cap noi wooldid y unit t ttre, att safety(`s')
  dis "LEG_plain_`s'_RC=" _rc
  dis "LEGEND_plain_`s'"
}

* empty parens and quoted token
use "`P'", clear
dis "LEGSTART_empty"
cap noi wooldid y unit t ttre, att safety()
dis "LEG_empty_RC=" _rc
use "`P'", clear
dis "LEGSTART_quoted"
cap noi wooldid y unit t ttre, att safety("off")
dis "LEG_quoted_RC=" _rc

* ---- BITE 1: a safety-keyed REFUSAL fires before the domain screen ----
use "`P'", clear
dis "LEGSTART_bite1"
cap noi wooldid y unit t ttre, att_i unconditionalse safety(OFF)
dis "LEG_bite1_RC=" _rc

* control: same line with correct casing runs
use "`P'", clear
dis "LEGSTART_bite1ok"
cap noi wooldid y unit t ttre, att_i unconditionalse safety(off)
dis "LEG_bite1ok_RC=" _rc

* ---- BITE 2: a safety-keyed CAUTION prints, then the run dies on the spelling ----
use "`P'", clear
dis "LEGSTART_bite2"
cap noi wooldid y unit t ttre [aw=wt], att unconditionalse robust safety(OFF)
dis "LEG_bite2_RC=" _rc

* ---- the one `== "off"' site inside the window (ado:4831): is it conservative? ----
use "`P'", clear
dis "LEGSTART_cfx"
cap noi wooldid y unit t ttre, att_i contreat(dose) contreatcontrols(byCohort) makecfxplots safety(OFF)
dis "LEG_cfx_RC=" _rc
use "`P'", clear
dis "LEGSTART_cfxoff"
cap noi wooldid y unit t ttre, att_i contreat(dose) contreatcontrols(byCohort) makecfxplots safety(off)
dis "LEG_cfxoff_RC=" _rc

dis "ALLDONE"
