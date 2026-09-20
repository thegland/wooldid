* r3_ipos.do -- SELF-CONTAINED. Skeptic repro of hv-ctw-ipos.
* Does contreatwithin + contreat() skip the "i must consist of positive,
* non-zero integers" refusal at ado:4959?
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
replace ttre = 4 if _n <= 15
replace ttre = 6 if _n > 15 & _n <= 30
gen double ue = rnormal()
expand 8
bys id: gen double t = _n
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
tempfile P
save "`P'"

* three bad shapes for ONE treated unit's identifier
foreach bad in 0.5 -3 0 {
  local tag = subinstr("`bad'",".","p",.)
  local tag = subinstr("`tag'","-","m",.)

  * --- WITHOUT contreatwithin: the documented refusal ---
  use "`P'", clear
  qui replace id = `bad' if id == 1
  dis "LEGSTART_noctw_`tag'"
  cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort)
  dis "LEG_noctw_`tag'_RC=" _rc

  * --- WITH contreatwithin: the claim says the refusal never runs ---
  use "`P'", clear
  qui replace id = `bad' if id == 1
  dis "LEGSTART_ctw_`tag'"
  cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin
  dis "LEG_ctw_`tag'_RC=" _rc

  * --- binary design, explicit safety(off): the DOCUMENTED bypass, for contrast ---
  use "`P'", clear
  qui replace id = `bad' if id == 1
  dis "LEGSTART_sfoff_`tag'"
  cap noi wooldid y id t ttre, att safety(off)
  dis "LEG_sfoff_`tag'_RC=" _rc
}

dis "ALLDONE"
