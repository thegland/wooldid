* r4_izero.do -- SELF-CONTAINED. Does a TREATED cohort whose identifier is 0
* (which contreatwithin+contreat, and safety(off), both let through) change the
* reported numbers relative to the identical panel with a nonzero label?
* wooldid sets _i = i on estimable treated rows and 0 elsewhere (ado:928-929), so
* i == 0 aliases a treated cohort onto the never-treated code. ado:1602 builds the
* att_i / att_itime weight classifier with an explicit `_i != 0' condition.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

clear
set obs 40
gen double unit = _n
gen double ttre = .
replace ttre = 5 if _n <= 15
replace ttre = 6 if _n > 15 & _n <= 30
gen double ue = rnormal()
expand 8
bys unit: gen double t = _n
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
tempfile P
save "`P'"

* A: the treated unit 1 carries identifier 0
* B: the SAME panel with every identifier shifted +1 (no zero; identical partition)
foreach leg in A B {
  use "`P'", clear
  if "`leg'" == "A" gen double id = cond(unit == 1, 0, unit)
  if "`leg'" == "B" gen double id = unit + 1

  dis "LEGSTART_`leg'_ctw"
  cap noi wooldid y id t ttre, att att_i att_itime att_it esprelength(3) espostlength(2) contreat(dose) contreatcontrols(byCohort) contreatwithin
  dis "LEG_`leg'_ctw_RC=" _rc
  dis "MATSTART_`leg'_ctw_b"
  cap noi mat li e(b), format(%21x)
  dis "MATEND_`leg'_ctw_b"
  dis "MATSTART_`leg'_ctw_main"
  cap noi mat li e(wdidmainresults), format(%21x)
  dis "MATEND_`leg'_ctw_main"

  * binary twin under the DOCUMENTED bypass safety(off)
  use "`P'", clear
  if "`leg'" == "A" gen double id = cond(unit == 1, 0, unit)
  if "`leg'" == "B" gen double id = unit + 1
  dis "LEGSTART_`leg'_bin"
  cap noi wooldid y id t ttre, att att_i att_itime att_it esprelength(3) espostlength(2) safety(off)
  dis "LEG_`leg'_bin_RC=" _rc
  dis "MATSTART_`leg'_bin_b"
  cap noi mat li e(b), format(%21x)
  dis "MATEND_`leg'_bin_b"
  dis "MATSTART_`leg'_bin_main"
  cap noi mat li e(wdidmainresults), format(%21x)
  dis "MATEND_`leg'_bin_main"
}

dis "ALLDONE"
