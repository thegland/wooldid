* r5_izero_strong.do -- SELF-CONTAINED.
* A treated cohort whose identifier is literally 0 is silently DOWNWEIGHTED in every
* pre-treatment att_i / att_itime row. Reachable two ways: the documented safety(off)
* bypass, and -- with no safety() option typed at all -- contreatwithin + contreat(),
* which rewrites the working safety to off above the i-positivity screen (ado:4829/4959).
*
* Leg A and leg B are THE SAME PANEL. B's identifiers are A's plus one. Nothing but the
* labels differs, so every reported number must be identical.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* 20 cohorts x 3 units x 8 periods. cohort 1 carries a distinctive pre-trend.
clear
set obs 60
gen double cohort = ceil(_n/3)
gen double ttre = .
replace ttre = 5 if cohort <= 7
replace ttre = 6 if cohort > 7 & cohort <= 14
gen double ue = rnormal()
expand 8
bys cohort ue: gen double t = _n
replace t = mod(_n-1,8)+1
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
* a distinctive pre-trend in cohort 1 only -- the cohort whose label leg A sets to 0
replace y = y + 0.9*(t-5) if cohort == 1 & t < 5
tempfile P
save "`P'"

foreach leg in A B {
  * ---- binary design under the DOCUMENTED safety(off) bypass ----
  use "`P'", clear
  if "`leg'" == "A" gen double id = cond(cohort == 1, 0, cohort)
  if "`leg'" == "B" gen double id = cohort + 1
  dis "LEGSTART_`leg'_bin"
  cap noi wooldid y id t ttre, att att_i att_itime att_it esprelength(3) espostlength(2) safety(off)
  dis "LEG_`leg'_bin_RC=" _rc
  dis "MATSTART_`leg'_bin"
  cap noi mat li e(wdidmainresults), format(%21x)
  dis "MATEND_`leg'_bin"

  * ---- contreatwithin + contreat(): NO safety() option typed by the user ----
  use "`P'", clear
  if "`leg'" == "A" gen double id = cond(cohort == 1, 0, cohort)
  if "`leg'" == "B" gen double id = cohort + 1
  dis "LEGSTART_`leg'_ctw"
  cap noi wooldid y id t ttre, att att_i att_itime att_it esprelength(3) espostlength(2) contreat(dose) contreatcontrols(byCohort) contreatwithin
  dis "LEG_`leg'_ctw_RC=" _rc
  dis "MATSTART_`leg'_ctw"
  cap noi mat li e(wdidmainresults), format(%21x)
  dis "MATEND_`leg'_ctw"

  * ---- the SAME binary design with NO safety(off): the guard that should have fired ----
  use "`P'", clear
  if "`leg'" == "A" gen double id = cond(cohort == 1, 0, cohort)
  if "`leg'" == "B" gen double id = cohort + 1
  dis "LEGSTART_`leg'_guard"
  cap noi wooldid y id t ttre, att att_i att_itime att_it esprelength(3) espostlength(2)
  dis "LEG_`leg'_guard_RC=" _rc
}

dis "ALLDONE"
