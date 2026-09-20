*=====================================================================
* repro_FN03.do   SELF-CONTAINED.  No project data, no project generator.
*
* CLAIM: with contreat() + subgroups() + an event-study window, `fastnose'
*        REPORTS EVENT-STUDY ROWS THE FULL RUN DROPS.  The extra rows are the
*        regression's structural zeros, printed as average marginal effects
*        (value exactly 0), and they also enter e(b_es_ame_s#) / e(numpre_ame_s#).
*        sthlp:292 says the opposite: "fast mode therefore reports the same
*        event-study rows as the full run".
*
* MECHANISM: ado:2912, the BINARY subgroup-by-horizon margins call, threads
*        `sgmarginsdof', which BS-26 forces to df(`df_r') under fastnose precisely so
*        margins' not-estimable marker survives for the SI-7 drop.  Its CONTINUOUS
*        twin at ado:3290 still threads `marginsdof', which IS `nose' under fastnose,
*        so wooldid_sglbl has no marker to filter on and nothing is dropped.
*
* RUN:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FN03.do
*       then read repro_FN03.log
*=====================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

set seed 8181
set sortseed 8181

*---------------- a plain staggered panel with a continuous treatment intensity -------
* 90 units in 3 treated cohorts + never-treated; T = 11; dose varies across units.
clear
set obs 120
gen long unit = _n
gen long cohort = 1 + mod(_n-1, 4)          // cohort 4 = never treated
gen double ttre = .
replace ttre = 2003 if cohort == 1
replace ttre = 2006 if cohort == 2
replace ttre = 2009 if cohort == 3
gen double a_i = rnormal()
gen double d_i = 0.5 + runiform()           // unit-level dose intensity
expand 11
bys unit: gen int t = 1999 + _n
gen int h = t - ttre
gen byte post = (ttre < . & t >= ttre)
gen double dose = 0
replace dose = d_i if post == 1
gen double yr = 0.4*(t - 1999) + rnormal()*0.3
gen double y = a_i + yr + 1.2*dose*post
* a REALISTIC ragged subgroup: a policy variant only the 2009 cohort ever adopted
gen byte sgvar = 1
replace sgvar = 2 if ttre == 2009
tab sgvar, m

tempfile P
qui save `P'

*---------------- ARM 1: subgroups() -- the divergence --------------------------------
foreach mode in ANALYTIC FASTNOSE {
  local ex
  if "`mode'" == "FASTNOSE" local ex fastnose
  use `P', clear
  ereturn clear
  dis _n "==================== subgroups(), `mode' ===================="
  cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) ///
      fe(unit) subgroups(sgvar) esprelength(4) espostlength(4) `ex'
  dis "SG_`mode'_RC=" _rc
  foreach s in 1 2 {
    tempname B
    cap matrix `B' = e(b_es_ame_s`s')
    if _rc == 0 {
      dis "SG_`mode'_b_es_ame_s`s'  colnames: `: colnames `B''"
      matrix list `B', format(%12.0g)
    }
    else dis "SG_`mode'_b_es_ame_s`s'_ABSENT"
    dis "SG_`mode'_numpre_ame_s`s'=" e(numpre_ame_s`s')
    cap matrix drop `B'
    tempname R
    cap matrix `R' = e(wdidsgresults_`s')
    if _rc == 0 dis "SG_`mode'_rows_s`s': `: rowfullnames `R''"
    cap matrix drop `R'
  }
}

*---------------- ARM 2: same spec WITHOUT subgroups() (control) ----------------------
foreach mode in ANALYTIC FASTNOSE {
  local ex
  if "`mode'" == "FASTNOSE" local ex fastnose
  use `P', clear
  ereturn clear
  dis _n "==================== no subgroups(), `mode' ===================="
  cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) ///
      fe(unit) esprelength(4) espostlength(4) `ex'
  dis "NOSG_`mode'_RC=" _rc
  tempname B
  cap matrix `B' = e(b_es_ame)
  if _rc == 0 {
    dis "NOSG_`mode'_b_es_ame colnames: `: colnames `B''"
    matrix list `B', format(%12.0g)
  }
  else dis "NOSG_`mode'_b_es_ame_ABSENT"
  dis "NOSG_`mode'_numpre_ame=" e(numpre_ame)
  cap matrix drop `B'
  tempname R
  cap matrix `R' = e(wdidmainresults)
  if _rc == 0 dis "NOSG_`mode'_rows: `: rowfullnames `R''"
  cap matrix drop `R'
}

*---------------- ARM 3: same, but contreatcontrols(basic) ---------------------------
foreach mode in ANALYTIC FASTNOSE {
  local ex
  if "`mode'" == "FASTNOSE" local ex fastnose
  use `P', clear
  ereturn clear
  dis _n "==================== subgroups() + contreatcontrols(basic), `mode' ===================="
  cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
      fe(unit) subgroups(sgvar) esprelength(4) espostlength(4) `ex'
  dis "BAS_`mode'_RC=" _rc
  foreach s in 1 2 {
    tempname B
    cap matrix `B' = e(b_es_ame_s`s')
    if _rc == 0 {
      dis "BAS_`mode'_b_es_ame_s`s'  colnames: `: colnames `B''"
      matrix list `B', format(%12.0g)
    }
    else dis "BAS_`mode'_b_es_ame_s`s'_ABSENT"
    dis "BAS_`mode'_numpre_ame_s`s'=" e(numpre_ame_s`s')
    cap matrix drop `B'
  }
}
dis "DONE_REPRO_FN03"
