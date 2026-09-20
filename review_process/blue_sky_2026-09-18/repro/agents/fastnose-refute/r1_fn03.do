*=====================================================================
* r1_fn03.do  -- SKEPTIC's INDEPENDENT test of FASTNOSE-03.
*
* CLAIM UNDER TEST: with contreat() + subgroups() + an event-study window,
*   `fastnose' PUBLISHES event-study rows in the ame family that the identical
*   analytic run DROPS as non-estimable, with the extra rows carrying exactly 0.
*   sthlp:292: "fast mode therefore reports the same event-study rows as the full run".
*
* I do NOT run the prober's do-file.  My own panel, my own geometry:
*   105 units, T = 2010..2019, onsets 2013 / 2016 / 2019 + never-treated,
*   dose varies WITHIN cohort, and the subgroup is a late-adopter region
*   (sg==2 <=> the 2019 cohort) -- a shape an applied economist would have.
*   esprelength(3) espostlength(2).
*
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r1_fn03.do
*=====================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 31337
set sortseed 31337

*--------------------------- my panel ---------------------------------
clear
set obs 105
gen long unit = _n
gen byte grp  = 1 + mod(_n-1, 5)          // grps 4,5 -> never treated
gen double ttre = .
replace ttre = 2013 if grp == 1
replace ttre = 2016 if grp == 2
replace ttre = 2019 if grp == 3
gen long cohort = grp
gen double a_i  = rnormal()*0.8
gen double d_i  = 0.4 + 1.1*runiform()     // dose intensity varies WITHIN cohort
expand 10
bys unit: gen int t = 2009 + _n
gen byte post = (ttre < . & t >= ttre)
gen double dose = 0
replace dose = d_i if post == 1
gen double y = a_i + 0.35*(t-2009) + 0.9*dose*post + rnormal()*0.4
* late-adopter region: only the 2019 cohort ever adopted the policy variant
gen byte region = 1
replace region = 2 if ttre == 2019
tab region grp, m

tempfile P
qui save `P'

*--------------------------- driver -----------------------------------
capture program drop dumpsg
program define dumpsg
  args tag
  foreach s in 1 2 {
    tempname B
    cap matrix `B' = e(b_es_ame_s`s')
    if _rc == 0 {
      dis "`tag'_bES_ame_s`s'_COLS: `: colnames `B''"
      matrix list `B', format(%12.0g)
    }
    else dis "`tag'_bES_ame_s`s'_ABSENT"
    dis "`tag'_numpre_ame_s`s'=" e(numpre_ame_s`s')
    dis "`tag'_ndegen_ame_s`s'=" e(es_ndegen_ame_s`s')
    tempname R
    cap matrix `R' = e(wdidsgresults_`s')
    if _rc == 0 dis "`tag'_sgrows_s`s': `: rowfullnames `R''"
  }
  dis "`tag'_bES_att_PRESENT=" (0==_rc)
  tempname A
  cap matrix `A' = e(b_es_att)
  if _rc == 0 dis "`tag'_bES_att_COLS: `: colnames `A''"
  else dis "`tag'_bES_att_ABSENT"
  dis "`tag'_numpre_att=" e(numpre_att)
  dis "`tag'_vcvstate=" "`e(wdid_vcvstate)'"
end

*========== ARM B : contreatcontrols(basic) ===========================
foreach m in ANA FST {
  local ex
  if "`m'" == "FST" local ex fastnose
  use `P', clear
  ereturn clear
  dis _n "############ ARM basic / `m' ############"
  cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
      fe(unit) subgroups(region) esprelength(3) espostlength(2) `ex'
  dis "BASIC_`m'_RC=" _rc
  dumpsg BASIC_`m'
}

*========== ARM U : contreatcontrols(untreatedZero) ===================
foreach m in ANA FST {
  local ex
  if "`m'" == "FST" local ex fastnose
  use `P', clear
  ereturn clear
  dis _n "############ ARM uzero / `m' ############"
  cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) ///
      fe(unit) subgroups(region) esprelength(3) espostlength(2) `ex'
  dis "UZERO_`m'_RC=" _rc
  dumpsg UZERO_`m'
}

*========== ARM C : control -- SAME spec, NO subgroups() ==============
foreach m in ANA FST {
  local ex
  if "`m'" == "FST" local ex fastnose
  use `P', clear
  ereturn clear
  dis _n "############ ARM nosg basic / `m' ############"
  cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
      fe(unit) esprelength(3) espostlength(2) `ex'
  dis "NOSG_`m'_RC=" _rc
  tempname B
  cap matrix `B' = e(b_es_ame)
  if _rc == 0 {
    dis "NOSG_`m'_bES_ame_COLS: `: colnames `B''"
    matrix list `B', format(%12.0g)
  }
  else dis "NOSG_`m'_bES_ame_ABSENT"
  dis "NOSG_`m'_numpre_ame=" e(numpre_ame)
}

dis "DONE_R1_FN03"
