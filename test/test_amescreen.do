* test_amescreen.do — certification of the AMESE R1/R3/R4 screens (standing-issues
* family A; design note review_process/amese/R1_design_note.md).
* Scenarios: the structural cell-containment predicate fires (demo design, the
* referee's partial-containment mixed-cluster acid data, subgroups, permute observed
* run, byCohort-with-R4-refit) and must NOT fire (robust, cluster(unit), impvalues,
* fastnose, binary-only). Assertions are on the e() surface: fired runs carry missing
* ame SEs with points intact, e(b_es_ame) WITHOUT e(V_es_ame), and no ame joint/
* adversarial scalars; unfired runs carry real ame SEs and a posted V_es_ame.
* Run from test/:  "<stata>" -e do test_amescreen.do   -> expect AMESCREEN_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 2010
set sortseed 2010

local REPO ".."
do "`REPO'/wooldidinstall/wooldid.ado"

global pf 0

cap program drop mkdemo
program define mkdemo
  clear
  set seed 2010
  qui set obs 30
  gen i = _n
  gen ttre = cond(mod(_n,2)==0, 6, .)
  qui expand 6
  bysort i: gen unit = _n
  gen unitid = (i-1)*6 + unit
  gen dose = exp(0.8*runiform() + 0.2*mod(i,5))
  qui expand 12
  bysort i unitid: gen t = _n
  gen y = 0.1*i + 0.05*t + rnormal()
  sort i unitid t
  ereturn clear
end

cap program drop assertfire
program define assertfire
  args tag wantfire
  local bad 0
  tempname W
  cap matrix `W' = e(wdidmainresults)
  if _rc != 0 {
    dis as error "`tag' FAIL: wdidmainresults absent"
    global pf = $pf + 1
    exit
  }
  local ra = rownumb(`W', "main:att")
  local rm = rownumb(`W', "main:ame")
  local attse = cond(`ra' < ., `W'[`ra', 2], .)
  local amese = cond(`rm' < ., `W'[`rm', 2], .)
  cap confirm matrix e(b_es_ame)
  local brc = _rc
  cap confirm matrix e(V_es_ame)
  local vrc = _rc
  if `wantfire' {
    if !missing(`attse') & `attse' > 0 & missing(`amese') & `brc' == 0 & `vrc' != 0 & "`e(joint_Choriz)'" == "" & "`e(adversarial_pval_ame)'" == "" dis "`tag' ok"
    else {
      dis as error "`tag' FAIL (want fire): attse=`attse' amese=`amese' b_es rc=`brc' V_es rc=`vrc' jointC=[`e(joint_Choriz)'] advame=[`e(adversarial_pval_ame)']"
      global pf = $pf + 1
    }
  }
  else {
    if !missing(`attse') & !missing(`amese') & `amese' > 0 & `brc' == 0 & `vrc' == 0 dis "`tag' ok"
    else {
      dis as error "`tag' FAIL (want no fire): attse=`attse' amese=`amese' b_es rc=`brc' V_es rc=`vrc'"
      global pf = $pf + 1
    }
  }
end

* A1: demo design + cluster(i) — fires, tests suppressed, points intact
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esf espre(2) espost(2) adversarial jointtests safety(off)
if _rc != 0 {
  dis as error "A1 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire A1 1

* A2: same + robust — no fire, adversarial ame p posts
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) robust esf espre(2) espost(2) adversarial jointtests
if _rc != 0 {
  dis as error "A2 FAIL rc=" _rc
  global pf = $pf + 1
}
else {
  assertfire A2 0
  if "`e(adversarial_pval_ame)'" == "" {
    dis as error "A2b FAIL: adversarial ame p missing under robust"
    global pf = $pf + 1
  }
}

* A3: cluster(unitid) strictly finer than cells — no fire
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(unitid) esf espre(2) espost(2)
if _rc != 0 {
  dis as error "A3 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire A3 0

* A4: the referee's mixed-cluster data (PARTIAL containment; full-rank V pre-R1,
* spurious p 1.2e-08) — the ANY-cell semantics must fire here
use "`REPO'/review_process/amese/probes/refA.dta", clear
ereturn clear
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(cmix) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "A4 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire A4 1

* A5: byCohort + cluster(i) — Stata's V-zeroing screen fires, so the whole covariance
* comes back zero. The explicit-dummies refit that used to rescue the surviving att SE
* was removed by author ruling on 2026-09-08 (it re-solved the specification over
* hundreds of collinear dummy columns to salvage part of an inference the containment
* had already invalidated, delaying the user's discovery that the specification was
* broken). Expect now: rc 0, the run completes, EVERY standard error is missing, and
* the new caution says so in as many words.
mkdemo
qui log using ames_a5.txt, text replace name(ames)
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(byCohort) cluster(i) esf espre(2) espost(2)
local a5rc = _rc
qui log close ames
cap program drop scana5
program define scana5, rclass
  args lf needle
  tempname fh
  local found 0
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') != 0 local found 1
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
end
scana5 ames_a5.txt "all-zero variance-covariance matrix, so no standard error in this run is identified"
local a5msg = r(found)
scana5 ames_a5.txt "re-estimating through the explicit-dummies fallback"
local a5old = r(found)
cap erase ames_a5.txt
tempname A5W
local a5attse .
local a5amese .
if `a5rc' == 0 {
  cap matrix `A5W' = e(wdidmainresults)
  if _rc == 0 {
    local a5ra = rownumb(`A5W', "main:att")
    local a5rm = rownumb(`A5W', "main:ame")
    local a5attse = cond(`a5ra' < ., `A5W'[`a5ra', 2], .)
    local a5amese = cond(`a5rm' < ., `A5W'[`a5rm', 2], .)
  }
}
dis "A5 rc=`a5rc' attse=`a5attse' amese=`a5amese' newcaution=`a5msg' oldrefit=`a5old'"
if `a5rc' == 0 & missing(`a5attse') & missing(`a5amese') & `a5msg' == 1 & `a5old' == 0 {
  dis "A5 ok"
}
else {
  dis as error "A5 FAIL: the all-zero-VCE run should complete with every SE missing and the new caution, and must not refit"
  global pf = $pf + 1
}

* A6: fastnose — predicate must stay off (no SEs anywhere by design, rc 0)
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esf espre(2) espost(2) fastnose
if _rc != 0 {
  dis as error "A6 FAIL rc=" _rc
  global pf = $pf + 1
}
else dis "A6 ok"

* A7: permute end-to-end — perm outputs intact while the observed run's ame display
* is scrubbed
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esf espre(2) permute(shuffle) permreps(25)
if _rc != 0 {
  dis as error "A7 FAIL rc=" _rc
  global pf = $pf + 1
}
else {
  if e(perm_completed) == 25 & !missing(e(perm_p_att)) dis "A7 ok"
  else {
    dis as error "A7 FAIL: perm_completed=" e(perm_completed) " perm_p=" e(perm_p_att)
    global pf = $pf + 1
  }
}

* A8: impvalues — never suppressed (cross-cluster dispersion inference is valid here)
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esf espre(2) impvalues
if _rc != 0 {
  dis as error "A8 FAIL rc=" _rc
  global pf = $pf + 1
}
else {
  tempname W
  matrix `W' = e(wdidmainresults)
  local rm = rownumb(`W', "main:ame")
  if `rm' < . & !missing(`W'[`rm', 2]) & `W'[`rm', 2] > 0 dis "A8 ok"
  else {
    dis as error "A8 FAIL: impvalues ame SE missing"
    global pf = $pf + 1
  }
}

* A9: binary-only cluster(i) — untouched (no contreat, no predicate)
mkdemo
cap noi wooldid y i t ttre, att cluster(i) esf espre(2) espost(2) jointtests safety(off)
if _rc != 0 {
  dis as error "A9 FAIL rc=" _rc
  global pf = $pf + 1
}
else {
  tempname W
  matrix `W' = e(wdidmainresults)
  local ra = rownumb(`W', "main:att")
  if `ra' < . & !missing(`W'[`ra', 2]) dis "A9 ok"
  else {
    dis as error "A9 FAIL: binary att SE missing"
    global pf = $pf + 1
  }
}

* A10: subgroups — fires cleanly (rc 0), sg ame surfaces scrubbed, att rows intact
mkdemo
gen sg = mod(i,2)+1
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests safety(off)
if _rc != 0 {
  dis as error "A10 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire A10 1

* ---- GROUP B: the byCohort point-identification screen ----

* B1: a treated cohort reduced to one unit (zero dose variation in its reference
* pool) must be REFUSED with r459, naming the cohort
mkdemo
drop if i == 2 & unit > 1
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(byCohort) esf espre(2) espost(2)
if _rc == 459 dis "B1 ok"
else {
  dis as error "B1 FAIL: rc=" _rc " (want 459)"
  global pf = $pf + 1
}

* B2: thin variation (two units, dose spread ~1e-7) -> fragility caution, run completes
mkdemo
drop if i == 2 & unit > 2
qui replace dose = 1 + 1e-7*(unit == 2) if i == 2
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(byCohort) esf espre(2) espost(2)
if _rc == 0 dis "B2 ok"
else {
  dis as error "B2 FAIL: rc=" _rc " (want 0 with the fragility caution)"
  global pf = $pf + 1
}

* B3: healthy byCohort + robust -> clean run, real ame inference, screen silent
mkdemo
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(byCohort) esf espre(2) espost(2)
if _rc != 0 {
  dis as error "B3 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire B3 0

* C1 (SI-A4/R6): cohort-constant dose (the d3 class) -> named refusal of the ame
* block: the run completes rc 0, ame POINT estimates are voided to missing (never
* the regression's fabricated structural zeros), att untouched with real SEs
mkdemo
qui replace dose = 0.5 + 0.1*i
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(basic) esf espre(2) espost(2)
if _rc != 0 {
  dis as error "C1 FAIL rc=" _rc
  global pf = $pf + 1
}
else {
  tempname W
  matrix `W' = e(wdidmainresults)
  local ra = rownumb(`W', "main:att")
  local rm = rownumb(`W', "main:ame")
  local attb = cond(`ra' < ., `W'[`ra', 1], -99)
  local attse = cond(`ra' < ., `W'[`ra', 2], -99)
  local ameb = cond(`rm' < ., `W'[`rm', 1], -99)
  if !missing(`attb') & !missing(`attse') & `attse' > 0 & missing(`ameb') dis "C1 ok"
  else {
    dis as error "C1 FAIL: attb=`attb' attse=`attse' ameb=`ameb' (want ame b missing with att real)"
    global pf = $pf + 1
  }
}

* C2 (SI-A4 untreatedZero exception): expected zero doses in the untreated periods
* must NOT trigger the void - the treated-cell slopes are live and the ame is real
mkdemo
qui replace dose = 0 if ttre >= . | t < ttre
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(untreatedZero) esf espre(2) espost(2)
if _rc != 0 {
  dis as error "C2 FAIL rc=" _rc
  global pf = $pf + 1
}
else {
  tempname W
  matrix `W' = e(wdidmainresults)
  local rm = rownumb(`W', "main:ame")
  local ameb = cond(`rm' < ., `W'[`rm', 1], .)
  local amese = cond(`rm' < ., `W'[`rm', 2], .)
  if !missing(`ameb') & !missing(`amese') & `amese' > 0 dis "C2 ok"
  else {
    dis as error "C2 FAIL: ame b=`ameb' se=`amese' (want real ame under untreatedZero)"
    global pf = $pf + 1
  }
}

if $pf == 0 dis "AMESCREEN_ALL_PASS"
else dis as error "AMESCREEN_FAILS=$pf"
exit, clear
