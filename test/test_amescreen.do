* test_amescreen.do — certification of the AMESE R1/R3/R4 screens (standing-issues
* family A; design note review_process/amese/R1_design_note.md).
* Scenarios: the structural cell-containment predicate fires (demo design, the
* referee's partial-containment mixed-cluster acid data, subgroups, permute observed
* run, byCohort-with-R4-refit) and must NOT fire (robust, cluster(unit), impvalues,
* fastnose, binary-only). Assertions are on the e() surface: fired runs carry missing
* ame SEs with points intact, e(b_es_ame) WITHOUT e(V_es_ame), and no ame joint/
* adversarial scalars; unfired runs carry real ame SEs and a posted V_es_ame.
* Run from a fresh subdirectory of test/:  "<stata>" -e do ../test_amescreen.do
*   -> expect AMESCREEN_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 2010
set sortseed 2010

local REPO "../.."
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
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esf espre(2) espost(2) adversarial jointtests
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

* A5: byCohort + cluster(i) — Stata's V-zeroing screen fires; the R4 explicit-dummies
* refit must recover a positive att SE (formerly a 20-of-20 blackout) with ame scrubbed
mkdemo
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(byCohort) cluster(i) esf espre(2) espost(2)
if _rc != 0 {
  dis as error "A5 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire A5 1

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
cap noi wooldid y i t ttre, att cluster(i) esf espre(2) espost(2) jointtests
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
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests
if _rc != 0 {
  dis as error "A10 FAIL rc=" _rc
  global pf = $pf + 1
}
else assertfire A10 1

if $pf == 0 dis "AMESCREEN_ALL_PASS"
else dis as error "AMESCREEN_FAILS=$pf"
exit, clear
