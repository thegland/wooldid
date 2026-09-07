* test_edges.do - edge-of-panel shapes (bug sweep BS-33, 2026-09-06).
* Run from test/: stata-mp -e do test_edges.do -> expect EDGES_ALL_PASS
*
* E1 panel A (the onset period dropped for the treated cohort, a common convention):
*    bare att -> rc 0 with att = 1.5 (the noise-free DGP effect); att espostlength(2) -> rc 0,
*    att = 1.5, and the only ES row is relative period 1 (0 is not realized); no pre-att row
* E2 panel B (biennial t = 1,3,..,11 with onset 6, so relative period 0 never occurs):
*    bare att -> rc 0, att = 1.5; espostlength(2) -> rc 0 with the single ES row at 1
* E3 panel B espostlength(1): the requested window is relative period 0 only, which nothing
*    realizes -> the informative error 459, not the bare r(198) from matrix ()
* E4 panel A esprelength(1) espostlength(2): rc 0 (relative period -1 survives; 4 rows)
* E5 regression pin: dta_binary att cluster(i) esprelength(2) espostlength(2) is unmoved
*    (att .8235285, se .0629643, six rows: att, pre att, pre2 pre1 post0 post1; e(wdid_espre) 2)
clear all
set more off
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop rowk
program rowk, rclass
  args mname rname
  local rn : rownames `mname'
  local k : list posof "`rname'" in rn
  return scalar k = `k'
end

* panel A: 12 units x 6 periods, onset 4 for units 1-4, onset row dropped
clear
set obs 12
gen i = _n
gen ttre = cond(i <= 4, 4, .)
expand 6
bysort i: gen t = _n
gen y = 1.5*(t >= ttre & ttre < .) + 0.3*i + 0.2*t
drop if i <= 4 & t == 4
tempfile pa
save `pa'

* panel B: biennial periods, onset 6 (never a realized period)
clear
set obs 12
gen i = _n
gen ttre = cond(i <= 4, 6, .)
expand 6
bysort i: gen t = 2*_n - 1
gen y = 1.5*(t >= ttre & ttre < .) + 0.3*i + 0.2*t
tempfile pb
save `pb'

* ---- E1
use `pa', clear
cap noi wooldid y i t ttre, att
local rc1a = _rc
local att1a = .
if `rc1a' == 0 {
  matrix E1A = e(wdidmainresults)
  rowk E1A att
  local att1a = E1A[r(k), 1]
}
use `pa', clear
cap noi wooldid y i t ttre, att espostlength(2)
local rc1b = _rc
local att1b = .
local rn1b
if `rc1b' == 0 {
  matrix E1B = e(wdidmainresults)
  rowk E1B att
  local att1b = E1B[r(k), 1]
  local rn1b : rownames E1B
}
di "E1 bare rc=`rc1a' att=`att1a' | espostlength(2) rc=`rc1b' att=`att1b' rows=[`rn1b']"
if `rc1a' != 0 | abs(`att1a' - 1.5) > 1e-6 | `rc1b' != 0 | abs(`att1b' - 1.5) > 1e-6 | rowsof(E1B) != 2 {
  di as error "E1 FAIL: the onset-period gap still breaks a bare att call"
  local ++fails
}

* ---- E2
use `pb', clear
cap noi wooldid y i t ttre, att
local rc2a = _rc
local att2a = .
if `rc2a' == 0 {
  matrix E2A = e(wdidmainresults)
  rowk E2A att
  local att2a = E2A[r(k), 1]
}
use `pb', clear
cap noi wooldid y i t ttre, att espostlength(2)
local rc2b = _rc
local att2b = .
local rn2b
if `rc2b' == 0 {
  matrix E2B = e(wdidmainresults)
  rowk E2B att
  local att2b = E2B[r(k), 1]
  local rn2b : rownames E2B
}
di "E2 bare rc=`rc2a' att=`att2a' | espostlength(2) rc=`rc2b' att=`att2b' rows=[`rn2b']"
if `rc2a' != 0 | abs(`att2a' - 1.5) > 1e-6 | `rc2b' != 0 | abs(`att2b' - 1.5) > 1e-6 | rowsof(E2B) != 2 {
  di as error "E2 FAIL: the off-grid onset still breaks att"
  local ++fails
}

* ---- E3
use `pb', clear
cap noi wooldid y i t ttre, att espostlength(1)
local rc3 = _rc
di "E3 espostlength(1) on the biennial panel rc=`rc3' (expect the informative 459)"
if `rc3' != 459 {
  di as error "E3 FAIL: an unrealizable window should give the informative error 459, got `rc3'"
  local ++fails
}

* ---- E4
use `pa', clear
cap noi wooldid y i t ttre, att esprelength(1) espostlength(2)
local rc4 = _rc
local att4 = .
local rn4
if `rc4' == 0 {
  matrix E4 = e(wdidmainresults)
  rowk E4 att
  local att4 = E4[r(k), 1]
  local rn4 : rownames E4
}
di "E4 rc=`rc4' att=`att4' rows=[`rn4']"
if `rc4' != 0 | abs(`att4' - 1.5) > 1e-6 | rowsof(E4) != 4 {
  di as error "E4 FAIL: the surviving pre period should still run"
  local ++fails
}

* ---- E5
use data/dta_binary.dta, clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
local rc5 = _rc
local b5 = .
local se5 = .
local rn5
if `rc5' == 0 {
  matrix E5 = e(wdidmainresults)
  rowk E5 att
  local b5 = E5[r(k), 1]
  local se5 = E5[r(k), 2]
  local rn5 : rownames E5
  local espre5 = e(wdid_espre)
}
di "E5 rc=`rc5' att=" %10.7f `b5' " se=" %10.7f `se5' " rows=[`rn5'] e(wdid_espre)=`espre5'"
if `rc5' != 0 | abs(`b5' - .8235285) > 5e-8 | abs(`se5' - .0629643) > 5e-8 | rowsof(E5) != 6 | `espre5' != 2 {
  di as error "E5 FAIL: dta_binary moved"
  local ++fails
}


cap program drop scanE
program define scanE, rclass
  args lf needle
  tempname fh
  local found 0
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') != 0 local ++found
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
end
set linesize 255

* ---- E6 (bug sweep BS-2): the forced-reference caution fires at esprelength(1)
* panel A's treated units with onset 2: the requested pre1 horizon is the panel's first
* period, so the reference is forced onto a requested horizon. The caution must print at
* esprelength(1) (it was gated off there), not at esprelength(0), and not on an
* untruncated onset; the fast branch (safety(off)) gives the identical att
clear
set obs 12
gen i = _n
gen ttre = cond(i <= 4, 2, .)
expand 6
bysort i: gen t = _n
gen y = 1.5*(t >= ttre & ttre < .) + 0.3*i + 0.2*t
tempfile pe6
save `pe6'
local e6c
local e6att
foreach leg in "esprelength(1)" "esprelength(1) safety(off)" "esprelength(0)" {
  use `pe6', clear
  qui log using edges_e6.txt, text replace name(ed)
  cap noi wooldid y i t ttre, att `leg' espostlength(1)
  local rc = _rc
  qui log close ed
  scanE edges_e6.txt "reference period for some units forced"
  local e6c `e6c' `r(found)'
  matrix E6 = e(wdidmainresults)
  rowk E6 att
  local e6att `e6att' `=string(E6[r(k),1], "%16.0g")'
  if `rc' != 0 {
    di as error "E6 FAIL: rc=`rc' under `leg'"
    local ++fails
  }
}
use `pe6', clear
replace ttre = 4 if ttre < .
qui log using edges_e6.txt, text replace name(ed)
cap noi wooldid y i t ttre, att esprelength(1) espostlength(1)
qui log close ed
scanE edges_e6.txt "reference period for some units forced"
local e6c `e6c' `r(found)'
cap erase edges_e6.txt
local a1 : word 1 of `e6att'
local a2 : word 2 of `e6att'
di "E6 caution counts (pre1, pre1 safety(off), pre0, onset-4) = `e6c'; att pre1 = `a1' vs fast branch `a2'"
if "`e6c'" != "1 0 0 0" | "`a1'" != "`a2'" {
  di as error "E6 FAIL: the forced-reference caution should fire exactly at esprelength(1) with the att unchanged"
  local ++fails
}

* ---- E7 (BS-14): unconditionalse with a multiway cluster() is refused by name
use data/dta_binary.dta, clear
gen c2 = mod(i, 4)
qui log using edges_e7.txt, text replace name(ed)
cap noi wooldid y i t ttre, att cluster(i c2) unconditionalse
local rc7 = _rc
qui log close ed
scanE edges_e7.txt "unconditionalse supports a single cluster() variable"
local e7a = r(found)
cap erase edges_e7.txt
use data/dta_binary.dta, clear
cap noi wooldid y i t ttre, att cluster(i) unconditionalse
local rc7b = _rc
di "E7 multiway+unconditionalse rc=`rc7' named=`e7a' | single cluster rc=`rc7b'"
if `rc7' != 198 | `e7a' != 1 | `rc7b' != 0 {
  di as error "E7 FAIL: the unconditionalse multiway gate"
  local ++fails
}

* ---- E8 (BS-29): the levels exemption's consequence - the subgroup adversarial p is
* identical collapsed and on the full data, with no fallback
use data/dta_binary.dta, clear
gen cl8 = mod(i, 4)
cap noi wooldid y i t ttre, att cluster(cl8) subgroups(sg) esfixedbaseperiod esprelength(4) espostlength(4) adversarial
local rc8 = _rc
local p8a = e(adversarial_pval_att_s1)
local g8 `e(wdid_fcguard)'
use data/dta_binary.dta, clear
gen cl8 = mod(i, 4)
cap noi wooldid y i t ttre, att cluster(cl8) subgroups(sg) esfixedbaseperiod esprelength(4) espostlength(4) adversarial fullprecisionslow
local rc8b = _rc
local p8b = e(adversarial_pval_att_s1)
di "E8 adversarial p (sg1) default=`p8a' fullprecisionslow=`p8b' fcguard=[`g8'] rc=`rc8'/`rc8b'"
if `rc8' != 0 | `rc8b' != 0 | "`p8a'" != "`p8b'" | "`g8'" != "" | "`p8a'" == "" {
  di as error "E8 FAIL: the subgroup adversarial p should be identical collapsed and uncollapsed with no fallback"
  local ++fails
}

* ---- E9 (BS-32): a failure after the post (a saveplots() name clash in the plot stage)
* leaves e(cmd) posted, so the documented replay refusal is reachable
use `pe6', clear
replace ttre = 4 if ttre < .
cap erase v1_ES.dta
cap erase v1_ES.gph
qui wooldid y i t ttre, att esprelength(1) espostlength(1) makeplots saveplots(v1)
cap noi wooldid y i t ttre, att esprelength(1) espostlength(1) makeplots saveplots(v1)
local rc9 = _rc
local e9cmd `e(cmd)'
qui log using edges_e9.txt, text replace name(ed)
cap noi wooldid
local rc9b = _rc
qui log close ed
scanE edges_e9.txt "predate that failed call"
local e9r = r(found)
cap erase edges_e9.txt
local e9files : dir . files "v1*"
foreach f of local e9files {
  cap erase "`f'"
}
di "E9 second saveplots(v1) rc=`rc9' e(cmd)=[`e9cmd'] replay rc=`rc9b' refusal named=`e9r'"
if `rc9' == 0 | "`e9cmd'" != "wooldid" | `rc9b' != 198 | `e9r' != 1 {
  di as error "E9 FAIL: e(cmd) must survive a post-estimation failure so replay refuses by name"
  local ++fails
}
* a fresh successful call clears the refusal
use `pe6', clear
replace ttre = 4 if ttre < .
qui wooldid y i t ttre, att
cap wooldid
if _rc != 0 {
  di as error "E9b FAIL: replay after a successful call rc=`=_rc'"
  local ++fails
}

* ---- E10 (BS-35): one cluster stops with the cluster count, not r(125) under a
* continuous-treatment caution; fastnose, robust and two clusters run
use `pe6', clear
replace ttre = 4 if ttre < .
gen c1 = 1
gen c2 = (i <= 6) + 1
local e10
foreach leg in "cluster(c1)" "cluster(c1) fastnose" "robust" "cluster(c2)" {
  use `pe6', clear
  replace ttre = 4 if ttre < .
  gen c1 = 1
  gen c2 = (i <= 6) + 1
  qui log using edges_e10.txt, text replace name(ed)
  cap noi wooldid y i t ttre, att `leg'
  local rc = _rc
  qui log close ed
  scanE edges_e10.txt "cluster-robust inference needs at least 2"
  local m = r(found)
  scanE edges_e10.txt "continuous"
  local c = r(found)
  local e10 `e10' [`leg': rc=`rc' msg=`m' continuous=`c']
  if "`leg'" == "cluster(c1)" & (`rc' != 459 | `m' != 1) local ++fails
  if "`leg'" != "cluster(c1)" & `rc' != 0 local ++fails
  if `c' > 0 local ++fails
}
cap erase edges_e10.txt
di "E10 `e10'"

* ---- E11 (bug sweep BS-3): interactivecontrols() inert -> named caution
* fine i (the unit) with a unit-level covariate: every cell slope drops, the caution names
* xtv, and att equals the controls(i.t#c.xtv) run; coarse i (pairs) gives within-cell
* variation, no caution and a different att; coarse i with a covariate constant across
* units within a period (period-level, so time-VARYING) still has no within-cell
* variation and must draw the caution - the screen keys on cell variation, not on time
clear
set obs 12
gen unit = _n
gen g = ceil(unit/2)
gen ttre = cond(unit <= 4, 4, .)
expand 8
bysort unit: gen t = _n
set seed 311
gen xtv = 0.3*unit + 0.1*t + rnormal()
bysort t (unit): gen xt = 0.2*t + rnormal() if _n == 1
bysort t (unit): replace xt = xt[1]
gen y = 1.5*(t >= ttre & ttre < .) + 0.3*unit + 0.2*t + 0.4*xtv + rnormal()
tempfile pe11
save `pe11'
local e11
local e11att
foreach leg in "unit interactivecontrols(xtv)" "unit controls(i.t#c.xtv)" "g interactivecontrols(xtv)" "g controls(i.t#c.xtv)" "g interactivecontrols(xt)" {
  use `pe11', clear
  local iv : word 1 of `leg'
  local opts : subinstr local leg "`iv' " ""
  qui log using edges_e11.txt, text replace name(ed)
  cap noi wooldid y `iv' t ttre, att `opts'
  local rc = _rc
  qui log close ed
  scanE edges_e11.txt "slope term for the interactive control"
  local c = r(found)
  cap erase edges_e11.txt
  local a .
  if `rc' == 0 {
    matrix E11 = e(wdidmainresults)
    rowk E11 att
    local a = E11[r(k), 1]
  }
  local e11 `e11' `c'
  local e11att `e11att' `=string(`a', "%16.0g")'
  if `rc' != 0 {
    di as error "E11 FAIL: rc=`rc' under [`leg']"
    local ++fails
  }
}
tokenize `e11att'
di "E11 cautions (fine ic, fine ctrl, coarse ic, coarse ctrl, coarse period-level ic) = `e11'; att fine ic=`1' vs ctrl=`2' | coarse ic=`3' vs ctrl=`4'"
if "`e11'" != "1 0 0 0 1" | abs(`1' - `2') > 1e-8 | abs(`3' - `4') < 1e-6 {
  di as error "E11 FAIL: the interactivecontrols inertness caution"
  local ++fails
}

* ---- E12 (BS-4): a ppmlhdfe separation is reported as separation, not multicollinearity
clear
set obs 12
gen i = _n
gen ttre = cond(i <= 4, 4, .)
expand 6
bysort i: gen t = _n
set seed 412
gen y = rpoisson(exp(0.5 + 0.3*(t >= ttre & ttre < .) + 0.05*i + 0.02*t))
replace y = 0 if i == 1 & t == 5
tempfile pe12
save `pe12'
qui log using edges_e12.txt, text replace name(ed)
cap noi wooldid y i t ttre, att poisson
local rc12 = _rc
qui log close ed
scanE edges_e12.txt "perfect separation"
local e12a = r(found)
scanE edges_e12.txt "dropped for multicollinearity"
local e12b = r(found)
scanE edges_e12.txt "contreatwithin"
local e12c = r(found)
cap erase edges_e12.txt
use `pe12', clear
qui log using edges_e12.txt, text replace name(ed)
cap noi wooldid y i t ttre, att poisson safety(off)
local rc12b = _rc
qui log close ed
scanE edges_e12.txt "treated observation(s) from the estimation sample (separation"
local e12d = r(found)
cap erase edges_e12.txt
di "E12 poisson separation: rc=`rc12' separation=`e12a' multicollinearity=`e12b' contreatwithin=`e12c' | safety(off) rc=`rc12b' excluded-caution=`e12d'"
if `rc12' != 322 | `e12a' != 1 | `e12b' != 0 | `e12c' != 0 | `rc12b' != 0 | `e12d' != 1 {
  di as error "E12 FAIL: the ppmlhdfe separation branch"
  local ++fails
}

if `fails' == 0 di as result "EDGES_ALL_PASS"
else di as error "EDGES_FAILS=`fails'"
