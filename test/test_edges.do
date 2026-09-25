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
* E14 the if() qualifier's exact-one rule (bug sweep BS-37): "if flag" with flag in {1,2}
*    keeps only the flag == 1 rows, now with a caution naming the count and those rows
*    folded into e(N_dropped); "if flag != 0" and "if flag == 1" draw no caution; a
*    missing-valued flag leaves the sample and is deliberately NOT counted
* E15 the worker's markout on factor-variable terms: controls(ib2.g) and fe(ib2.g) run (they
*    used to die on Stata's bare "ib2:  operator invalid", r(198)), mark out the rows where g
*    is missing exactly as controls(i.g) and fe(i.g) do, and give the same att
* E16 a call failing after the stage-1 fit (the r(322) collinearity screen) leaves the caller's
*    e() in place (a prior regress), not wooldid's internal reghdfe fit
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
scanE edges_e7.txt "does not support multiway clustering"
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
cap noi wooldid y i t ttre, att cluster(cl8) subgroups(sg) esfixedbaseperiod esprelength(4) espostlength(4) adversarial nofastcollapse
local rc8b = _rc
local p8b = e(adversarial_pval_att_s1)
di "E8 adversarial p (sg1) default=`p8a' nofastcollapse=`p8b' fcguard=[`g8'] rc=`rc8'/`rc8b'"
if `rc8' != 0 | `rc8b' != 0 | "`p8a'" != "`p8b'" | "`g8'" != "" | "`p8a'" == "" {
  di as error "E8 FAIL: the subgroup adversarial p should be identical collapsed and uncollapsed with no fallback"
  local ++fails
}

* ---- E9 (BS-32): a failure after the post (a saveplots() name clash in the plot stage).
* Since the author's round-2 ruling (2026-09-23, "keep the promise") the failed call leaves e()
* as it was before the call: the first run's results, not the failed call's own post, which
* used to stay in e(). The second call adds jointtests so the two command lines differ and
* e(cmdline) says whose results e() holds
use `pe6', clear
replace ttre = 4 if ttre < .
cap erase v1_ES.dta
cap erase v1_ES.gph
qui wooldid y i t ttre, att esprelength(1) espostlength(1) makeplots saveplots(v1)
local e9cl0 `"`e(cmdline)'"'
cap noi wooldid y i t ttre, att esprelength(1) espostlength(1) makeplots saveplots(v1) jointtests
local rc9 = _rc
local e9cmd `e(cmd)'
if `"`e(cmdline)'"' != `"`e9cl0'"' {
  di as error `"E9a FAIL: after the failed call e(cmdline) is [`e(cmdline)'], not the earlier run's [`e9cl0']"'
  local ++fails
}
qui log using edges_e9.txt, text replace name(ed)
cap noi wooldid
local rc9b = _rc
qui log close ed
scanE edges_e9.txt "refit wooldid before replaying"
local e9r = r(found)
scanE edges_e9.txt "Wooldid Estimation"
local e9w = r(found)
cap erase edges_e9.txt
local e9files : dir . files "v1*"
foreach f of local e9files {
  cap erase "`f'"
}
di "E9 second saveplots(v1) rc=`rc9' e(cmd)=[`e9cmd'] replay rc=`rc9b' refusal named=`e9r' replayed=`e9w'"
* since the author's round-2 ruling ("replay a restored store") the earlier run's results, back
* in e(), replay after the failure; the replay used to refuse by name
if `rc9' == 0 | "`e9cmd'" != "wooldid" | `rc9b' != 0 | `e9r' != 0 | `e9w' != 1 {
  di as error "E9 FAIL: after a post-estimation failure the earlier results must be back in e() and replay"
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
  scanE edges_e10.txt "too few clusters"
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
scanE edges_e12.txt "dropped every observation of the treated"
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

* ---- E13 (BS-33 follow-up): the two causes of an empty horizon scan get different
* messages. E3 above covers the data cause (the requested period is never realized).
* Here the periods ARE realized but an absorbed fixed effect absorbs the cohort-by-period
* treatment interactions, so the fitted model carries no coefficient to find; the run must say
* so rather than blaming the data (it used to die with a bare invalid-syntax r(198)). The
* trigger was ccc_absorb(x2) on this unit-level i until the centred coarse controls
* (2026-09-25): its uncentred absorbed slope spanned each one-row cell's dummy. The centred
* slope is zero in a one-row cell, so that run now completes (E13b below), and the absorbing
* term here is an fe() with one level per treated window row
use suite/data/dta_s_covars.dta, clear
gen double cellfe = cond(ttre < . & t - ttre >= -1 & t - ttre <= 1, i*100 + t, 0)
qui log using edges_e13.txt, text replace name(ed)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) fe(cellfe)
local rc13 = _rc
qui log close ed
scanE edges_e13.txt "fitted model carries no estimable treatment-effect coefficient"
local e13a = r(found)
scanE edges_e13.txt "no treated unit is observed at any of the requested"
local e13b = r(found)
cap erase edges_e13.txt
di "E13 fe(cellfe) rc=`rc13' absorbed-message=`e13a' data-message=`e13b'"
if `rc13' != 459 | `e13a' != 1 | `e13b' != 0 {
  di as error "E13 FAIL: an absorbed treatment term must not be reported as missing data"
  local ++fails
}

* ---- E13b (centred coarse controls, ruling 2026-09-25): on this unit-level i every treated
* cell is one row, so the centred copy of the control is zero in every cell and is left out,
* and the control enters once as a pooled regressor. ccc_absorb(x2) and
* coarsecohortcontrols(c.x2) must then run and equal controls(c.x2) exactly
tempname e13m0 e13m1 e13m2
use suite/data/dta_s_covars.dta, clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) controls(c.x2)
local rc13c = _rc
if `rc13c' == 0 matrix `e13m0' = e(wdidmainresults)
use suite/data/dta_s_covars.dta, clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) ccc_absorb(x2)
local rc13a = _rc
local e13d1 .
if `rc13a' == 0 & `rc13c' == 0 local e13d1 = mreldif(e(wdidmainresults), `e13m0')
use suite/data/dta_s_covars.dta, clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) coarsecohortcontrols(c.x2)
local rc13k = _rc
local e13d2 .
if `rc13k' == 0 & `rc13c' == 0 local e13d2 = mreldif(e(wdidmainresults), `e13m0')
di "E13b one-row cells: rc controls=`rc13c' ccc_absorb=`rc13a' coarse=`rc13k' | mreldif ccc_absorb vs controls=`e13d1' coarse vs controls=`e13d2'"
if `rc13c' != 0 | `rc13a' != 0 | `rc13k' != 0 | `e13d1' > 1e-10 | `e13d2' > 1e-10 {
  di as error "E13b FAIL: with one-row cells ccc_absorb(x2) and coarsecohortcontrols(c.x2) must run and equal controls(c.x2)"
  local ++fails
}

* ---- E14 (BS-37): the exact-one if() rule is named and counted ----
cap program drop scanEJ
program define scanEJ, rclass
  args lf needle
  * Stata wraps a long display line in the log and continues it after "> "; rejoin
  * continuations before searching (the test_si7 lesson)
  tempname fh
  local found 0
  local acc
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if substr(`"`macval(line)'"', 1, 2) == "> " {
      local acc `"`macval(acc)'`=substr(`"`macval(line)'"', 3, .)'"'
    }
    else {
      if strpos(`"`macval(acc)'"', `"`needle'"') != 0 local found 1
      local acc `"`macval(line)'"'
    }
    file read `fh' line
  }
  if strpos(`"`macval(acc)'"', `"`needle'"') != 0 local found 1
  file close `fh'
  return scalar found = `found'
end

cap program drop mkbs37
program define mkbs37
  clear
  set obs 12
  qui gen i = _n
  qui gen ttre = cond(i <= 6, 4, .)
  qui expand 6
  qui bysort i: gen t = _n
  qui gen y = 1.5*(t >= ttre & ttre < .) + 0.3*i + 0.2*t
  qui gen flag = 1
  qui replace flag = 2 if inlist(i, 11, 12)
end

local needle37 nonzero but not exactly 1 on 12 observation(s)

mkbs37
qui log using edges_e14a.txt, text replace name(ed)
cap noi wooldid y i t ttre if flag, att
local rc14a = _rc
qui log close ed
local n14a = e(N)
local d14a = e(N_dropped)
scanEJ edges_e14a.txt "`needle37'"
local c14a = r(found)
cap erase edges_e14a.txt

mkbs37
qui log using edges_e14b.txt, text replace name(ed)
cap noi wooldid y i t ttre if flag != 0, att
local rc14b = _rc
qui log close ed
local n14b = e(N)
local d14b = e(N_dropped)
scanEJ edges_e14b.txt "`needle37'"
local c14b = r(found)
cap erase edges_e14b.txt

mkbs37
qui log using edges_e14c.txt, text replace name(ed)
cap noi wooldid y i t ttre if flag == 1, att
local rc14c = _rc
qui log close ed
local n14c = e(N)
local d14c = e(N_dropped)
scanEJ edges_e14c.txt "`needle37'"
local c14c = r(found)
cap erase edges_e14c.txt

mkbs37
qui replace flag = . if i == 10
qui log using edges_e14d.txt, text replace name(ed)
cap noi wooldid y i t ttre if flag, att
local rc14d = _rc
qui log close ed
local n14d = e(N)
local d14d = e(N_dropped)
scanEJ edges_e14d.txt "`needle37'"
local c14d = r(found)
cap erase edges_e14d.txt

di "E14 BS-37: bare rc=`rc14a' caution=`c14a' N=`n14a' drop=`d14a' | !=0 rc=`rc14b' caution=`c14b' N=`n14b' drop=`d14b' | ==1 rc=`rc14c' caution=`c14c' N=`n14c' drop=`d14c' | missing rc=`rc14d' caution=`c14d' N=`n14d' drop=`d14d'"
if `rc14a' != 0 | `c14a' != 1 | `n14a' != 60 | `d14a' != 12 {
  di as error "E14a FAIL: the bare-flag run does not name and count the rows the exact-one rule removed"
  local ++fails
}
if `rc14b' != 0 | `c14b' != 0 | `n14b' != 72 | `d14b' != 0 {
  di as error "E14b FAIL: a comparison if() must draw no caution and keep every row"
  local ++fails
}
if `rc14c' != 0 | `c14c' != 0 | `n14c' != 60 | `d14c' != 0 {
  di as error "E14c FAIL: an ordinary == comparison must contribute nothing to the fold"
  local ++fails
}
if `rc14d' != 0 | `c14d' != 1 | `n14d' != 54 | `d14d' != 12 {
  di as error "E14d FAIL: missing-valued rows must leave the sample and stay out of the count"
  local ++fails
}

* ---- E15 (change audit 2026-09-22, markout-fv-prefix): a factor-variable term whose operator
* the markout's three strips do not know (ib2.g) is marked out on its base variable ----
use "data/dta_binary", clear
qui gen g = 1 + mod(i, 3)
qui replace g = . in 1/40
local e15ok 1
foreach e15opt in controls fe {
  cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) `e15opt'(i.g)
  local rc15i = _rc
  local n15i = e(N)
  local a15i = e(b)[1,1]
  cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) `e15opt'(ib2.g)
  local rc15b = _rc
  local n15b = e(N)
  local a15b = e(b)[1,1]
  if `rc15b' != 0 local a15b = .
  di "E15 `e15opt': i.g rc=`rc15i' N=`n15i' att=`a15i' | ib2.g rc=`rc15b' N=`n15b' att=`a15b'"
  if `rc15i' != 0 | `rc15b' != 0 | `n15i' != 920 | `n15b' != 920 | reldif(`a15i', `a15b') > 1e-10 {
    di as error "E15 FAIL: `e15opt'(ib2.g) must run on the same sample and give the same att as `e15opt'(i.g)"
    local ++fails
    local e15ok 0
  }
}
if `e15ok' di as result "E15 PASS: ib2.g is marked out on g in controls() and fe()"

* ---- E16 (change audit round 2, S1 hyg-ereturn-internal-fit / U12 guide-failed-call-e-claim;
* the author's ruling 2026-09-23, "keep the promise"): a call that fails after the stage-1 fit
* leaves the caller's e() in place. A region-by-period fixed effect whose region 1 is exactly
* the 2005 cohort absorbs that cohort's cell dummies, and the post-fit collinearity screen
* refuses with r(322). e() used to hold wooldid's internal reghdfe fit on tempvar regressors
* (e(cmd) reghdfe); it must hold the regress run before it, e(sample) included
use data/dta_binary.dta, clear
gen byte region = 1 if ttre == 2005
replace region = 2 + mod(i, 2) if missing(region)
egen long rt = group(region t)
qui regress y t if i <= 40
local e16n = e(N)
local e16cl `"`e(cmdline)'"'
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fe(rt)
local rc16 = _rc
qui count if e(sample)
local e16s = r(N)
di "E16 rc=`rc16' e(cmd)=[`e(cmd)'] e(N)=" e(N) " e(sample)=`e16s' (want regress, `e16n', `e16n')"
if `rc16' != 322 | "`e(cmd)'" != "regress" | `"`e(cmdline)'"' != `"`e16cl'"' | e(N) != `e16n' | `e16s' != `e16n' {
  di as error "E16 FAIL: a call failing after the stage-1 fit must leave the caller's e() (regress) in place"
  local ++fails
}
else di as result "E16 PASS: the failed call left the caller's regress results in e()"

if `fails' == 0 di as result "EDGES_ALL_PASS"
else di as error "EDGES_FAILS=`fails'"
