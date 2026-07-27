* Targeted suite for the B16 adversarial degenerate-fit guards: the closed-form
* parse-time guard (esf + altmethod + esrelativeto(-1) + a single pre horizon) and
* the count backstop at the reg sites (data gaps hollowing out the fit sample).
* Run from test/:  "<stata>" -e do test_advguard.do   -> expect ADVGUARD_ALL_PASS
* On the pre-guard build, T1 and T2 abort with r(2001) instead of skipping cleanly.
clear all
set more off
set maxvar 20000
set graphics off
set seed 4242
set sortseed 4242
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* ---- T1: closed-form cell -> clean skip with Caution, estimation completes ----
* the deterministic 1-row fit needs min_t == 1 (a single pre period before the earliest
* onset): validate widens esprelength by one under esf to reserve the reference slot, so
* on longer pre histories the fit always has a real coefficient besides the planted row
clear
set obs 30
gen i = _n
gen ttre = cond(_n <= 10, 2, .)
expand 8
bysort i: gen t = _n
gen te = 2*(t >= ttre & ttre < .)
gen y = te + 0.05*i + rnormal()
drop te
qui log using advg_t1.txt, text replace name(advg)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(2) espostlength(2) adversarial altmethod_adversarial
local rc1 = _rc
qui log close advg
if `rc1' != 0 {
  dis as error "T1 FAIL: rc=" `rc1' " (expected clean skip, rc 0)"
  local ++fails
}
else {
  if "`e(adversarial_pval_att)'" != "" {
    dis as error "T1 FAIL: adversarial_pval_att posted despite the guard"
    local ++fails
  }
  cap confirm matrix e(b)
  if _rc != 0 {
    dis as error "T1 FAIL: e(b) missing - estimation did not complete"
    local ++fails
  }
  if strpos(fileread("advg_t1.txt"), "Caution: adversarial pre-trend testing was cancelled") == 0 {
    dis as error "T1 FAIL: no cancellation Caution in the console output"
    local ++fails
  }
}

* ---- T1b: esrelativeto(-1) + esprelength(1) on a long pre history is VIABLE ----
* (validate's esf reference-slot widening leaves a 2-row fit; the guard must not fire)
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(2) adversarial altmethod_adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T1b FAIL: rc=" _rc " pval=[`e(adversarial_pval_att)'] (viable fit was skipped)"
  local ++fails
}

* ---- T2: gap-data backstop cell -> clean skip with Caution ----
* panel where the treated cohort is missing exactly relyr -2 (one calendar year absent);
* the reported pre horizons are then the data hole and the reference hole, so the
* altmethod fit sample is a single planted row: only the runtime count can catch this
clear
set obs 30
gen i = _n
gen ttre = cond(_n <= 10, 5, .)
expand 10
bysort i: gen t = _n
gen te = 2*(t >= ttre & ttre < .)
gen y = te + 0.05*i + rnormal()
drop te
drop if ttre == 5 & t == 3
tempfile hole2
qui save `hole2'
qui log using advg_t2.txt, text replace name(advg)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(2) adversarial altmethod_adversarial
local rc2 = _rc
qui log close advg
if `rc2' != 0 {
  dis as error "T2 FAIL: rc=" `rc2' " (expected clean skip, rc 0)"
  local ++fails
}
else {
  if "`e(adversarial_pval_att)'" != "" {
    dis as error "T2 FAIL: adversarial_pval_att posted despite the backstop"
    local ++fails
  }
  cap confirm matrix e(b)
  if _rc != 0 {
    dis as error "T2 FAIL: e(b) missing - estimation did not complete"
    local ++fails
  }
  if strpos(fileread("advg_t2.txt"), "Caution: adversarial pre-trend testing was cancelled") == 0 {
    dis as error "T2 FAIL: no cancellation Caution in the console output"
    local ++fails
  }
}

* ---- T2b: same gap panel, plain adversarial -> must still produce a p ----
use `hole2', clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(2) adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T2b FAIL: rc=" _rc " pval=[`e(adversarial_pval_att)'] (plain path must be untouched)"
  local ++fails
}

* ---- T3: viable altmethod cell -> guard must not overfire ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-3) esprelength(2) espostlength(2) adversarial altmethod_adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T3 FAIL: rc=" _rc " pval=[`e(adversarial_pval_att)'] (viable fit was skipped)"
  local ++fails
}

* ---- T4: plain adversarial in the closed-form cell -> untouched ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(2) adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T4 FAIL: rc=" _rc " pval=[`e(adversarial_pval_att)'] (plain adversarial must be untouched)"
  local ++fails
}

* ---- T5: subgroups adversarial normal path survives the SG-site edits ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(2) espostlength(2) subgroups(sg) adversarial altmethod_adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T5 FAIL: rc=" _rc " main pval=[`e(adversarial_pval_att)']"
  local ++fails
}
else {
  local sgpcount 0
  local sc : e(scalars)
  foreach s of local sc {
    if strpos("`s'", "adversarial_pval_att_s") == 1 local ++sgpcount
  }
  if `sgpcount' == 0 {
    dis as error "T5 FAIL: no per-subgroup adversarial pvals posted"
    local ++fails
  }
}

* ---- T6: no spurious cancellation Caution when adversarial was never requested ----
* (the finaladvtest clearing branch runs unconditionally; its Caution must not)
use "`DATA'/dta_binary", clear
qui log using advg_t6.txt, text replace name(advg)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(-999)
local rc6 = _rc
qui log close advg
if `rc6' != 0 {
  dis as error "T6 FAIL: rc=" `rc6'
  local ++fails
}
else if strpos(fileread("advg_t6.txt"), "adversarial pre-trend testing was cancelled") != 0 {
  dis as error "T6 FAIL: cancellation Caution printed on a run that never requested adversarial"
  local ++fails
}

* ---- T7: timetrends + esf: the two planted zero reference rows alone must not pass ----
* treated units observed only from the second reference onward -> no real pre coefficient;
* pre-tightening, the backstop counted the planted rows and fit a fabricated (0,0) line
clear
set obs 30
gen i = _n
gen ttre = cond(_n <= 10, 5, .)
expand 10
bysort i: gen t = _n
gen te = 2*(t >= ttre & ttre < .)
gen y = te + 0.05*i + 0.1*t*(ttre < .) + rnormal()
drop te
drop if ttre == 5 & t < 3
qui log using advg_t7.txt, text replace name(advg)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(2) timetrends adversarial altmethod_adversarial
local rc7 = _rc
qui log close advg
if `rc7' != 0 {
  dis as error "T7 FAIL: rc=" `rc7'
  local ++fails
}
else {
  cap confirm matrix e(adversarial_line_att)
  if _rc == 0 {
    dis as error "T7 FAIL: fabricated reference line posted from planted rows alone"
    local ++fails
  }
  if "`e(adversarial_pval_att)'" != "" {
    dis as error "T7 FAIL: adversarial_pval_att posted despite the planted-rows-only fit"
    local ++fails
  }
  if strpos(fileread("advg_t7.txt"), "Caution: adversarial pre-trend testing was cancelled") == 0 {
    dis as error "T7 FAIL: no cancellation Caution in the console output"
    local ++fails
  }
}

if `fails' == 0 dis "ADVGUARD_ALL_PASS"
else dis as error "ADVGUARD_FAILS=`fails'"
exit, clear
