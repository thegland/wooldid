* Targeted suite for the option-domain screen batch (issue files VAL-G19, VAL-G21A,
* VAL-G21B, VAL-G21C, VAL-G21E): five screens appended at the end of wooldid_validate
* close silent-accept paths (orphan contreatwithin; unrecognized contreatelasticitytype
* collapsing to dydx; contreatpoly(0/neg) running as order 1; unrecognized safety()
* tokens behaving as on; esrelativeto(>-1)+esf silently clamped to -1 in the worker).
* Run from test/: "<stata>" -e do test_valdomains.do -> VALDOMAINS_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 13579
set sortseed 13579
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

* ---- D1: VAL-G19 contreatwithin without contreat -> r(198), formerly rc 0 with the
* binary-design estimability checks silently disabled ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) contreatwithin
if _rc != 198 {
  dis as error "D1 FAIL: orphan contreatwithin rc=" _rc " (expected 198)"
  local ++fails
}

* ---- D2: contreatwithin with contreat still runs ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatwithin cluster(i) esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "D2 FAIL: legitimate contreatwithin rc=" _rc
  local ++fails
}

* ---- D3: VAL-G21A unrecognized/wrong-case elasticity tokens -> r(198); explicit valid
* literal passes ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatelasticitytype(bogus) cluster(i) esprelength(2) espostlength(2)
if _rc != 198 {
  dis as error "D3 FAIL: contreatelasticitytype(bogus) rc=" _rc " (expected 198)"
  local ++fails
}
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatelasticitytype(EYEX) cluster(i) esprelength(2) espostlength(2)
if _rc != 198 {
  dis as error "D3b FAIL: contreatelasticitytype(EYEX) rc=" _rc " (expected 198; tokens are exact-case)"
  local ++fails
}
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatelasticitytype(dydx) cluster(i) esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "D3c FAIL: explicit contreatelasticitytype(dydx) rc=" _rc
  local ++fails
}

* ---- D4: VAL-G21B contreatpoly must be >= 1; valid higher order still runs ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatpoly(0) cluster(i) esprelength(2) espostlength(2)
if _rc != 198 {
  dis as error "D4 FAIL: contreatpoly(0) rc=" _rc " (expected 198)"
  local ++fails
}
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatpoly(-1) cluster(i) esprelength(2) espostlength(2)
if _rc != 198 {
  dis as error "D4b FAIL: contreatpoly(-1) rc=" _rc " (expected 198)"
  local ++fails
}
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatpoly(2) cluster(i) esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "D4c FAIL: contreatpoly(2) rc=" _rc
  local ++fails
}

* ---- D5: VAL-G21C safety() token domain {on, off}; typos and case variants error ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(of)
if _rc != 198 {
  dis as error "D5 FAIL: safety(of) rc=" _rc " (expected 198)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(OFF)
if _rc != 198 {
  dis as error "D5b FAIL: safety(OFF) rc=" _rc " (expected 198; tokens are exact-case)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(on)
if _rc != 0 {
  dis as error "D5c FAIL: explicit safety(on) rc=" _rc
  local ++fails
}
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "D5d FAIL: safety(off) rc=" _rc
  local ++fails
}

* ---- D5e: the screen reads the RAW user value (safetyinput): contreat+contreatwithin
* internally rewrites the working safety copy to off, which would mask a bogus token ----
use "`DATA'/dta_cont", clear
qui log using vald_d5e.txt, text replace name(vald)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) contreatwithin safety(bogus) cluster(i) esprelength(2) espostlength(2)
local rc5e = _rc
qui log close vald
if `rc5e' != 198 {
  dis as error "D5e FAIL: contreatwithin+safety(bogus) rc=" `rc5e' " (expected 198)"
  local ++fails
}
else if strpos(fileread("vald_d5e.txt"), "safety() must be") == 0 {
  dis as error "D5e FAIL: expected safety() domain error text missing"
  local ++fails
}

* ---- D6: VAL-G21E esrelativeto() >= 0 under esfixedbaseperiod -> r(198), formerly a
* silent clamp to -1 in the worker; negative values and the no-esf caution unchanged ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) esrelativeto(0)
if _rc != 198 {
  dis as error "D6 FAIL: esf+esrelativeto(0) rc=" _rc " (expected 198)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) esrelativeto(2)
if _rc != 198 {
  dis as error "D6b FAIL: esf+esrelativeto(2) rc=" _rc " (expected 198)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) esrelativeto(-2)
if _rc != 0 {
  dis as error "D6c FAIL: esf+esrelativeto(-2) rc=" _rc
  local ++fails
}
use "`DATA'/dta_binary", clear
qui log using vald_d6d.txt, text replace name(vald)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) esrelativeto(2)
local rc6d = _rc
qui log close vald
if `rc6d' != 0 {
  dis as error "D6d FAIL: esrelativeto(2) without esf rc=" `rc6d' " (expected 0 + caution, unchanged)"
  local ++fails
}
else if strpos(fileread("vald_d6d.txt"), "no effect without esfixedbaseperiod") == 0 {
  dis as error "D6d FAIL: pre-existing inert-esrelativeto caution missing"
  local ++fails
}

* ---- D7: plain valid run untouched (screens are appended-only; valid inputs cannot
* reach them) ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
if _rc != 0 | e(b)[1,1] >= . {
  dis as error "D7 FAIL: plain valid run rc=" _rc
  local ++fails
}

* ---- D8 P6 (author ruling 2026-09-07): a change to the ESTIMAND WINDOW is disclosed whatever
* safety() says. subgroups() with a pre-treatment window and espostlength(0) silently raises
* espostlength to 1, so the run reports a post-treatment horizon the user did not ask for. The
* announcement used to be suppressed under safety(off). The fact of the change must now print
* either way; safety(off) still suppresses the advice beneath it. ----
use "`DATA'/dta_binary", clear
qui log using vald_d8.txt, text replace name(vd8)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(0) safety(off)
qui log close vd8
local d8fact 0
local d8advice 0
tempname fh8
file open `fh8' using vald_d8.txt, read text
file read `fh8' line
while r(eof) == 0 {
  if strpos(`"`line'"', "espostlength set to 1 due to use of subgroups") > 0 local d8fact 1
  if strpos(`"`line'"', "subgroups may crash if pre-treatment effects are not estimable") > 0 local d8advice 1
  file read `fh8' line
}
file close `fh8'
if `d8fact' != 1 | `d8advice' != 0 {
  dis as error "D8 FAIL: under safety(off) the espostlength raise must be announced (fact=`d8fact') and the advice suppressed (advice=`d8advice')"
  local ++fails
}
else dis as result "D8 PASS: safety(off) suppresses the advice, not the disclosure that the estimand window changed"

* ---- D9 (change audit 2026-09-22, val-i-assert-domain): the two i-domain screens run on the
* rows estimation can use. an i of 0, or a treated i past 2,147,483,647, confined to rows the
* markout drops (here: a missing outcome) is not refused, and the run equals the one that drops
* the same rows with if(); the same code on a kept row is still refused r(411). the screens
* used to read only the if() domain and refused the first two runs ----
use "`DATA'/dta_binary", clear
qui sum i if ttre < .
local d9tr = r(max)
cap noi wooldid y i t ttre if i != `d9tr', att cluster(i) esprelength(1) espostlength(1)
local rc9r = _rc
local n9r = e(N)
local a9r = e(b)[1,1]
qui gen double ymo = y
qui replace ymo = . if i == `d9tr'
qui gen double izero = i
qui replace izero = 0 if i == `d9tr'
qui gen double ibig = i
qui replace ibig = 3000000000 + i if i == `d9tr'
foreach d9v in izero ibig {
  cap noi wooldid ymo `d9v' t ttre, att cluster(`d9v') esprelength(1) espostlength(1)
  local rc9`d9v' = _rc
  local n9`d9v' = cond(_rc == 0, e(N), .)
  local a9`d9v' = cond(_rc == 0, e(b)[1,1], .)
  cap noi wooldid y `d9v' t ttre, att cluster(`d9v') esprelength(1) espostlength(1)
  local rk9`d9v' = _rc
}
dis "D9: reference rc=`rc9r' N=`n9r' att=`a9r' | i=0 on dropped rows rc=`rc9izero' N=`n9izero' att=`a9izero', on kept rows rc=`rk9izero' | big i on dropped rows rc=`rc9ibig' N=`n9ibig' att=`a9ibig', on kept rows rc=`rk9ibig'"
if `rc9r' != 0 | `rc9izero' != 0 | `rc9ibig' != 0 | `n9izero' != `n9r' | `n9ibig' != `n9r' | reldif(`a9izero', `a9r') > 1e-10 | reldif(`a9ibig', `a9r') > 1e-10 | `rk9izero' != 411 | `rk9ibig' != 411 {
  dis as error "D9 FAIL: the i-domain screens must read the estimation domain (markout folded in)"
  local ++fails
}
else dis as result "D9 PASS: an invalid i on markout-dropped rows runs; on a kept row it is refused"

if `fails' == 0 dis as result "VALDOMAINS_ALL_PASS"
else dis as error "VALDOMAINS_FAILS=`fails'"
