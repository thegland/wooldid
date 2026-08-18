* Targeted suite for the S6 cleanup-epilogue completion: the _ppmlhdfe_d obligation is
* registered at the creation moment (before est sto, whose failure is reachable with the
* estimates directory at Stata's 300-store limit), and the epilogue's drop is
* frame-qualified so it can never act in the wrong frame.
* Run from test/:  "<stata>" -e do test_s6.do   -> expect S6_ALL_PASS
* On the pre-fix build, S1 leaks _ppmlhdfe_d into the user's data.
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* ---- S1: est dir at the 300-store limit -> est sto fails mid-run; no variable leak ----
use "`DATA'/dta_binary", clear
qui reg y w
forval k = 1/300 {
  qui est sto s6fill`k'
}
cap wooldid ycount i t ttre, att cluster(i) poisson
local rc1 = _rc
if `rc1' == 0 {
  dis as error "S1 FAIL: expected the full estimates directory to abort the run"
  local ++fails
}
cap confirm variable _ppmlhdfe_d
if _rc == 0 {
  dis as error "S1 FAIL: _ppmlhdfe_d leaked into the data through the unregistered window"
  local ++fails
}
* the user's own stores must have survived the failed run
cap est restore s6fill300
if _rc != 0 {
  dis as error "S1 FAIL: user store s6fill300 destroyed (rc=" _rc ")"
  local ++fails
}
est drop _all

* ---- S2: recovery - the next clean poisson run works and cleans up after itself ----
use "`DATA'/dta_binary", clear
cap noi wooldid ycount i t ttre, att cluster(i) poisson
if _rc != 0 {
  dis as error "S2 FAIL: clean poisson run rc=" _rc
  local ++fails
}
cap confirm variable _ppmlhdfe_d
if _rc == 0 {
  dis as error "S2 FAIL: _ppmlhdfe_d left after a clean run"
  local ++fails
}

* ---- S3: verbose contract - _ppmlhdfe_d and the full model are retained on purpose ----
use "`DATA'/dta_binary", clear
cap noi wooldid ycount i t ttre, att cluster(i) poisson verbose
if _rc != 0 {
  dis as error "S3 FAIL: verbose poisson run rc=" _rc
  local ++fails
}
else {
  cap confirm variable _ppmlhdfe_d
  if _rc != 0 {
    dis as error "S3 FAIL: verbose no longer retains _ppmlhdfe_d"
    local ++fails
  }
  cap est restore ___wooldidfullmodel
  if _rc != 0 {
    dis as error "S3 FAIL: verbose no longer retains ___wooldidfullmodel"
    local ++fails
  }
  cap est drop ___wooldidfullmodel
  cap drop _ppmlhdfe_d
}

* ---- S4: failing run from a user frame - drop acts in that frame, decoys elsewhere survive ----
frame create s6deco
frame s6deco {
  qui set obs 5
  qui gen _ppmlhdfe_d = 42
}
use "`DATA'/dta_binary", clear
frame put *, into(s6work)
frame change s6work
qui reg y w
forval k = 1/300 {
  qui est sto s6fill`k'
}
cap wooldid ycount i t ttre, att cluster(i) poisson
qui frame
if "`r(currentframe)'" != "s6work" {
  dis as error "S4 FAIL: not restored to the calling frame (in `r(currentframe)')"
  local ++fails
}
cap confirm variable _ppmlhdfe_d
if _rc == 0 {
  dis as error "S4 FAIL: _ppmlhdfe_d leaked in the calling frame"
  local ++fails
}
est drop _all
frame s6deco {
  qui count if _ppmlhdfe_d == 42
  if r(N) != 5 {
    dis as error "S4 FAIL: decoy _ppmlhdfe_d in another frame was modified"
    local ++fails
  }
}
frame change default
frame drop s6work
frame drop s6deco

* ---- S5: cleanmatrices vs a validation-rejected call - user matrices and mata survive ----
* The S5m fix moved all parse-time validation (wooldid_validate) BEFORE the
* cleanmatrices wipe, so a call refused at validation must leave user matrices and
* mata objects untouched. fastnose+jointtests is a documented parse-time refusal (198).
* A clean call with cleanmatrices must still succeed (its wipe is the documented one).
use "`DATA'/dta_binary", clear
matrix s5umat = (1,2\3,4)
mata: s5uobj = (5,6\7,8)
cap wooldid y i t ttre, att cluster(i) fastnose jointtests cleanmatrices
if _rc == 0 {
  dis as error "S5 FAIL: fastnose+jointtests+cleanmatrices was not refused"
  local ++fails
}
cap confirm matrix s5umat
if _rc != 0 {
  dis as error "S5 FAIL: user matrix wiped by the rejected cleanmatrices call (rc=" _rc ")"
  local ++fails
}
cap mata: assert(s5uobj[2,2] == 8)
if _rc != 0 {
  dis as error "S5 FAIL: user mata object wiped by the rejected cleanmatrices call (rc=" _rc ")"
  local ++fails
}
* recovery: a clean run with cleanmatrices succeeds (note: it wipes matrices/mata by
* design, so this is the last block of the suite)
cap noi wooldid y i t ttre, att cluster(i) cleanmatrices
if _rc != 0 {
  dis as error "S5 FAIL: clean cleanmatrices run rc=" _rc
  local ++fails
}

if `fails' == 0 dis "S6_ALL_PASS"
else dis as error "S6_FAILS=`fails'"
exit, clear
