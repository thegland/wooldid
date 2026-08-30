* test_prefixharden.do - commit-A hardening: worker coeftabresults defense +
* shell properties(svyr) display routing. Run from test/: stata-mp -e do test_prefixharden.do
* Gates: T1 plain run + replay unchanged; T2 coeftabresults-off cure (parent r503);
* T3 naked bootstrap prefix end-to-end; T4 failure-marker semantics unchanged;
* T5 svy remains a loud failure. Prints PREFIXHARDEN_ALL_PASS on column 1 iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"

local fails 0

* build a small deterministic staggered panel
clear
set obs 8
gen i = _n
gen ttre = cond(i<=2, 4, cond(i<=4, 6, .))
expand 6
bysort i: gen t = _n
set seed 7
gen y = 2*(t>=ttre & ttre<.) + 0.1*i + 0.05*t + rnormal()
tempfile panel
save `panel'

* T1: plain run + replay both rc 0
cap noi qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
if _rc != 0 {
  di as error "T1a FAIL: plain run rc " _rc
  local fails = `fails' + 1
}
scalar t1att = _b[main_att]
cap qui wooldid
if _rc != 0 {
  di as error "T1b FAIL: replay rc " _rc
  local fails = `fails' + 1
}

* T2: the cure - wooldid must survive coeftabresults off (parent build: r503)
set coeftabresults off
use `panel', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
local t2rc = _rc
set coeftabresults on
if `t2rc' != 0 {
  di as error "T2 FAIL: run under coeftabresults off rc " `t2rc'
  local fails = `fails' + 1
}
else if reldif(_b[main_att], t1att) > 1e-14 {
  di as error "T2 FAIL: coeftabresults-off att differs: " _b[main_att] " vs " t1att
  local fails = `fails' + 1
}

* T3: naked official bootstrap prefix end-to-end (unweighted; idcluster form)
use `panel', clear
cap noi qui bootstrap _b, reps(8) seed(3) cluster(i) idcluster(ni) nowarn: wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) fastnose
if _rc != 0 {
  di as error "T3a FAIL: naked prefix rc " _rc
  local fails = `fails' + 1
}
else {
  if e(N_reps) != 8 {
    di as error "T3b FAIL: N_reps " e(N_reps) " != 8"
    local fails = `fails' + 1
  }
  cap qui estat bootstrap, all
  if _rc != 0 {
    di as error "T3c FAIL: estat bootstrap rc " _rc
    local fails = `fails' + 1
  }
  * bootstrap V posted and nondegenerate
  cap noi qui {
    matrix t3v = e(V)
    if t3v[1,1] >= . | t3v[1,1] <= 0 error 9
  }
  if _rc != 0 {
    di as error "T3d FAIL: bootstrap e(V) missing/degenerate"
    local fails = `fails' + 1
  }
}

* T4: failure-marker semantics unchanged - a failing call sets the marker and
* blocks replay; the next successful call clears it
use `panel', clear
cap qui wooldid y i t, att
if "${WOOLDID_LASTRC}" == "" {
  di as error "T4a FAIL: failed call did not set WOOLDID_LASTRC"
  local fails = `fails' + 1
}
cap qui wooldid
if _rc == 0 {
  di as error "T4b FAIL: replay after failure did not refuse"
  local fails = `fails' + 1
}
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
if _rc != 0 | "${WOOLDID_LASTRC}" != "" {
  di as error "T4c FAIL: successful call did not clear the marker (rc " _rc ")"
  local fails = `fails' + 1
}

* T5: svy stays a LOUD failure (properties(svyr) must not open a silent path)
use `panel', clear
cap qui svyset i
cap qui svy: wooldid y i t ttre, att esprelength(1) espostlength(1)
if _rc == 0 {
  di as error "T5 FAIL: svy: wooldid ran rc 0 (must fail loudly)"
  local fails = `fails' + 1
}

if `fails' == 0 {
  di as result "PREFIXHARDEN_ALL_PASS"
}
else {
  di as error "PREFIXHARDEN_FAILS=`fails'"
}
