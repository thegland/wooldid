* test_prefixharden.do - commit-A hardening: worker coeftabresults defense +
* shell properties(svyr) display routing. Run from test/: stata-mp -e do test_prefixharden.do
* Gates: T1 plain run + replay unchanged; T2 coeftabresults-off cure (parent r503);
* T3 naked bootstrap prefix end-to-end; T4 failure-marker semantics (a complete store replays
* after a failed call);
* T5 svy remains a loud failure; T6 a bare wooldid after the bootstrap or jackknife prefix
* prints the prefix's own table. Prints PREFIXHARDEN_ALL_PASS on column 1 iff all pass.
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

* T4: the failure marker. A failing call still sets it and the next successful call clears it
* (T4a, T4c). Since the author's round-2 rulings (2026-09-23: a failed call keeps the caller's
* e(); "replay a restored store") it no longer blocks the replay of a complete store: the bare
* replay after a failed call shows the earlier results (T4b, which used to require a refusal),
* an option typed on the replay line, itself a failed call, does not block the next bare replay
* (T4d), and a store selected with estimates restore after a failure replays (T4e)
use `panel', clear
qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
est sto t4good
cap qui wooldid y i t, att
if "${WOOLDID_LASTRC}" == "" {
  di as error "T4a FAIL: failed call did not set WOOLDID_LASTRC"
  local fails = `fails' + 1
}
qui log using prefixharden_t4.txt, text replace name(ph4)
cap noi wooldid
local t4rc = _rc
qui log close ph4
if `t4rc' != 0 | strpos(fileread("prefixharden_t4.txt"), "Wooldid Estimation") == 0 {
  di as error "T4b FAIL: the replay after a failed call did not show the earlier results (rc `t4rc')"
  local fails = `fails' + 1
}
cap qui wooldid, level(90)
local t4drc = _rc
cap qui wooldid
if `t4drc' == 0 | _rc != 0 {
  di as error "T4d FAIL: an option typed on the replay line (rc `t4drc') blocked the next bare replay (rc " _rc ")"
  local fails = `fails' + 1
}
cap qui wooldid y i t, att
qui est restore t4good
cap qui wooldid
if _rc != 0 {
  di as error "T4e FAIL: a store selected with estimates restore after a failed call did not replay (rc " _rc ")"
  local fails = `fails' + 1
}
cap erase prefixharden_t4.txt
cap est drop t4good
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

* T6 (change audit round 2, S1 prefix-replay-analytic; the author's ruling 2026-09-23): after
* Stata's bootstrap or jackknife prefix, e() is the prefix's post (e(prefix) set, e(V) the
* prefix's covariance), and a bare wooldid hands the display to Stata's _prefix_display, as a
* bare regress does: the prefix's own table, not wooldid's analytic tables with their analytic
* standard errors. A plain run's replay keeps wooldid's tables. The logs hold only the bare
* replay line, so neither needle can come from an echoed command
clear
set obs 24
gen i = _n
gen ttre = cond(i<=6, 4, cond(i<=12, 6, .))
expand 6
bysort i: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*i + 0.05*t + rnormal()
tempfile panel24
save `panel24'
local t6 0
foreach pfx in bootstrap jackknife {
  use `panel24', clear
  if "`pfx'" == "bootstrap" cap qui bootstrap _b, reps(8) seed(3) cluster(i) nowarn: wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
  else cap qui jackknife _b, cluster(i): wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
  local t6rc = _rc
  local t6pfx `e(prefix)'
  qui log using prefixharden_t6.txt, text replace name(ph6)
  cap noi wooldid
  local t6rc2 = _rc
  qui log close ph6
  local t6lg = fileread("prefixharden_t6.txt")
  local t6want = cond("`pfx'" == "bootstrap", "Bootstrap results", "Jackknife results")
  di "T6 `pfx': run rc=`t6rc' e(prefix)=[`t6pfx'] replay rc=`t6rc2' prefix table=" (strpos(`"`t6lg'"', "`t6want'") > 0) " wooldid banner=" (strpos(`"`t6lg'"', "Wooldid Estimation") > 0)
  if `t6rc' != 0 | "`t6pfx'" != "`pfx'" | `t6rc2' != 0 | strpos(`"`t6lg'"', "`t6want'") == 0 | strpos(`"`t6lg'"', "Wooldid Estimation") > 0 {
    di as error "T6 FAIL: a bare wooldid after the `pfx' prefix does not print the prefix's own table"
    local ++t6
  }
}
use `panel24', clear
qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
qui log using prefixharden_t6.txt, text replace name(ph6)
cap noi wooldid
local t6rc2 = _rc
qui log close ph6
local t6lg = fileread("prefixharden_t6.txt")
if `t6rc2' != 0 | strpos(`"`t6lg'"', "Wooldid Estimation") == 0 | strpos(`"`t6lg'"', "Bootstrap results") > 0 | strpos(`"`t6lg'"', "Jackknife results") > 0 {
  di as error "T6c FAIL: a plain run's replay no longer prints wooldid's own tables"
  local ++t6
}
cap erase prefixharden_t6.txt
local fails = `fails' + `t6'

if `fails' == 0 {
  di as result "PREFIXHARDEN_ALL_PASS"
}
else {
  di as error "PREFIXHARDEN_FAILS=`fails'"
}
