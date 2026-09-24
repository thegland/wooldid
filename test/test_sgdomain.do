* Targeted suite for the subgroup sign-screen domain harmonization (VAL-G15A): the
* negative-sg check moved from wooldid_validate (raw if-domain, no markout knowledge)
* into the worker beside the integer screen, on `ifstatement & internal-ttre < .`.
* S7-S9 pin the partial-labelling screen in wooldid_sgdom.
* Run from test/: "<stata>" -e do test_sgdomain.do -> SGDOM_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 24680
set sortseed 24680
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop mksg
program define mksg
  use "data/dta_binary", clear
end

* ---- S1: false-positive cure - sg=-1 confined to one missing-y treated-period row:
* was r(411); now rc 0 and bit-identical to the run with that row deleted ----
tempname B1 B2
mksg
qui su i if ttre < . & t >= ttre, meanonly
local i0 = r(min)
qui su t if ttre < . & t >= ttre & i == `i0', meanonly
local t0 = r(min)
qui replace y = . if i == `i0' & t == `t0'
qui replace sg = -1 if i == `i0' & t == `t0'
qui count if sg == -1
local npoison = r(N)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 | `npoison' == 0 {
  dis as error "S1 FAIL: rc=" _rc " npoison=" `npoison' " (expected rc 0 on poisoned markout row)"
  local ++fails
}
else {
  matrix `B1' = e(b)
  mksg
  qui drop if i == `i0' & t == `t0'
  qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
  matrix `B2' = e(b)
  if mreldif(`B1', `B2') != 0 {
    dis as error "S1 FAIL: poisoned-markout-row twin not bit-identical, mreldif=" mreldif(`B1', `B2')
    local ++fails
  }
}

* ---- S2: sg=-1 on a live treated-post row still errors r(411) ----
mksg
qui replace sg = -1 if ttre < . & t >= ttre & _n < 50
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 411 {
  dis as error "S2 FAIL: live treated-post negative sg rc=" _rc " (expected 411)"
  local ++fails
}

* ---- S3: sg=-1 on a live PRE-onset row of a treated unit still errors r(411)
* (the SG pre-trend margins consume pre-onset rows; the t>=ttre domain would miss this) ----
mksg
qui replace sg = -1 if ttre < . & t < ttre & _n < 300
qui count if sg == -1
local npre = r(N)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 411 | `npre' == 0 {
  dis as error "S3 FAIL: pre-onset negative sg rc=" _rc " npre=" `npre' " (expected 411)"
  local ++fails
}

* ---- S4: newly-tolerated class certification - sg=-1 confined to never-treated rows
* runs rc 0 bit-identical to the same run with sg=. on those rows ----
tempname B4 B5
mksg
qui replace sg = -1 if ttre >= . & mod(i, 7) == 0
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
matrix `B4' = e(b)
mksg
qui replace sg = . if ttre >= . & mod(i, 7) == 0
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
matrix `B5' = e(b)
if mreldif(`B4', `B5') != 0 {
  dis as error "S4 FAIL: never-treated sg=-1 vs sg=. twins differ, mreldif=" mreldif(`B4', `B5')
  local ++fails
}

* ---- S5: doubly-bad value (sg = -1.5) on a live treated-post row: negative message
* fires first (sign screen precedes the integer screen) ----
mksg
qui log using sgdom_s5.txt, text replace name(sgd)
qui replace sg = -1.5 if ttre < . & t >= ttre & _n < 50
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
local rc5 = _rc
qui log close sgd
if `rc5' != 411 {
  dis as error "S5 FAIL: doubly-bad sg rc=" `rc5' " (expected 411)"
  local ++fails
}
else if strpos(fileread("sgdom_s5.txt"), "contains negative values") == 0 {
  dis as error "S5 FAIL: negative-message-first precedence lost for doubly-bad values"
  local ++fails
}

* ---- S6: clean sg run unchanged (screen is domain-relaxation only) ----
mksg
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 {
  dis as error "S6 FAIL: clean sg run rc=" _rc
  local ++fails
}

* ---- S7-S9: the partial-labelling screen. A treated unit labelled on some of its
* in-sample rows and missing on others that carry a treatment dummy is refused r(459)
* under default safety (S7) and proceeds with the caution under safety(off) (S8). Unit 1
* (onset 2005, sg 2) loses its label on 2005 and 2006, relative periods 0 and 1, so 2 rows.
* Passing shapes (S9), all under default safety: fully labelled, a wholly unlabelled
* treated unit, and labels missing on some rows of never-treated units only ----
cap program drop sglog
program define sglog, rclass
  * run wooldid ${SGCMD} on the data in memory under a named log; return the rc, and leave
  * the log text, batch-wrap continuations removed, in the string scalar SGLOG
  args logf
  qui log using `logf', text replace name(sgd)
  cap noi wooldid ${SGCMD}
  local rc = _rc
  qui log close sgd
  scalar SGLOG = subinstr(subinstr(fileread("`logf'"), char(13), "", .), char(10) + "> ", "", .)
  cap erase `logf'
  return scalar rc = `rc'
end
local perr "Error: the subgroup variable sg is missing for some in-sample row(s) of treated units of i that carry a subgroup value on their other rows"
local pcau "Caution: the subgroup variable sg is missing on 2 in-sample row(s) of treated units of i that carry a subgroup value on their other rows; proceeding under safety(off) using share weighting"
local pany "that carry a subgroup value on their other rows"

mksg
qui replace sg = . if i == 1 & inlist(t, 2005, 2006)
global SGCMD y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) subgroups(sg)
sglog sgdom_s7.txt
if r(rc) != 459 | strpos(SGLOG, "`perr'") == 0 {
  dis as error "S7 FAIL: partially labelled treated unit rc=" r(rc) " (expected 459 with the partial-labelling error)"
  local ++fails
}

mksg
qui replace sg = . if i == 1 & inlist(t, 2005, 2006)
global SGCMD y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) subgroups(sg) safety(off)
sglog sgdom_s8.txt
if r(rc) != 0 | strpos(SGLOG, "`pcau'") == 0 {
  dis as error "S8 FAIL: partially labelled treated unit under safety(off) rc=" r(rc) " (expected 0 with the partial-labelling caution on 2 rows)"
  local ++fails
}

foreach shape in full unlabelled nevertreated {
  mksg
  if "`shape'" == "unlabelled" qui replace sg = . if i == 1
  if "`shape'" == "nevertreated" qui replace sg = . if ttre >= . & mod(i, 7) == 0 & t >= 2008
  global SGCMD y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) subgroups(sg)
  sglog sgdom_s9.txt
  if r(rc) != 0 | strpos(SGLOG, "`pany'") > 0 {
    dis as error "S9 FAIL: passing shape `shape' rc=" r(rc) " (expected 0 with no partial-labelling message)"
    local ++fails
  }
}
global SGCMD
scalar drop SGLOG

* ---- S10: the within-cell and split-unit screens run once, on the observed run, not in every
* replication (change audit round 2, U11 sgdom-guards-per-replication). Like the onset and si8
* screens (f8b4a64) and the partial-labelling screen above, they now skip replication mode
* (`rm'): a replication cannot newly trip them (the bootstrap re-id leg makes i finer, and
* permute()'s own screen refuses every design unit whose subgroup varies within a cluster).
* wooldid_sgdom is called directly, as the worker calls it, on a split unit (unit 1 moves from
* its subgroup to another from 2007) and on a within-cell split (unit 1's post-onset rows
* doubled with the copy in another subgroup): r(459) on the observed run (S10a, S10c), rc 0
* in replication mode (S10b, S10d). S10e-f: end to end, both engines still refuse the
* split-unit data on the observed run
mksg
qui gen relt = cond(ttre < ., t - ttre, 0)
qui replace sg = cond(sg == 1, 2, 1) if i == 1 & t >= 2007
cap noi wooldid_sgdom , sg(sg) ivar(i) tvar(t) ttrevar(ttre) ifs(1) relt(relt)
local s10a = _rc
cap noi wooldid_sgdom , sg(sg) ivar(i) tvar(t) ttrevar(ttre) ifs(1) relt(relt) rm(1)
local s10b = _rc
mksg
qui gen relt = cond(ttre < ., t - ttre, 0)
qui expand 2 if i == 1 & t >= ttre, generate(s10dup)
qui replace sg = cond(sg == 1, 2, 1) if s10dup
cap noi wooldid_sgdom , sg(sg) ivar(i) tvar(t) ttrevar(ttre) ifs(1) relt(relt)
local s10c = _rc
cap noi wooldid_sgdom , sg(sg) ivar(i) tvar(t) ttrevar(ttre) ifs(1) relt(relt) rm(1)
local s10d = _rc
mksg
qui replace sg = cond(sg == 1, 2, 1) if i == 1 & t >= 2007
cap wooldid y i t ttre, att cluster(i) subgroups(sg) bootstrap(standard) reps(5) rseed(1) dots(0)
local s10e = _rc
cap wooldid y i t ttre, att cluster(i) subgroups(sg) permute(shuffle) permreps(20) rseed(1) dots(0)
local s10f = _rc
dis "S10 info: split unit observed=`s10a' replication=`s10b'; within-cell observed=`s10c' replication=`s10d'; bootstrap=`s10e' permute=`s10f'"
if `s10a' != 459 | `s10c' != 459 | `s10e' != 459 | `s10f' != 459 {
  dis as error "S10 FAIL: a subgroup screen no longer refuses on the observed run"
  local ++fails
}
if `s10b' != 0 | `s10d' != 0 {
  dis as error "S10 FAIL: the within-cell or split-unit screen still runs in replication mode"
  local ++fails
}

if `fails' == 0 dis as result "SGDOM_ALL_PASS"
else dis as error "SGDOM_FAILS=`fails'"
