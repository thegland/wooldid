* test_jointclus.do — certification of the SI-A1 joint-test screen (author-ruled
* 2026-08-25): jointtests + analytic clustered SEs clustered at the cohort-i level
* or above -> error 198 under default safety (the joint p-values are severely
* oversized there); safety(off) downgrades to a caution and proceeds with the
* p-values computed; finer clustering, robust, and jointtests-free runs untouched.
* customvce cluster-type counts as clustered; two-way cluster fires on ANY
* i-or-above member. Adversarial deliberately not gated (separate ruling).
* Run from test/:  "<stata>" -e do test_jointclus.do   -> expect JOINTCLUS_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 4411
set sortseed 4411

local REPO ".."
do "`REPO'/wooldidinstall/wooldid.ado"

global pf 0

cap program drop mkjc
program define mkjc
  clear
  set seed 4411
  qui set obs 12
  gen i = _n
  gen state = ceil(i/3)
  gen ttre = cond(_n <= 4, 5, cond(_n <= 8, 7, .))
  qui expand 4
  bysort i: gen unit = _n
  gen unitid = (i-1)*4 + unit
  qui expand 10
  bysort i unitid: gen t = _n
  gen y = 0.05*i + 0.04*t + 0.6*(ttre < . & t >= ttre) + rnormal()
  sort i unitid t
  ereturn clear
end

cap program drop scanlog
program define scanlog, rclass
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

cap program drop runjc
program define runjc
  args tag opts wantrc wanterr wantcau
  mkjc
  qui log using jc_`tag'.txt, text replace name(jc)
  cap noi wooldid y i t ttre, att esprelength(2) espostlength(2) `opts'
  local rc = _rc
  qui log close jc
  scanlog jc_`tag'.txt "Error: jointtests is not compatible with wooldid's analytic clustered SE"
  local fe = r(found)
  scanlog jc_`tag'.txt "Caution: jointtests was combined with analytic clustered SEs clustered at the level"
  local fc = r(found)
  cap erase jc_`tag'.txt
  if `rc' == `wantrc' & `fe' == `wanterr' & `fc' == `wantcau' dis as result "`tag' PASS"
  else {
    dis as error "`tag' FAIL: rc=`rc' (want `wantrc') err=`fe' (want `wanterr') cau=`fc' (want `wantcau')"
    global pf = $pf + 1
  }
end

* T1: cluster(i) + jointtests -> refused r198 with the error text
runjc T1 "cluster(i) jointtests" 198 1 0

* T2: safety(off) -> caution + proceed, joint p's actually computed
mkjc
qui log using jc_T2.txt, text replace name(jc)
cap noi wooldid y i t ttre, att esprelength(2) espostlength(2) cluster(i) jointtests safety(off)
local rc2 = _rc
qui log close jc
scanlog jc_T2.txt "Caution: jointtests was combined with analytic clustered SEs clustered at the level"
local fc2 = r(found)
cap erase jc_T2.txt
if `rc2' == 0 & `fc2' == 1 & "`e(joint_horiz)'" != "" dis as result "T2 PASS (joint_horiz = `e(joint_horiz)')"
else {
  dis as error "T2 FAIL: rc=`rc2' cau=`fc2' joint_horiz=[`e(joint_horiz)']"
  global pf = $pf + 1
}

* T3: robust + jointtests -> clean (no screen text)
runjc T3 "robust jointtests" 0 0 0

* T4: cluster(unitid) + jointtests -> clean (clustering FINER than i)
runjc T4 "cluster(unitid) jointtests" 0 0 0

* T5: cluster(state) + jointtests -> refused (clustering ABOVE i)
runjc T5 "cluster(state) jointtests" 198 1 0

* T6: two-way cluster(i unitid) + jointtests -> refused (ANY member i-or-above)
runjc T6 "cluster(i unitid) jointtests" 198 1 0

* T7: cluster(i) WITHOUT jointtests -> screen entirely dormant
runjc T7 "cluster(i)" 0 0 0

* T8: customvce cluster-type at the i level counts as clustered
runjc T8 "customvce(cluster i) jointtests" 198 1 0

if $pf == 0 dis as result "JOINTCLUS_ALL_PASS"
else dis as error "JOINTCLUS_FAILURES: $pf"
exit, clear
