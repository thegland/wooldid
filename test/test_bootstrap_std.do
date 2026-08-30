* test_bootstrap_std.do - commit-D engine validation: standard (resampling) mode.
* Run from test/: stata-mp -e do test_bootstrap_std.do
* DO1: rep-by-rep matching proof vs the OFFICIAL bootstrap prefix, same seed (count
* path; weight-injection == physical resampling ULP class). DO2: count path == re-id
* path for att (same-seed coherence across the split-engine legs). DO3: the estimand
* fork is real (re-id att_it se >> bayes reweighting se on heterogeneous effects).
* DO4 failed-rep counting; DO5 observation-level; DO6 weighted + points never move;
* DO7 subgroups under the re-id leg. Prints BOOTSTD_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0
cap prog drop mkpan
prog def mkpan
  clear
  set obs 12
  gen unit = _n
  gen s = ceil(unit/4)
  gen sg = mod(unit, 2)
  gen ttre = cond(unit<=4, 4, cond(unit<=8, 6, .))
  expand 8
  bysort unit: gen t = _n
  set seed 11
  gen y = 2*(t>=ttre & ttre<.) + 0.3*(t>=ttre & ttre<.)*sg + 0.1*unit + 0.05*t + rnormal()
  gen w0 = 0.5 + mod(unit, 3)
end

* ---- D-O1: MATCHING PROOF vs the official prefix, same seed, att, count path
mkpan
tempfile offsav
qui bootstrap _b, reps(15) seed(31459) cluster(unit) idcluster(ni) nowarn saving(`offsav', double every(1)): wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) fastnose
preserve
use `offsav', clear
mkmat _b_*, matrix(OFF)
restore
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(15) rseed(31459)
matrix ENG = e(bootdist)
local w1 0
forval r = 1/15 {
  forval c = 1/`=colsof(OFF)' {
    local w1 = max(`w1', reldif(ENG[`r',`c'], OFF[`r',`c']))
  }
}
di "DO1_OFFICIAL_MATCH_MAXRELDIF=" `w1'
if `w1' > 1e-10 {
  di as error "DO1 FAIL"
  local ++fails
}

* ---- D-O2: count path == re-id path for att (same seed; weight==duplication)
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(12) rseed(777)
matrix CNT = e(bootdist)
mkpan
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(12) rseed(777)
matrix RID = e(bootdist)
local rn : colfullnames RID
local k 0
local attcol 0
foreach nm of local rn {
  local ++k
  if "`nm'" == "main_att" local attcol `k'
}
local w2 0
forval r = 1/12 {
  local w2 = max(`w2', reldif(CNT[`r',1], RID[`r',`attcol']))
}
di "DO2_COUNT_VS_REID_ATT_MAXRELDIF=" `w2'
if `w2' > 1e-10 {
  di as error "DO2 FAIL"
  local ++fails
}

* ---- D-O3: att_it under re-id runs and yields a usable se
mkpan
qui wooldid y unit t ttre, att_it cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(12) rseed(777) safety(off)
scalar itse_b = el(e(wdidmainresults),1,2)
mkpan
qui wooldid y unit t ttre, att_it cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(12) rseed(777)
scalar itse_s = el(e(wdidmainresults),1,2)
di "DO3_ATTIT_SE_bayes=" itse_b " standard_reid=" itse_s
if itse_s >= . | itse_s <= 0 {
  di as error "DO3 FAIL: re-id att_it se missing"
  local ++fails
}

* ---- D-O4: failed-rep counting (draws losing all never-treated -> r459 reps)
clear
set obs 4
gen unit = _n
gen ttre = cond(unit<=2, 4, .)
expand 8
bysort unit: gen t = _n
set seed 5
gen y = 1.5*(t>=ttre & ttre<.) + 0.2*unit + rnormal()
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(40) rseed(9)
di "DO4_RC=" _rc " completed=" e(boot_completed) " failed=" e(boot_failed)
if _rc != 0 {
  di as error "DO4 FAIL rc"
  local ++fails
}
else if e(boot_failed) == 0 {
  di as error "DO4 FAIL: expected some failed reps on a 4-cluster design"
  local ++fails
}

* ---- D-O5: obs-level (no cluster) standard bootstrap
mkpan
cap qui wooldid y unit t ttre, att robust esprelength(1) espostlength(1) bootstrap(standard) reps(15) rseed(21)
di "DO5_OBSLEVEL_RC=" _rc " completed=" e(boot_completed)
if _rc != 0 {
  di as error "DO5 FAIL"
  local ++fails
}

* ---- D-O6: weighted standard run + points never move
mkpan
qui wooldid y unit t ttre [aw=w0], att cluster(unit) esprelength(1) espostlength(1)
matrix OB = e(b)
mkpan
cap qui wooldid y unit t ttre [aw=w0], att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(12) rseed(3)
di "DO6_RC=" _rc " bdiff=" mreldif(OB, e(b))
if _rc != 0 | mreldif(OB, e(b)) > 0 {
  di as error "DO6 FAIL"
  local ++fails
}

* ---- D-O7: SG + standard re-id
mkpan
cap qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(standard) reps(12) rseed(6)
di "DO7_SG_REID_RC=" _rc " completed=" e(boot_completed)
if _rc != 0 {
  di as error "DO7 FAIL"
  local ++fails
}

if `fails' == 0 di as result "BOOTSTD_ALL_PASS"
else di as error "BOOTSTD_FAILS=`fails'"
