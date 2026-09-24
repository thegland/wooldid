* test_bootstrap_std.do - commit-D engine validation: standard (resampling) mode.
* Run from test/: stata-mp -e do test_bootstrap_std.do
* DO1: rep-by-rep matching proof vs the OFFICIAL bootstrap prefix, same seed (count
* path; weight-injection == physical resampling ULP class). DO2: count path == re-id
* path for att (same-seed coherence across the split-engine legs). DO3: att_it under the
* re-id leg yields a usable se, and bootstrap(bayes) with att_it is refused under safety(off).
* DO4 failed-rep counting; DO5 observation-level; DO6 weighted + points never move;
* DO7 subgroups under the re-id leg; DO8 latticeignoreweights bins each replication's
* resample on the count and bayes legs. Prints BOOTSTD_ALL_PASS iff all pass.
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

* ---- D-O3: att_it under re-id runs and yields a usable se; bootstrap(bayes) with att_it,
* whose reweighting cannot resample the cells, is refused, safety(off) included
mkpan
cap qui wooldid y unit t ttre, att_it cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(12) rseed(777) safety(off)
local do3b = _rc
mkpan
qui wooldid y unit t ttre, att_it cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(12) rseed(777)
scalar itse_s = el(e(wdidmainresults),1,2)
di "DO3_ATTIT_SE_standard_reid=" itse_s " bayes_safetyoff_rc=`do3b'"
if itse_s >= . | itse_s <= 0 | `do3b' != 198 {
  di as error "DO3 FAIL: re-id att_it se missing (" itse_s "), or bootstrap(bayes) with att_it under safety(off) was not refused (rc `do3b')"
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

* ---- D-O8 (change audit round 2, S4 liw-countleg-bins-not-resampled; the author's ruling
* 2026-09-23): latticeignoreweights ignores the SAMPLE weight, not the resample. The count leg
* and bayes carry a replicate as [aw = base weight x multiplicity], and the lattice under
* latticeignoreweights used to drop the whole weight, so each replication binned the distinct
* drawn rows (bayes: the original sample) rather than the resample. It now keeps the
* multiplicity. On unweighted data the option is then a no-op on both legs, draw for draw, as
* it is analytically; with a sample weight the count leg equals the re-id leg (which bins the
* physically expanded resample) as DO2's att column does
use "data/dta_cont", clear
local do8L contreat(intensity) contreatcontrols(basic) cluster(i) lattice(all 3)
qui wooldid y i t ttre, `do8L' bootstrap(standard) reps(12) rseed(4242) dots(0)
matrix D8CN = e(bootdist)
qui wooldid y i t ttre, `do8L' latticeignoreweights bootstrap(standard) reps(12) rseed(4242) dots(0)
matrix D8CL = e(bootdist)
qui wooldid y i t ttre, `do8L' bootstrap(bayes) reps(12) rseed(4242) dots(0)
matrix D8BN = e(bootdist)
qui wooldid y i t ttre, `do8L' latticeignoreweights bootstrap(bayes) reps(12) rseed(4242) dots(0)
matrix D8BL = e(bootdist)
qui wooldid y i t ttre [aw=w], `do8L' latticeignoreweights bootstrap(standard) reps(12) rseed(4242) dots(0)
matrix D8WC = e(bootdist)
qui wooldid y i t ttre [aw=w], `do8L' att att_it latticeignoreweights bootstrap(standard) reps(12) rseed(4242) dots(0)
matrix D8WR = e(bootdist)
local w8c = mreldif(D8CN, D8CL)
local w8b = mreldif(D8BN, D8BL)
local w8w 0
foreach nm in main_att main_ame {
  local kc = colnumb(D8WC, "`nm'")
  local kr = colnumb(D8WR, "`nm'")
  if `kc' >= . | `kr' >= . local w8w = 1
  else {
    forval r = 1/12 {
      local w8w = max(`w8w', reldif(D8WC[`r', `kc'], D8WR[`r', `kr']))
    }
  }
}
di "DO8_LIW count noLIW-vs-LIW=" `w8c' " bayes noLIW-vs-LIW=" `w8b' " weighted count-vs-reid=" `w8w'
if `w8c' != 0 | `w8b' != 0 | `w8w' > 1e-10 {
  di as error "DO8 FAIL: under latticeignoreweights a replication did not bin its own resample"
  local ++fails
}
matrix drop D8CN D8CL D8BN D8BL D8WC D8WR

if `fails' == 0 di as result "BOOTSTD_ALL_PASS"
else di as error "BOOTSTD_FAILS=`fails'"
