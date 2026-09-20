*=====================================================================
* f2_jointwobs.do -- SELF-CONTAINED repro.
* sthlp:532 (bootstrap() entry) states, for the joint pre-trend and
* across-subgroup equality tests: "e(boot_jt_family) records which
* (studentized or fixedmetric), the replication statistics are in
* e(boot_jointW) and e(boot_sgW) with the observed statistics in
* e(boot_jointWobs) and e(boot_sgWobs)".
* Under boottstat (family = studentized) three of those four are never
* posted: ado:6946 posts boot_jointW alone, while ado:6871-6872 and
* ado:7152-7153 post all four on the fixedmetric branch.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 60
gen long unit = _n
gen byte arm  = mod(_n,3)
gen int ttre  = .
replace ttre  = 2004 if arm==1
replace ttre  = 2006 if arm==2
gen byte sg   = mod(_n,2)
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double a_i = rnormal()
gen double y   = a_i + 0.10*(t-2000) + 1.50*(t>=ttre & ttre<.) + rnormal()
keep unit t ttre y sg
tempfile base
save `base'

* --- studentized family (boottstat) ---------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) subgroups(sg) esprelength(2) espostlength(2) jointtests bootstrap(standard) reps(20) boottstat
dis "STUD_RC=" _rc
dis "STUD_FAMILY=" e(boot_jt_family)
foreach m in boot_jointW boot_jointWobs boot_sgW boot_sgWobs {
  cap confirm matrix e(`m')
  dis "STUD_HAS_`m'=" (_rc==0)
}

* --- fixed-metric family (no boottstat) -----------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) subgroups(sg) esprelength(2) espostlength(2) jointtests bootstrap(standard) reps(20)
dis "FIX_RC=" _rc
dis "FIX_FAMILY=" e(boot_jt_family)
foreach m in boot_jointW boot_jointWobs boot_sgW boot_sgWobs {
  cap confirm matrix e(`m')
  dis "FIX_HAS_`m'=" (_rc==0)
}
dis "F2_DONE"
