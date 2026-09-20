*=====================================================================
* p3_replay.do
*  A. sthlp:542 promises: "a run without rseed() can still be replayed
*     exactly by set rngstate to that value before repeating the call".
*     Test that claim for BOTH engines, draw by draw, via the save files.
*  B. sthlp:532 promises e(boot_jointWobs)/e(boot_sgW)/e(boot_sgWobs)
*     alongside e(boot_jointW).  Check both joint-test families.
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
compress
tempfile base
save `base'

foreach f in r1.dta r2.dta b1.dta b2.dta s1.dta s2.dta {
  cap erase "`f'"
}

*=====================================================================
* A1  permute() with NO rseed(): record e(perm_rngstate), replay
*=====================================================================
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(r1)
dis "A1a_RC=" _rc
local pstate = e(perm_rngstate)
dis "A1_STATELEN=" length("`pstate'")
scalar pA = e(perm_p_att)
cap noi dis "A1_PATT=" %21x e(perm_p_att)

set rngstate `pstate'
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(r2)
dis "A1b_RC=" _rc
cap noi dis "A1_PATT2=" %21x e(perm_p_att)

use r1, clear
cap noi cf _all using r2
dis "A1_CF_RC=" _rc
* an explicit numeric comparison too
gen double a1 = main_att
merge 1:1 _rep using r2, keepusing(main_att) nogen keep(3) update replace
* (cf is the real test; this is belt and braces)

*=====================================================================
* A2  permute() WITH rseed(): trivially reproducible?  (control)
*=====================================================================
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(4242) permsave(s1)
dis "A2a_RC=" _rc
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(4242) permsave(s2)
dis "A2b_RC=" _rc
use s1, clear
cap noi cf _all using s2
dis "A2_CF_RC=" _rc

*=====================================================================
* A3  bootstrap() with NO rseed(): record e(boot_rngstate), replay
*=====================================================================
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(b1)
dis "A3a_RC=" _rc
local bstate = e(boot_rngstate)
dis "A3_STATELEN=" length("`bstate'")
cap noi dis "A3_SE=" %21x e(boot_se_att)

set rngstate `bstate'
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(b2)
dis "A3b_RC=" _rc
cap noi dis "A3_SE2=" %21x e(boot_se_att)

use b1, clear
cap noi cf _all using b2
dis "A3_CF_RC=" _rc

*=====================================================================
* B  the joint-test e() promise, both families
*=====================================================================
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) subgroups(sg) esprelength(2) espostlength(2) jointtests bootstrap(standard) reps(20) boottstat
dis "B1_RC=" _rc
dis "B1_FAMILY=" e(boot_jt_family)
foreach m in boot_jointW boot_jointWobs boot_sgW boot_sgWobs {
  cap confirm matrix e(`m')
  dis "B1_HAS_`m'=" (_rc==0)
}

use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) subgroups(sg) esprelength(2) espostlength(2) jointtests bootstrap(standard) reps(20)
dis "B2_RC=" _rc
dis "B2_FAMILY=" e(boot_jt_family)
foreach m in boot_jointW boot_jointWobs boot_sgW boot_sgWobs {
  cap confirm matrix e(`m')
  dis "B2_HAS_`m'=" (_rc==0)
}

dis "P3_DONE"
