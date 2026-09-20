*=====================================================================
* p2_twins.do -- perm*/boot* twin asymmetries beyond the comma check.
*  E1 permsave under permtstat: does the _rep==0 se_ row survive?
*     (it is read as e(perm_ase_obs) INSIDE a new frame)
*  E2 bootsave under boottstat+jointtests: do jointW_ columns survive?
*     (read as e(boot_jointW) INSIDE a new frame)
*  E3 bootsave under boottstat+adversarial: advW column
*  E4 layout diff permsave vs bootsave (subgroups + t-stat legs)
*  E5 dots(-5): silent under permute(), r(198) under bootstrap()?
*  E6 permsave + permonsetcohorts: onset_ column count vs cluster count
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
gen long cohid = cond(missing(ttre), 9999, ttre)
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double a_i = rnormal()
gen double y   = a_i + 0.10*(t-2000) + 1.50*(t>=ttre & ttre<.) + rnormal()
keep unit t ttre y sg cohid
compress
tempfile base
save `base'

foreach f in e1.dta e2.dta e3.dta e4p.dta e4b.dta e6.dta {
  cap erase "`f'"
}

*---------------------------------------------------------------------
* E1  permsave under permtstat
*---------------------------------------------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permtstat permsave(e1)
dis "E1_RC=" _rc
cap noi dis "E1_HAS_ASE_OBS=" (e(perm_ase_obs) != "")
cap noi mat list e(perm_ase_obs)
cap noi use e1, clear
dis "E1_USE_RC=" _rc
cap noi ds
cap noi list _rep main_att se_main_att in 1/4
cap noi count if _rep==0 & !missing(se_main_att)
cap noi dis "E1_OBSROW_SE_NONMISS=" r(N)

*---------------------------------------------------------------------
* E2  bootsave under boottstat + jointtests
*---------------------------------------------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) jointtests bootstrap(standard) reps(20) boottstat bootsave(e2)
dis "E2_RC=" _rc
cap noi mat list e(boot_jointW)
cap noi mat list e(boot_jointWobs)
cap noi use e2, clear
dis "E2_USE_RC=" _rc
cap noi ds
cap noi ds jointW*
dis "E2_JOINTW_VARS=" r(varlist)
cap noi ds w_*
cap noi dis "E2_HAS_W=" ("`r(varlist)'" != "")

*---------------------------------------------------------------------
* E3  bootsave under boottstat + adversarial
*---------------------------------------------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) adversarial bootstrap(standard) reps(20) boottstat bootsave(e3)
dis "E3_RC=" _rc
cap noi use e3, clear
dis "E3_USE_RC=" _rc
cap noi ds advW*
dis "E3_ADVW_VARS=" r(varlist)

*---------------------------------------------------------------------
* E4  layout diff with subgroups()
*---------------------------------------------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) subgroups(sg) permute(shuffle) permreps(20) permtstat permsave(e4p)
dis "E4P_RC=" _rc
cap noi use e4p, clear
dis "E4P_USE_RC=" _rc
cap noi ds

use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) subgroups(sg) bootstrap(standard) reps(20) boottstat bootsave(e4b)
dis "E4B_RC=" _rc
cap noi use e4b, clear
dis "E4B_USE_RC=" _rc
cap noi ds

*---------------------------------------------------------------------
* E5  dots(-5)
*---------------------------------------------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) dots(-5)
dis "E5_BOOT_DOTSNEG_RC=" _rc
use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) dots(-5)
dis "E5_PERM_DOTSNEG_RC=" _rc

*---------------------------------------------------------------------
* E6  permsave + permonsetcohorts (i = onset grouping, cluster finer)
*---------------------------------------------------------------------
use `base', clear
cap noi wooldid y cohid t ttre, att cluster(unit) permute(shuffle) permreps(20) permonsetcohorts permsave(e6)
dis "E6_RC=" _rc
cap noi use e6, clear
dis "E6_USE_RC=" _rc
cap noi ds onset_*
local ov `r(varlist)'
dis "E6_NONSET=" `: word count `ov''
cap noi list _rep main_att onset_1 onset_2 onset_3 in 1/4

dis "P2_DONE"
