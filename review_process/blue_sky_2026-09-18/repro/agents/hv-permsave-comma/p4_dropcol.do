*=====================================================================
* p4_dropcol.do
*  C. D-10 (author ruling 2026-09-11, "the two failure modes are
*     announced rather than left to be discovered in the saved file")
*     shipped on permsave only.  Same data, same collision:
*     permsave PRINTS a Caution and counts the dropped column;
*     bootsave drops its w_ column SILENTLY (ado:14912 `continue').
*     Collision forced with cluster ids .1 and -1, both of which
*     strtoname() renders "w__1" / "onset__1".
*  D. what row 0 of jointW_att holds in the e2.dta file p2 wrote
*  E. permsave(draws.dta, replace) -- the comma bug's uglier twin
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
* cluster id: distinct per unit, but unit 1 -> .1 and unit 2 -> -1,
* which strtoname() maps to the SAME name
gen double cl = unit
replace cl = .1 if unit==1
replace cl = -1 if unit==2
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double a_i = rnormal()
gen double y   = a_i + 0.10*(t-2000) + 1.50*(t>=ttre & ttre<.) + rnormal()
keep unit t ttre y cl
compress
tempfile base
save `base'
qui levelsof cl, local(cls)
dis "NCLUS=" `: word count `cls''

foreach f in c1.dta c2.dta d1.dta {
  cap erase "`f'"
}

*--- C1  permsave on the colliding design -----------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(cl) permute(shuffle) permreps(20) permsave(c1)
dis "C1_RC=" _rc
cap noi use c1, clear
cap noi ds onset_*
dis "C1_NONSET=" `: word count `r(varlist)''

*--- C2  bootsave on the SAME colliding design ------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(cl) bootstrap(standard) reps(20) bootsave(c2)
dis "C2_RC=" _rc
cap noi use c2, clear
cap noi ds w_*
dis "C2_NW=" `: word count `r(varlist)''

*--- D  jointW_att row 0 in the file p2 wrote --------------------------
cap noi use e2, clear
dis "D_USE_RC=" _rc
cap noi list _rep jointW_att in 1/4
cap noi summ jointW_att if _rep==0
cap noi dis "D_OBSW=" %21x r(mean)

*--- E  permsave(draws.dta, replace) ----------------------------------
use `base', clear
cap noi wooldid y unit t ttre, att cluster(cl) permute(shuffle) permreps(20) permsave(d1.dta, replace)
dis "E_RC=" _rc
cap confirm file "d1.dta, replace.dta"
dis "E_UGLY_EXISTS_RC=" _rc
cap confirm file "d1.dta"
dis "E_PLAIN_EXISTS_RC=" _rc

dis "P4_DONE"
