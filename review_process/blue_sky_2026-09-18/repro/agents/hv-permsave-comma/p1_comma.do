*=====================================================================
* p1_comma.do -- SELF-CONTAINED repro: permsave() lacks bootsave()'s
*                comma check, so permsave(f, replace) writes a file
*                literally named "f, replace.dta" at rc 0.
* Run:  cd <this dir>
*       /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do p1_comma.do
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* --- a plain staggered panel: 60 units, 8 years, two treated cohorts -----
clear
set obs 60
gen long unit = _n
gen byte arm  = mod(_n,3)
gen int ttre  = .
replace ttre  = 2004 if arm==1
replace ttre  = 2006 if arm==2
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double a_i = rnormal()
gen double y   = a_i + 0.10*(t-2000) + 1.50*(t>=ttre & ttre<.) + rnormal()
keep unit t ttre y
compress

* wipe any leftovers from a previous run
foreach f in "draws, replace.dta" "draws.dta" "pdraws.dta" "bdraws.dta" "bd2, replace.dta" {
  cap erase "`f'"
}

dis "PWD=" c(pwd)

*---------------------------------------------------------------------
* T1  permsave(draws, replace)  -- the Stata-idiomatic call
*---------------------------------------------------------------------
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(draws, replace)
dis "T1_PERMSAVE_COMMA_RC=" _rc

cap confirm file "draws, replace.dta"
dis "T1_FILE_COMMA_EXISTS_RC=" _rc
cap confirm file "draws.dta"
dis "T1_FILE_PLAIN_EXISTS_RC=" _rc


*---------------------------------------------------------------------
* T1b  the SAME command a second time.  The pre-flight now finds
*      "draws, replace.dta" and (because `replace' was never parsed as
*      an option) tells the user to add the replace option they wrote.
*---------------------------------------------------------------------
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(draws, replace)
dis "T1b_PERMSAVE_COMMA_AGAIN_RC=" _rc

*---------------------------------------------------------------------
* T2  bootsave(bd2, replace) -- the twin.  Expect a helpful r(198).
*---------------------------------------------------------------------
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(bd2, replace)
dis "T2_BOOTSAVE_COMMA_RC=" _rc
cap confirm file "bd2, replace.dta"
dis "T2_FILE_COMMA_EXISTS_RC=" _rc

*---------------------------------------------------------------------
* T3  replace semantics, permsave: run once plain, then again plain.
*---------------------------------------------------------------------
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(pdraws)
dis "T3a_PERMSAVE_FIRST_RC=" _rc
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(pdraws)
dis "T3b_PERMSAVE_SECOND_NOREPLACE_RC=" _rc
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(pdraws) replace
dis "T3c_PERMSAVE_SECOND_REPLACE_RC=" _rc

*---------------------------------------------------------------------
* T4  replace semantics, bootsave (pre-flight, cheap)
*---------------------------------------------------------------------
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(bdraws)
dis "T4a_BOOTSAVE_FIRST_RC=" _rc
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(bdraws)
dis "T4b_BOOTSAVE_SECOND_NOREPLACE_RC=" _rc

*---------------------------------------------------------------------
* T5  what does the mis-named file contain?  is it a real draws file?
*---------------------------------------------------------------------
cap noi use "draws, replace.dta", clear
dis "T5_USE_RC=" _rc
cap noi dis "T5_N=" _N
cap noi ds
cap noi char list _dta[]

*---------------------------------------------------------------------
* T6  the "obvious next step" a user takes: use draws
*---------------------------------------------------------------------
cap noi use draws, clear
dis "T6_USE_PLAIN_RC=" _rc

dis "P1_DONE"
