* =================================================================
* probe_br1.do  -- SKEPTIC repro of BR-1: permsave() has no comma check
* Self-contained: generates its own data, needs nothing else.
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_br1.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242

which wooldid

* ---- data ----
set obs 40
gen long unit = _n
gen int  ttre = .
replace  ttre = 2004 if unit <= 13
replace  ttre = 2006 if unit > 13 & unit <= 26
gen double ufe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double te  = cond(missing(ttre), 0, cond(year >= ttre, 0.5, 0))
gen double y   = 1 + ufe + 0.12*(year-1999) + te + rnormal()*0.7

* ---- housekeeping: make sure no stale files ----
cap erase "skperm.dta"
cap erase "skperm, replace.dta"
cap erase "skboot.dta"
cap erase "skboot, replace.dta"

dis "===== LEG 1: bootsave(skboot, replace)  -- expect r(198) comma refusal ====="
cap noi wooldid y unit year ttre, att bootstrap(standard) reps(20) rseed(7) cluster(unit) bootsave(skboot, replace) dots(0)
dis "B1_BOOTCOMMA_RC=" _rc

dis "===== LEG 2: permsave(skperm, replace) -- claim says rc 0 and a comma filename ====="
cap noi wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(7) permsave(skperm, replace) dots(0)
dis "B1_PERMCOMMA_RC=" _rc

cap confirm file "skperm, replace.dta"
dis "B1_FILE_WITH_COMMA_RC=" _rc
cap confirm file "skperm.dta"
dis "B1_FILE_PLAIN_RC=" _rc

dis "===== LEG 3: control -- permsave(skperm2) plain name should work ====="
cap erase "skperm2.dta"
cap noi wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(7) permsave(skperm2) dots(0)
dis "B1_PERMPLAIN_RC=" _rc
cap confirm file "skperm2.dta"
dis "B1_FILE_PERM2_RC=" _rc

dis "===== LEG 4: does the comma file actually contain the draws? ====="
cap noi use "skperm, replace.dta", clear
dis "B1_USE_RC=" _rc
cap noi describe, short
cap noi dis "B1_NOBS=" _N

dis "PROBE_BR1_DONE"
