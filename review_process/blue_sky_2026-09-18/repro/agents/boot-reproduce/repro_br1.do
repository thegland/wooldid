* repro_br1.do -- permsave() has no comma check: it writes a file named "<name>, replace.dta".
* bootsave() refuses the same input with r(198) (RF-6, 2026-09-13).  Self-contained.
clear all
set more off
set processors 1
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 1
set sortseed 1
* --- a tiny staggered panel, built here, no repo data ---
clear
set obs 40
gen long unit = _n
gen long ttre = cond(unit <= 12, 2003, cond(unit <= 24, 2005, .))
expand 8
bysort unit: gen int t = 1999 + _n
set seed 7
gen double y = 0.4*unit/40 + 0.1*t + 1.5*(t >= ttre) + rnormal()
cap erase "brtest.dta"
cap erase "brtest, replace.dta"
dis "--- bootsave(brtest, replace): expect r(198), the RF-6 comma message ---"
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(9) dots(0) bootsave(brtest, replace)
dis "BOOTSAVE_RC=" _rc
dis "--- permsave(brtest, replace): what happens? ---"
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(9) dots(0) permsave(brtest, replace)
dis "PERMSAVE_RC=" _rc
cap confirm file "brtest, replace.dta"
dis "FILE_NAMED_WITH_COMMA_EXISTS_RC=" _rc  // 0 means the badly named file was written
cap confirm file "brtest.dta"
dis "FILE_NAMED_CORRECTLY_EXISTS_RC=" _rc   // 601 means it was NOT written under the intended name
dis "REPRO_BR1_DONE"
