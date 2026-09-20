* =================================================================
* probe_rng.do -- a probe NOBODY suggested: does the documented
* e(boot_rngstate) replay actually replay?
* sthlp (rseed()): "Whether or not rseed() is given, the random-number
* state in force before the first draw is stored in e(boot_rngstate) ...
* so a run without rseed() can still be replayed exactly by set rngstate
* to that value before repeating the call".
* Also: same rseed() twice must give identical draws; and the bootsave()
* characteristics _dta[wooldid_engine]/[wooldid_seed]/[wooldid_rngstate].
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_rng.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 111222
set sortseed 111222

set obs 30
gen long unit = _n
gen int ttre = .
replace ttre = 2004 if unit <= 10
replace ttre = 2006 if unit > 10 & unit <= 20
gen double ufe = rnormal()
expand 8
bysort unit: gen int year = 2000 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.5, 0))
gen double y  = 1 + ufe + 0.1*(year-2000) + te + rnormal()*0.7
save "sk_rng.dta", replace

* ---- A: run with NO rseed(); capture e(boot_rngstate) ----
use "sk_rng.dta", clear
set seed 998877
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(60) dots(0)
dis "RNG_A_RC=" _rc
local st0 = e(boot_rngstate)
matrix A = e(bootdist)
dis "RNG_A_SE=" %21x sqrt(el(e(V),1,1))
dis "RNG_A_STATELEN=" length("`st0'")
dis "RNG_A_SEED=" e(boot_seed)

* ---- B: set rngstate to that value and repeat the SAME call ----
use "sk_rng.dta", clear
set rngstate `st0'
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(60) dots(0)
dis "RNG_B_RC=" _rc
matrix B = e(bootdist)
dis "RNG_B_SE=" %21x sqrt(el(e(V),1,1))
dis "RNG_REPLAY_MRELDIF=" mreldif(A, B)
forvalues r = 1/4 {
  dis "RNG_REP`r'  A=" %12.9f el(A,`r',1) "   B=" %12.9f el(B,`r',1)
}

* ---- C: control -- same rseed() twice must be identical ----
use "sk_rng.dta", clear
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(60) rseed(31415) dots(0)
matrix C1 = e(bootdist)
use "sk_rng.dta", clear
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(60) rseed(31415) dots(0)
matrix C2 = e(bootdist)
dis "RNG_SEEDREPEAT_MRELDIF=" mreldif(C1, C2)

* ---- D: bootsave() characteristics ----
use "sk_rng.dta", clear
cap erase "skrng.dta"
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(30) rseed(31415) bootsave(skrng) dots(0)
use "skrng.dta", clear
dis "RNG_CHAR_ENGINE=" `"`: char _dta[wooldid_engine]'"'
dis "RNG_CHAR_SEED="   `"`: char _dta[wooldid_seed]'"'
local cs : char _dta[wooldid_rngstate]
dis "RNG_CHAR_STATELEN=" length(`"`cs'"')

dis "PROBE_RNG_DONE"
