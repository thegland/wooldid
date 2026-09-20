*==========================================================================
* r11_remedy.do -- is there a user-side remedy for HV-BOOT-01/04/05?
* The bootstrap frame DOES honour the user `if' (smp = ((ifr)==1)), so an
* explicit `if y < .' should restore the correct cluster universe without
* the user deleting any rows.  Tested here rather than asserted.
* Also tested: the always-treated trigger (N1), where the analogous `if'
* is NOT obvious and wooldid's own advice (truncate the window) does not
* remove the phantom clusters.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- A: missing-outcome trigger --------------------------------------
clear
set seed 4242
qui set obs `=20*2*12'
qui gen int state = ceil(_n / 24)
qui gen int unit  = ceil(_n / 12)
qui bysort unit: gen int year = 1999 + _n
qui gen int onset = .
qui replace onset = 2003 + mod(state - 1, 4) if state <= 8
qui gen double a_i = 0
qui bysort unit: replace a_i = rnormal() if _n == 1
qui bysort unit: replace a_i = a_i[1]
qui gen byte post = (onset < . & year >= onset)
qui gen double y = 1 + a_i + .6*post + rnormal()*.3
qui replace y = . if state <= 4
tempfile A
qui save `A'
dis "--- A1. no if ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(80) rseed(9) dots(0)
dis "A1_RC=" _rc
dis "RESULT A1 N=" e(N)
dis "RESULT A1 bootdf=" e(boot_df)
dis "--- A2. if y < . ---"
use `A', clear
cap noi wooldid y unit year onset if y < . , cluster(state) bootstrap(standard) reps(80) rseed(9) dots(0)
dis "A2_RC=" _rc
dis "RESULT A2 N=" e(N)
dis "RESULT A2 bootdf=" e(boot_df)

* ---- B: always-treated trigger ---------------------------------------
clear
set seed 55
qui set obs `=16*2*12'
qui gen int state = ceil(_n / 24)
qui gen int unit  = ceil(_n / 12)
qui bysort unit: gen int year = 1999 + _n
qui gen int onset = .
qui replace onset = 1996 if state <= 4
qui replace onset = 2004 + mod(state - 1, 3) if state >= 5 & state <= 8
qui gen double a_i = 0
qui bysort unit: replace a_i = rnormal() if _n == 1
qui bysort unit: replace a_i = a_i[1]
qui gen byte post = (onset < . & year >= onset)
qui gen double y = 1 + a_i + .6*post + rnormal()*.3
tempfile B
qui save `B'
dis "--- B1. no if ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(80) rseed(9) dots(0)
dis "B1_RC=" _rc
dis "RESULT B1 N=" e(N)
dis "RESULT B1 bootdf=" e(boot_df)
dis "--- B2. if y < . (the obvious hygiene if -- does NOT help here) ---"
use `B', clear
cap noi wooldid y unit year onset if y < . , cluster(state) bootstrap(standard) reps(80) rseed(9) dots(0)
dis "B2_RC=" _rc
dis "RESULT B2 N=" e(N)
dis "RESULT B2 bootdf=" e(boot_df)
dis "--- B3. if missing(onset) | onset >= 1999 (the remedy that works) ---"
use `B', clear
cap noi wooldid y unit year onset if missing(onset) | onset >= 1999 , cluster(state) bootstrap(standard) reps(80) rseed(9) dots(0)
dis "B3_RC=" _rc
dis "RESULT B3 N=" e(N)
dis "RESULT B3 bootdf=" e(boot_df)
dis "DONE_R11"
