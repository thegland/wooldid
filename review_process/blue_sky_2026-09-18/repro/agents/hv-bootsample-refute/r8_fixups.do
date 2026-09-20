*==========================================================================
* r8_fixups.do
*  PART 1 (HV-BOOT-07): the BCa acceleration, printed directly.  The
*    acceleration is a pure function of the delete-one-cluster jackknife
*    (ado:15786-15790) and NOT of the bootstrap draws, so full-vs-oracle
*    is a Monte-Carlo-FREE comparison.  Control: the same full data with a
*    different rseed must give a bit-identical acceleration.
*  PART 2 (HV-VAL-01 siblings): the poisson y>=0 assert and the weight
*    assert, on fixtures that actually have a usable design.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

*==========================================================================
* PART 1: BCa acceleration
*==========================================================================
clear
set seed 1234
qui set obs `=18*2*10'
qui gen int state = ceil(_n / 20)
qui gen int unit  = ceil(_n / 10)
qui bysort unit: gen int year = 1999 + _n
qui gen int onset = .
qui replace onset = 2003 + mod(state - 1, 3) if state <= 9
qui gen double a_i = 0
qui bysort unit: replace a_i = rnormal() if _n == 1
qui bysort unit: replace a_i = a_i[1]
qui gen double tau_i = 0
qui bysort unit: replace tau_i = .6 + rnormal()*.4 if _n == 1
qui bysort unit: replace tau_i = tau_i[1]
qui gen double g_t = 0
qui bysort year: replace g_t = rnormal()*.4 if _n == 1
qui bysort year (unit): replace g_t = g_t[1]
qui gen byte post = (onset < . & year >= onset)
qui gen double y = 1 + a_i + g_t + tau_i*post + rnormal()*.3
qui replace y = . if inlist(state, 2, 6, 11, 14, 17)
tempfile full
qui save `full'

dis "--- P1.full (18 frame clusters, 13 live) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    bootdisplay(bca) reps(300) rseed(24680) dots(0)
dis "P1F_RC=" _rc
dis "RESULT P1F N=" e(N)
dis "RESULT P1F Nclust=" e(N_clust)
dis "RESULT P1F bootdf=" e(boot_df)
dis "RESULT P1F jkdone=" e(boot_jk_completed)
mat list e(boot_accel), format(%21x)
mat list e(boot_ci_bca), format(%12.8f)
mat list e(boot_p_bca), format(%12.8f)

dis "--- P1.full, DIFFERENT rseed (acceleration must be bit-identical) ---"
use `full', clear
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    bootdisplay(bca) reps(300) rseed(13579) dots(0)
dis "P1F2_RC=" _rc
dis "RESULT P1F2 jkdone=" e(boot_jk_completed)
mat list e(boot_accel), format(%21x)

dis "--- P1.oracle (13 clusters) ---"
use `full', clear
qui drop if y >= .
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    bootdisplay(bca) reps(300) rseed(24680) dots(0)
dis "P1O_RC=" _rc
dis "RESULT P1O N=" e(N)
dis "RESULT P1O Nclust=" e(N_clust)
dis "RESULT P1O bootdf=" e(boot_df)
dis "RESULT P1O jkdone=" e(boot_jk_completed)
mat list e(boot_accel), format(%21x)
mat list e(boot_ci_bca), format(%12.8f)
mat list e(boot_p_bca), format(%12.8f)

*==========================================================================
* PART 2: the poisson and weight asserts, on a valid staggered design
*==========================================================================
cap prog drop mkp2
prog def mkp2
  syntax , [seed(integer 31)]
  clear
  set seed `seed'
  qui set obs `=12*2*10'
  qui gen int state = ceil(_n / 20)
  qui gen int unit  = ceil(_n / 10)
  qui bysort unit: gen int year = 1999 + _n
  qui gen int onset = .
  qui replace onset = 2003 + mod(state - 1, 3) if state <= 6
  qui gen double a_i = 0
  qui bysort unit: replace a_i = rnormal() if _n == 1
  qui bysort unit: replace a_i = a_i[1]
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 4 + a_i + .6*post + rnormal()*.3
  qui gen double w = 1 + runiform()
end

dis "===== PART 2A: poisson y>=0 assert ====="
mkp2
tempfile base2
qui save `base2'
* five unusable records: unit id MISSING (they can never enter), y negative
qui set obs `=_N+5'
qui replace unit  = .    if unit >= .
qui replace state = 99   if state >= .
qui replace year  = 2001 if year >= .
qui replace y     = -3   if y >= .
qui replace w     = 1    if w >= .
tempfile negy
qui save `negy'
dis "--- 2A1. poisson, negative y on rows whose unit id is missing ---"
cap noi wooldid y unit year onset , cluster(state) poisson
dis "A1_RC=" _rc
dis "RESULT 2A1 N=" e(N)
dis "--- 2A2. same, safety(off) ---"
use `negy', clear
cap noi wooldid y unit year onset , cluster(state) poisson safety(off)
dis "A2_RC=" _rc
dis "RESULT 2A2 N=" e(N)
dis "--- 2A3. rows deleted ---"
use `negy', clear
qui drop if unit >= .
cap noi wooldid y unit year onset , cluster(state) poisson
dis "A3_RC=" _rc
dis "RESULT 2A3 N=" e(N)

dis "===== PART 2B: weight assert ====="
use `base2', clear
qui set obs `=_N+5'
qui replace unit  = 500  if unit >= .
qui replace state = 99   if state >= .
qui replace year  = 2001 if year >= .
qui replace y     = .    if y >= .
qui replace w     = -1   if w >= .
tempfile negw
qui save `negw'
dis "--- 2B1. aweight, negative weight on rows with NO outcome ---"
cap noi wooldid y unit year onset [aw=w] , cluster(state)
dis "B1_RC=" _rc
dis "RESULT 2B1 N=" e(N)
dis "--- 2B2. same, safety(off) ---"
use `negw', clear
cap noi wooldid y unit year onset [aw=w] , cluster(state) safety(off)
dis "B2_RC=" _rc
dis "RESULT 2B2 N=" e(N)
dis "--- 2B3. rows deleted ---"
use `negw', clear
qui drop if w < 0
cap noi wooldid y unit year onset [aw=w] , cluster(state)
dis "B3_RC=" _rc
dis "RESULT 2B3 N=" e(N)

dis "DONE_R8"
