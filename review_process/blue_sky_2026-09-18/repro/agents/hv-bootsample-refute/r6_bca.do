*==========================================================================
* r6_bca.do -- SKEPTIC's own repro of HV-BOOT-07, with the isolation the
* original probe skipped.
*
* The BCa acceleration a = d3/(6 d2^1.5) is computed from the delete-one-
* cluster JACKKNIFE alone (ado:15786-15790).  It does NOT depend on the
* bootstrap draws.  e(boot_accel) is posted.  So comparing e(boot_accel)
* between the full data and the oracle is a MONTE-CARLO-FREE measurement
* of the defect -- unlike the interval comparison, which is confounded by
* the differing draw streams.
*
* Prediction if the claim is right: the jackknife runs G_frame fits, the
* 5 fits that delete a wholly-dead cluster reproduce the full estimate
* exactly, and e(boot_accel) therefore differs from the oracle's.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

cap prog drop rep6
prog def rep6
  syntax , tag(string)
  tempname A B
  cap mat `A' = e(boot_accel)
  cap mat `B' = e(b)
  dis "RESULT `tag' N=" e(N)
  dis "RESULT `tag' Nclust=" e(N_clust)
  dis "RESULT `tag' bootdf=" e(boot_df)
  cap {
    local wa = el("`A'", 1, 1)
    dis "RESULT `tag' accel=" %21x `wa'
    dis "RESULT `tag' accel_dec=" %14.10f `wa'
  }
  cap {
    local wb = el("`B'", 1, 1)
    dis "RESULT `tag' b=" %21x `wb'
  }
  cap mat li `A'
end

*--- 18 states x 2 units x 10 years; 5 states report no outcome ----------
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
tempfile base
qui save `base'
* five states (a mix of treated and never-treated) report NO outcome
qui replace y = . if inlist(state, 2, 6, 11, 14, 17)
tempfile full
qui save `full'

dis "--- bca.full (18 frame clusters, 13 live) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    bootdisplay(bca) reps(300) rseed(24680) dots(0)
dis "BCAF_RC=" _rc
rep6 , tag(BCA_FULL)

use `full', clear
qui drop if y >= .
dis "--- bca.oracle (13 clusters) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    bootdisplay(bca) reps(300) rseed(24680) dots(0)
dis "BCAO_RC=" _rc
rep6 , tag(BCA_ORAC)

*--- second determinism control: a DIFFERENT rseed on the full data.
*    If the acceleration is jackknife-only, it must be bit-identical to
*    BCA_FULL; that proves the full-vs-oracle gap is not draw noise.
use `full', clear
dis "--- bca.full, different rseed ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    bootdisplay(bca) reps(300) rseed(13579) dots(0)
dis "BCAF2_RC=" _rc
rep6 , tag(BCA_FULL2)

dis "DONE_R6"
