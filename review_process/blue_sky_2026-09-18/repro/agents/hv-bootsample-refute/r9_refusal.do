*==========================================================================
* r9_refusal.do -- two questions HV-BOOT-01 raises but does not answer.
*  Q1. The claim says the r(459) "at least 2 clusters in the estimation
*      sample" REFUSAL fails to fire.  Does anything else catch the run?
*      Fixture: 12 frame clusters, exactly ONE with any outcome data (it
*      holds both a treated and a never-treated unit, so the design is
*      otherwise estimable).
*  Q2. Does the wrong e(boot_df) actually move a number a user sees?
*      Default bootdisplay is percentile (rank-based, df-free), so the
*      exposure is e(V)+e(df_r) downstream.  Run lincom by hand, exactly
*      as the run's own footer tells the user to.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

cap prog drop mkp
prog def mkp
  syntax , nst(integer) ntrt(integer) [nT(integer 10) t0(integer 2000) seed(integer 4242) nper(integer 4)]
  clear
  set seed `seed'
  qui set obs `=`nst' * `nper' * `nT''
  qui gen int state = ceil(_n / (`nper' * `nT'))
  qui gen int unit  = ceil(_n / `nT')
  qui gen int wu    = mod(unit - 1, `nper') + 1
  qui bysort unit: gen int year = `t0' + _n - 1
  qui gen int onset = .
  * WITHIN each state: half the units treated, half never -- so one state alone
  * still identifies the design
  qui replace onset = `t0' + 4 + mod(state - 1, 3) if wu <= 2
  qui gen double a_i = 0
  qui bysort unit: replace a_i = rnormal() if _n == 1
  qui bysort unit: replace a_i = a_i[1]
  qui gen double g_t = 0
  qui bysort year: replace g_t = rnormal()*.4 if _n == 1
  qui bysort year (unit): replace g_t = g_t[1]
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 1 + a_i + g_t + .6*post + rnormal()*.3
end

*==========================================================================
dis "===== Q1: ONE live cluster out of 12 ====="
mkp , nst(12) ntrt(12) seed(1212)
qui replace y = . if state != 4
tempfile Q1
qui save `Q1'
dis "--- Q1a. bootstrap on the full file (12 frame clusters, 1 live) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(50) rseed(2) dots(0)
dis "Q1A_RC=" _rc
dis "RESULT Q1A N=" e(N)
dis "--- Q1b. same file with the dead rows deleted ---"
use `Q1', clear
qui drop if y >= .
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(50) rseed(2) dots(0)
dis "Q1B_RC=" _rc
dis "RESULT Q1B N=" e(N)
dis "--- Q1c. ordinary (analytic) run on the full file ---"
use `Q1', clear
cap noi wooldid y unit year onset , cluster(state)
dis "Q1C_RC=" _rc
dis "RESULT Q1C N=" e(N)

*==========================================================================
dis "===== Q2: does e(boot_df) move a number the user reads? ====="
mkp , nst(20) ntrt(20) seed(3131)
qui replace y = . if state <= 5
tempfile Q2
qui save `Q2'
dis "--- Q2a. full (20 frame clusters, 15 live) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(400) rseed(606) dots(0)
dis "Q2A_RC=" _rc
dis "RESULT Q2A N="      e(N)
dis "RESULT Q2A Nclust=" e(N_clust)
dis "RESULT Q2A bootdf=" e(boot_df)
cap noi lincom main_att
dis "RESULT Q2A lincom_p=" %14.10f r(p)
dis "RESULT Q2A lincom_lb=" %14.10f r(lb)
dis "RESULT Q2A lincom_ub=" %14.10f r(ub)
dis "RESULT Q2A lincom_se=" %14.10f r(se)

dis "--- Q2b. oracle (15 clusters) ---"
use `Q2', clear
qui drop if y >= .
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(400) rseed(606) dots(0)
dis "Q2B_RC=" _rc
dis "RESULT Q2B N="      e(N)
dis "RESULT Q2B Nclust=" e(N_clust)
dis "RESULT Q2B bootdf=" e(boot_df)
cap noi lincom main_att
dis "RESULT Q2B lincom_p=" %14.10f r(p)
dis "RESULT Q2B lincom_lb=" %14.10f r(lb)
dis "RESULT Q2B lincom_ub=" %14.10f r(ub)
dis "RESULT Q2B lincom_se=" %14.10f r(se)

dis "--- Q2c. full data forced to the CORRECT df by hand: bootdf(14) ---"
use `Q2', clear
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(400) rseed(606) dots(0) bootdf(14)
dis "Q2C_RC=" _rc
dis "RESULT Q2C bootdf=" e(boot_df)
cap noi lincom main_att
dis "RESULT Q2C lincom_p=" %14.10f r(p)
dis "RESULT Q2C lincom_lb=" %14.10f r(lb)
dis "RESULT Q2C lincom_ub=" %14.10f r(ub)
dis "RESULT Q2C lincom_se=" %14.10f r(se)

dis "DONE_R9"
