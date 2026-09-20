*==========================================================================
* r10_alwaystreated.do -- deepens my own finding N1.
* A cohort treated BEFORE the estimation window is dropped by DOCUMENTED
* design (sthlp:73) with no missing data anywhere.  r7 showed it inflates
* e(boot_df) and the permute() design size.  Does it also:
*   (a) suppress the few-treated-cluster bootstrap caution?
*   (b) move the permutation p and the null SD?
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
  syntax , seed(integer) [tau(real 0.15)]
  clear
  set seed `seed'
  qui set obs `=20*2*12'
  qui gen int state = ceil(_n / 24)
  qui gen int unit  = ceil(_n / 12)
  qui bysort unit: gen int year = 1999 + _n
  qui gen int onset = .
  * states 1-6 treated BEFORE the window opens (1996); 7-10 treated inside;
  * 11-20 never treated.  EVERY variable present on EVERY row.
  qui replace onset = 1996 if state <= 6
  qui replace onset = 2004 + mod(state - 1, 3) if state >= 7 & state <= 10
  qui gen double a_i = 0
  qui bysort unit: replace a_i = rnormal() if _n == 1
  qui bysort unit: replace a_i = a_i[1]
  qui gen double g_t = 0
  qui bysort year: replace g_t = rnormal()*.4 if _n == 1
  qui bysort year (unit): replace g_t = g_t[1]
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 1 + a_i + g_t + `tau'*post + rnormal()*.5
  qui sort unit year
end

dis "===== A: BOOTSTRAP GUARDS ====="
mkp , seed(501)
tempfile A
qui save `A'
dis "--- A.full (6 always-treated states present) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(100) rseed(31) dots(0)
dis "A_FULL_RC=" _rc
dis "RESULT A_FULL N=" e(N)
dis "RESULT A_FULL Nclust=" e(N_clust)
dis "RESULT A_FULL bootdf=" e(boot_df)
dis "--- A.oracle (always-treated states deleted) ---"
use `A', clear
qui drop if onset == 1996
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(100) rseed(31) dots(0)
dis "A_ORAC_RC=" _rc
dis "RESULT A_ORAC N=" e(N)
dis "RESULT A_ORAC Nclust=" e(N_clust)
dis "RESULT A_ORAC bootdf=" e(boot_df)

dis "===== B: PERMUTATION p AND NULL SD ====="
foreach s in 501 602 703 {
  mkp , seed(`s')
  tempfile B
  qui save `B'
  cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(600) rseed(4004) dots(0)
  dis "B_FULL_`s'_RC=" _rc
  dis "RESULT B_FULL_`s' N=" e(N)
  dis "RESULT B_FULL_`s' npossible=" "`e(perm_npossible_str)'"
  dis "RESULT B_FULL_`s' p=" %12.8f e(perm_p_att)
  dis "RESULT B_FULL_`s' sd=" %12.8f e(perm_sd_att)
  use `B', clear
  qui drop if onset == 1996
  cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(600) rseed(4004) dots(0)
  dis "B_ORAC_`s'_RC=" _rc
  dis "RESULT B_ORAC_`s' N=" e(N)
  dis "RESULT B_ORAC_`s' npossible=" "`e(perm_npossible_str)'"
  dis "RESULT B_ORAC_`s' p=" %12.8f e(perm_p_att)
  dis "RESULT B_ORAC_`s' sd=" %12.8f e(perm_sd_att)
}
dis "DONE_R10"
