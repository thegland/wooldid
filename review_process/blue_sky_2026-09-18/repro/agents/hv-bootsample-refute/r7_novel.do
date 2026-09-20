*==========================================================================
* r7_novel.do -- MY OWN probes on the same axis, different triggers.
* Nobody suggested these.
*  A: ALWAYS-TREATED clusters (onset BEFORE the sample window).  wooldid
*     drops those cohorts by documented design (sthlp:73).  Are they still
*     bootstrap resampling units, and do they still silence the guards?
*  B: ZERO survey weights.  ado:686 drops w <= 0 from the sample; the
*     bootstrap frame keeps the rows.  Realistic: zero-weight records are
*     routine in weighted administrative extracts.
*  C: an EVENT-STUDY WINDOW that drops rows (ado:9367's own comment names
*     it as a trigger and the earlier probe could not reproduce it).
*  D: does the same frame-vs-sample gap reach impvalues' cluster count?
*  E: permstrata() -- is a stratum built out of clusters with no data?
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
  syntax , nst(integer) [nT(integer 12) t0(integer 2000) seed(integer 4242)]
  clear
  set seed `seed'
  qui set obs `=`nst' * 2 * `nT''
  qui gen int state = ceil(_n / (2 * `nT'))
  qui gen int unit  = ceil(_n / `nT')
  qui bysort unit: gen int year = `t0' + _n - 1
  qui gen int onset = .
  qui gen double a_i = 0
  qui bysort unit: replace a_i = rnormal() if _n == 1
  qui bysort unit: replace a_i = a_i[1]
  qui gen double g_t = 0
  qui bysort year: replace g_t = rnormal()*.4 if _n == 1
  qui bysort year (unit): replace g_t = g_t[1]
  qui gen double w = 1
end

cap prog drop fin
prog def fin
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 1 + a_i + g_t + .6*post + rnormal()*.3
end

cap prog drop rep7
prog def rep7
  syntax , tag(string)
  dis "RESULT `tag' N="      e(N)
  dis "RESULT `tag' Nclust=" e(N_clust)
  dis "RESULT `tag' bootdf=" e(boot_df)
  dis "RESULT `tag' dfr="    e(df_r)
end

*==========================================================================
* A: always-treated clusters (onset 1996, sample starts 2000)
*==========================================================================
dis "===== A: ALWAYS-TREATED CLUSTERS ====="
mkp , nst(16) seed(55)
* states 1-4 treated BEFORE the window opens; 5-8 treated inside; 9-16 never
qui replace onset = 1996 if state <= 4
qui replace onset = 2004 + mod(state - 1, 3) if state >= 5 & state <= 8
fin
tempfile A
qui save `A'

dis "--- A.full (always-treated states present) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(100) rseed(4242) dots(0)
dis "A_FULL_RC=" _rc
rep7 , tag(A_FULL)

use `A', clear
qui drop if onset == 1996
dis "--- A.oracle (always-treated states deleted) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(100) rseed(4242) dots(0)
dis "A_ORAC_RC=" _rc
rep7 , tag(A_ORAC)

dis "--- A.perm.full ---"
use `A', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(60) rseed(77) dots(0)
dis "AP_FULL_RC=" _rc
dis "RESULT AP_FULL N=" e(N)
dis "RESULT AP_FULL npossible=" "`e(perm_npossible_str)'"
dis "--- A.perm.oracle ---"
use `A', clear
qui drop if onset == 1996
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(60) rseed(77) dots(0)
dis "AP_ORAC_RC=" _rc
dis "RESULT AP_ORAC N=" e(N)
dis "RESULT AP_ORAC npossible=" "`e(perm_npossible_str)'"

*==========================================================================
* B: zero survey weights on two whole clusters
*==========================================================================
dis "===== B: ZERO WEIGHTS ====="
mkp , nst(14) seed(66)
qui replace onset = 2004 + mod(state - 1, 3) if state <= 7
fin
qui replace w = 0 if inlist(state, 2, 5)
tempfile B
qui save `B'

dis "--- B.full (two zero-weight states) ---"
cap noi wooldid y unit year onset [aw=w] , cluster(state) bootstrap(standard) reps(100) rseed(4242) dots(0)
dis "B_FULL_RC=" _rc
rep7 , tag(B_FULL)

use `B', clear
qui drop if w == 0
dis "--- B.oracle ---"
cap noi wooldid y unit year onset [aw=w] , cluster(state) bootstrap(standard) reps(100) rseed(4242) dots(0)
dis "B_ORAC_RC=" _rc
rep7 , tag(B_ORAC)

* the same data under bootalpha(size): is n_c the frame row count?
use `B', clear
dis "--- B.size.full ---"
cap noi wooldid y unit year onset [aw=w] , cluster(state) bootstrap(bayes) bootalpha(size) reps(200) rseed(4242) dots(0)
dis "BS_FULL_RC=" _rc
rep7 , tag(BS_FULL)

*==========================================================================
* C: an event-study window that removes rows from the estimation sample
*==========================================================================
dis "===== C: EVENT-STUDY WINDOW ====="
mkp , nst(12) nT(20) t0(1995) seed(77)
qui replace onset = 2010 if state <= 6
fin
tempfile C
qui save `C'
dis "--- C.nowindow ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(60) rseed(4242) dots(0)
dis "C_NOW_RC=" _rc
rep7 , tag(C_NOW)
use `C', clear
dis "--- C.window esprelength(2) espostlength(2) ---"
cap noi wooldid y unit year onset , cluster(state) esprelength(2) espostlength(2) bootstrap(standard) reps(60) rseed(4242) dots(0)
dis "C_WIN_RC=" _rc
rep7 , tag(C_WIN)

*==========================================================================
* D: impvalues with a wholly dead cluster
*==========================================================================
dis "===== D: IMPVALUES ====="
mkp , nst(14) seed(88)
qui replace onset = 2004 + mod(state - 1, 3) if state <= 7
fin
qui replace y = . if inlist(state, 3, 10)
tempfile D
qui save `D'
dis "--- D.full ---"
cap noi wooldid y unit year onset , cluster(state) impvalues
dis "D_FULL_RC=" _rc
dis "RESULT D_FULL N=" e(N)
dis "RESULT D_FULL Nclust=" e(N_clust)
dis "RESULT D_FULL dfr=" e(df_r)
use `D', clear
qui drop if y >= .
dis "--- D.oracle ---"
cap noi wooldid y unit year onset , cluster(state) impvalues
dis "D_ORAC_RC=" _rc
dis "RESULT D_ORAC N=" e(N)
dis "RESULT D_ORAC Nclust=" e(N_clust)
dis "RESULT D_ORAC dfr=" e(df_r)

*==========================================================================
* E: permstrata() with a stratum made only of clusters that carry no data
*==========================================================================
dis "===== E: PERMSTRATA ====="
mkp , nst(16) seed(99)
qui replace onset = 2004 + mod(state - 1, 3) if state <= 8
fin
qui gen byte region = 1 + (state > 8)
* region 1 holds states 1-8 (all treated); make states 1-4 carry no outcome
qui replace y = . if state <= 4
tempfile E
qui save `E'
dis "--- E.full ---"
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permstrata(region) permreps(60) rseed(77) dots(0)
dis "E_FULL_RC=" _rc
dis "RESULT E_FULL N=" e(N)
dis "RESULT E_FULL npossible=" "`e(perm_npossible_str)'"
use `E', clear
qui drop if y >= .
dis "--- E.oracle ---"
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permstrata(region) permreps(60) rseed(77) dots(0)
dis "E_ORAC_RC=" _rc
dis "RESULT E_ORAC N=" e(N)
dis "RESULT E_ORAC npossible=" "`e(perm_npossible_str)'"

dis "DONE_R7"
