*==========================================================================
* r5_size.do -- SKEPTIC's own repro of HV-BOOT-02 and HV-BOOT-03.
*  HV-BOOT-02: under bootstrap(bayes) bootalpha(size) the Gamma shape and
*     the row divisor are the cluster's UNSCREENED row count, so with a
*     live fraction f the reported SE is ~sqrt(f) of a bootstrap of the
*     estimation sample.
*  HV-BOOT-03: bootsave()'s _rep==0 row holds 1/N_frame, not 1/e(N).
*  CONTROL (mine, not theirs): bootalpha(equal) on the same data, where
*     Dirichlet aggregation says dead rows should NOT move the SE.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

cap prog drop rep5
prog def rep5
  syntax , tag(string)
  tempname B V
  cap mat `B' = e(b)
  cap mat `V' = e(V)
  local wb = el("`B'", 1, 1)
  local wv = el("`V'", 1, 1)
  dis "RESULT `tag' N="     e(N)
  dis "RESULT `tag' name1=" "`: colfullnames `B''"
  dis "RESULT `tag' b="     %21x `wb'
  dis "RESULT `tag' se="    %21x sqrt(`wv')
  dis "RESULT `tag' sedec=" %14.10f sqrt(`wv')
end

*--- 14 states x 2 units x 12 years; y MAR-missing on half the cells -------
clear
set seed 606
qui set obs `=14*2*12'
qui gen int state = ceil(_n / 24)
qui gen int unit  = ceil(_n / 12)
qui bysort unit: gen int year = 1999 + _n
qui gen int onset = .
qui replace onset = 2003 + mod(state - 1, 3) if state <= 7
qui gen double a_i = 0
qui bysort unit: replace a_i = rnormal() if _n == 1
qui bysort unit: replace a_i = a_i[1]
qui gen double g_t = 0
qui bysort year: replace g_t = rnormal()*.4 if _n == 1
qui bysort year (unit): replace g_t = g_t[1]
qui gen byte post = (onset < . & year >= onset)
qui gen double y = 1 + a_i + g_t + .6*post + rnormal()*.3
* missing at random on half the rows; NO cluster is wholly dead
qui gen double u = runiform()
qui replace y = . if u < .5
qui count
dis "RESULT FRAMEROWS=" r(N)
qui count if y < .
dis "RESULT LIVEROWS=" r(N)
tempfile D
qui save `D'

dis "--- size.full ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(bayes) bootalpha(size) ///
    reps(600) rseed(191919) dots(0) bootsave(szfull) replace
dis "SZF_RC=" _rc
rep5 , tag(SZ_FULL)

use `D', clear
qui drop if y >= .
dis "--- size.oracle ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(bayes) bootalpha(size) ///
    reps(600) rseed(191919) dots(0) bootsave(szorac) replace
dis "SZO_RC=" _rc
rep5 , tag(SZ_ORAC)

*--- CONTROL: bootalpha(equal), same data -------------------------------
use `D', clear
dis "--- equal.full ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(bayes) ///
    reps(600) rseed(191919) dots(0)
dis "EQF_RC=" _rc
rep5 , tag(EQ_FULL)

use `D', clear
qui drop if y >= .
dis "--- equal.oracle ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(bayes) ///
    reps(600) rseed(191919) dots(0)
dis "EQO_RC=" _rc
rep5 , tag(EQ_ORAC)

*--- CONTROL: bootstrap(standard), same data ----------------------------
use `D', clear
dis "--- std.full ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    reps(600) rseed(191919) dots(0)
dis "STF_RC=" _rc
rep5 , tag(ST_FULL)
use `D', clear
qui drop if y >= .
dis "--- std.oracle ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) ///
    reps(600) rseed(191919) dots(0)
dis "STO_RC=" _rc
rep5 , tag(ST_ORAC)

*--- HV-BOOT-03: the observed row of the bootsave file -------------------
dis "===== BOOTSAVE OBSERVED ROW ====="
cap noi {
  use "szfull.dta", clear
  qui ds w_*
  local wv1 : word 1 of `r(varlist)'
  qui sum `wv1' if _rep == 0
  local wobs = r(mean)
  dis "RESULT SAVE_FULL wcol=`wv1' obsrow=" %21x `wobs'
  dis "RESULT SAVE_FULL obsrow_dec=" %14.10f `wobs'
  dis "RESULT SAVE_FULL oneover=" %14.4f (1/`wobs')
  * the identity the help states: sum over rows of n_c*w_c = 1 in a replication
  qui ds w_*
  dis "RESULT SAVE_FULL ncols=" `: word count `r(varlist)''
}
dis "SAVE_FULL_RC=" _rc
cap noi {
  use "szorac.dta", clear
  qui ds w_*
  local wv1 : word 1 of `r(varlist)'
  qui sum `wv1' if _rep == 0
  local wobs = r(mean)
  dis "RESULT SAVE_ORAC wcol=`wv1' obsrow=" %21x `wobs'
  dis "RESULT SAVE_ORAC obsrow_dec=" %14.10f `wobs'
  dis "RESULT SAVE_ORAC oneover=" %14.4f (1/`wobs')
}
dis "SAVE_ORAC_RC=" _rc

dis "DONE_R5"
