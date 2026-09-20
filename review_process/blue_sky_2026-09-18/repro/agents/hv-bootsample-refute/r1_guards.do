*==========================================================================
* r1_guards.do -- SKEPTIC's own repro of HV-BOOT-01 / 04 / 05.
* Self-contained: generates its own data, no repo data, no repo generator.
* Question: does the bootstrap engine's cluster universe include clusters
* that contribute ZERO rows to e(sample), and does that (a) suppress the
* few-cluster guards, (b) set e(boot_df)/e(df_r) from the frame count?
* Control: does it also MOVE the standard error? (prober did not claim this)
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

*--------------------------------------------------------------- generator
* gN_states states, 2 units per state, gT years.  first gK states treated.
cap prog drop mkpanel
prog def mkpanel
  syntax , nst(integer) ntrt(integer) [nT(integer 12) t0(integer 2000) seed(integer 4242)]
  clear
  set seed `seed'
  local nrow = `nst' * 2 * `nT'
  qui set obs `nrow'
  qui gen int state = ceil(_n / (2 * `nT'))
  qui gen int unit  = (state - 1) * 2 + 1 + mod(ceil(_n / `nT') - 1, 2)
  qui bysort unit: gen int year = `t0' + _n - 1
  * onsets: treated states cycle through 4 onset dates; the rest never treated
  qui gen int onset = .
  qui replace onset = `t0' + 4 + mod(state - 1, 4) if state <= `ntrt'
  qui gen double a_i = 0
  qui bysort unit: replace a_i = rnormal() if _n == 1
  qui bysort unit: replace a_i = a_i[1]
  qui gen double g_t = 0
  qui bysort year: replace g_t = rnormal()*.4 if _n == 1
  qui bysort year (unit): replace g_t = g_t[1]
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 1 + a_i + g_t + .6*post + rnormal()*.3
end

*--------------------------------------------------------------- reporter
cap prog drop rep1
prog def rep1
  syntax , tag(string)
  tempname B V
  cap mat `B' = e(b)
  cap mat `V' = e(V)
  dis "RESULT `tag' N="       e(N)
  dis "RESULT `tag' Nclust="  e(N_clust)
  dis "RESULT `tag' bootdf="  e(boot_df)
  dis "RESULT `tag' dfr="     e(df_r)
  dis "RESULT `tag' col1=" "`: colfullnames `B''"
  cap dis "RESULT `tag' b="   %21x `B'[1,1]
  cap dis "RESULT `tag' se="  %21x sqrt(`V'[1,1])
end

*==========================================================================
* CASE A: 20 states (8 treated), states 1-4 (all treated) report NO outcome
*==========================================================================
dis "===== CASE A ====="
mkpanel , nst(20) ntrt(8) seed(4242)
qui replace y = . if state <= 4
tempfile A
qui save `A'

dis "--- A.full (dead rows present) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(400) rseed(778899) dots(0)
dis "A_FULL_RC=" _rc
rep1 , tag(A_FULL)

use `A', clear
qui drop if y >= .
dis "--- A.oracle (dead rows physically dropped) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(400) rseed(778899) dots(0)
dis "A_ORAC_RC=" _rc
rep1 , tag(A_ORAC)

* A.bayes equal -- advisor's bound: dead clusters should barely move this SE
use `A', clear
dis "--- A.full bayes(equal) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(bayes) reps(400) rseed(778899) dots(0) safety(off)
dis "AB_FULL_RC=" _rc
rep1 , tag(AB_FULL)
use `A', clear
qui drop if y >= .
dis "--- A.oracle bayes(equal) ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(bayes) reps(400) rseed(778899) dots(0) safety(off)
dis "AB_ORAC_RC=" _rc
rep1 , tag(AB_ORAC)

*==========================================================================
* CASE B: singleton unit, EVERY variable present.  No missing data at all.
*==========================================================================
dis "===== CASE B ====="
mkpanel , nst(10) ntrt(4) nT(8) seed(777)
tempfile Bbase
qui save `Bbase'
* append one extra state holding one unit observed in exactly one period
qui set obs `=_N+1'
qui replace state = 11       in L
qui replace unit  = 999      in L
qui replace year  = 2003     in L
qui replace onset = .        in L
qui replace y     = 1.25     in L
qui replace post  = 0        in L
tempfile B
qui save `B'

dis "--- B.with singleton row ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(200) rseed(31415) dots(0)
dis "B_WITH_RC=" _rc
rep1 , tag(B_WITH)

use `Bbase', clear
dis "--- B.without singleton row ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(200) rseed(31415) dots(0)
dis "B_WOUT_RC=" _rc
rep1 , tag(B_WOUT)

*==========================================================================
* CASE C: 12 states, only 2 of them report an outcome (1 treated, 1 never)
*==========================================================================
dis "===== CASE C ====="
mkpanel , nst(12) ntrt(6) nT(10) seed(999)
qui replace y = . if !inlist(state, 3, 12)
tempfile C
qui save `C'

dis "--- C.full ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(200) rseed(5150) dots(0) bootdisplay(wald)
dis "C_FULL_RC=" _rc
rep1 , tag(C_FULL)

use `C', clear
qui drop if y >= .
dis "--- C.oracle ---"
cap noi wooldid y unit year onset , cluster(state) bootstrap(standard) reps(200) rseed(5150) dots(0) bootdisplay(wald)
dis "C_ORAC_RC=" _rc
rep1 , tag(C_ORAC)

dis "DONE_R1"
