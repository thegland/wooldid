*==========================================================================
* r3_perm.do -- SKEPTIC's own repro of HV-PERM-01 / 02 / 03.
* Self-contained; own data.
*   A: does the permute() design range over clusters with no rows in
*      e(sample)?  (permsave column count, npossible, treated caution, p)
*   B: permexact -- is the "no Monte-Carlo error" claim printed WITHOUT
*      the reduced-reference-set qualification in the same footer?
*   C: do the three wooldid_permctl screens refuse r(459) on rows wooldid
*      itself drops?  Does safety(off) lift any of them?
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , nst(integer) ntrt(integer) [nT(integer 10) t0(integer 2000) seed(integer 4242) nper(integer 2)]
  clear
  set seed `seed'
  qui set obs `=`nst' * `nper' * `nT''
  qui gen int state = ceil(_n / (`nper' * `nT'))
  qui gen int unit  = ceil(_n / `nT')
  qui bysort unit: gen int year = `t0' + _n - 1
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

cap prog drop rep3
prog def rep3
  syntax , tag(string)
  dis "RESULT `tag' N="        e(N)
  dis "RESULT `tag' npossible=" "`e(perm_npossible_str)'"
  dis "RESULT `tag' completed=" e(perm_completed)
  dis "RESULT `tag' reps="     e(perm_reps)
  dis "RESULT `tag' exact="    "`e(perm_exact)'"
  dis "RESULT `tag' incomplete=" "`e(perm_incomplete)'"
  dis "RESULT `tag' p_att="    %12.8f e(perm_p_att)
  cap {
    tempname BB
    mat `BB' = e(b)
    local wk = colnumb(`BB', "main_att")
    local wb = el("`BB'", 1, `wk')
    dis "RESULT `tag' b=" %21x `wb'
  }
end

*==========================================================================
* A: 20 states, 8 treated; states 1-4 (all treated) report NO outcome
*==========================================================================
dis "===== PERM A ====="
mkpanel , nst(20) ntrt(8) seed(4242)
qui replace y = . if state <= 4
tempfile A
qui save `A'

dis "--- A.full ---"
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(300) ///
   rseed(246810) dots(0) permsave(pA_full) replace
dis "A_FULL_RC=" _rc
rep3 , tag(A_FULL)
cap noi {
  preserve
  use "pA_full.dta", clear
  qui ds onset_*
  dis "RESULT A_FULL onsetcols=" `: word count `r(varlist)''
  dis "RESULT A_FULL onsetnames=" "`r(varlist)'"
  restore
}

use `A', clear
qui drop if y >= .
dis "--- A.oracle ---"
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(300) ///
   rseed(246810) dots(0) permsave(pA_orac) replace
dis "A_ORAC_RC=" _rc
rep3 , tag(A_ORAC)
cap noi {
  preserve
  use "pA_orac.dta", clear
  qui ds onset_*
  dis "RESULT A_ORAC onsetcols=" `: word count `r(varlist)''
  dis "RESULT A_ORAC onsetnames=" "`r(varlist)'"
  restore
}

*==========================================================================
* B: permexact on a 7-cluster design; clusters 1 and 7 report no outcome
*==========================================================================
dis "===== PERM B ====="
clear
qui set obs `=7*10'
qui gen int unit = ceil(_n / 10)
qui bysort unit: gen int year = 1999 + _n
qui gen int onset = .
qui replace onset = 2003 if unit == 1
qui replace onset = 2004 if unit == 2
set seed 8080
qui gen double a_i = 0
qui bysort unit: replace a_i = rnormal() if _n == 1
qui bysort unit: replace a_i = a_i[1]
qui gen byte post = (onset < . & year >= onset)
qui gen double y = 1 + a_i + .3*(year - 1999) + .6*post + rnormal()*.3
tempfile Bd
qui save `Bd'

qui replace y = . if inlist(unit, 1, 7)
tempfile Bfull
qui save `Bfull'

dis "--- B.full (permexact, units 1 and 7 dead) ---"
cap noi wooldid y unit year onset , cluster(unit) permute(shuffle) permexact dots(0)
dis "B_FULL_RC=" _rc
rep3 , tag(B_FULL)

use `Bfull', clear
qui drop if y >= .
dis "--- B.oracle ---"
cap noi wooldid y unit year onset , cluster(unit) permute(shuffle) permexact dots(0)
dis "B_ORAC_RC=" _rc
rep3 , tag(B_ORAC)

*==========================================================================
* C: the three permctl screens, each tripped by ONE row wooldid drops
*==========================================================================
dis "===== PERM C ====="
mkpanel , nst(8) ntrt(4) nT(10) seed(1717)
tempfile Cbase
qui save `Cbase'

dis "--- C0. ordinary (no permute) run on the STRAY file must succeed ---"
* stray record: one row for a unit in state 3 carrying a DIFFERENT onset, y missing
use `Cbase', clear
qui set obs `=_N+1'
qui replace state = 3    in L
qui replace unit  = 5    in L
qui replace year  = 2006 in L
qui replace onset = 2007 in L
qui replace y     = .    in L
tempfile Cnest
qui save `Cnest'
cap noi wooldid y unit year onset , cluster(state)
dis "C0_RC=" _rc
dis "RESULT C0 N=" e(N)

dis "--- C1. permute() on the same file ---"
use `Cnest', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0)
dis "C1_RC=" _rc
dis "--- C1s. same, safety(off) ---"
use `Cnest', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0) safety(off)
dis "C1S_RC=" _rc
dis "--- C1o. stray row deleted ---"
use `Cnest', clear
qui drop if y >= .
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0)
dis "C1O_RC=" _rc
dis "RESULT C1O N=" e(N)

dis "--- C2. extended missing (.a) on a dead row ---"
use `Cbase', clear
qui set obs `=_N+1'
qui replace state = 3    in L
qui replace unit  = 5    in L
qui replace year  = 2006 in L
qui replace onset = .a   in L
qui replace y     = .    in L
tempfile Cext
qui save `Cext'
cap noi wooldid y unit year onset , cluster(state)
dis "C2ORD_RC=" _rc
dis "RESULT C2ORD N=" e(N)
use `Cext', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0)
dis "C2_RC=" _rc
use `Cext', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0) safety(off)
dis "C2S_RC=" _rc

dis "--- C3. permonsetcohorts with one dead stray-onset row ---"
* i = cohort code that IS the onset grouping; cluster BELOW i
clear
qui set obs `=6*4*10'
qui gen int coh  = ceil(_n / 40)
qui gen int hh   = ceil(_n / 10)
qui bysort hh: gen int year = 1999 + _n
qui gen int onset = .
qui replace onset = 2003 if coh == 1
qui replace onset = 2005 if coh == 2
set seed 4321
qui gen double a_i = 0
qui bysort hh: replace a_i = rnormal() if _n == 1
qui bysort hh: replace a_i = a_i[1]
qui gen byte post = (onset < . & year >= onset)
qui gen double y = 1 + a_i + .2*(year - 1999) + .6*post + rnormal()*.3
tempfile C3base
qui save `C3base'
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0)
dis "C3CLEAN_RC=" _rc
dis "RESULT C3CLEAN N=" e(N)

use `C3base', clear
qui set obs `=_N+1'
qui replace coh   = 3    in L
qui replace hh    = 9    in L
qui replace year  = 2006 in L
qui replace onset = 2008 in L
qui replace y     = .    in L
tempfile C3s
qui save `C3s'
cap noi wooldid y coh year onset , cluster(hh)
dis "C3ORD_RC=" _rc
dis "RESULT C3ORD N=" e(N)
use `C3s', clear
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0)
dis "C3_RC=" _rc
use `C3s', clear
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0) safety(off)
dis "C3S_RC=" _rc
use `C3s', clear
qui drop if y >= .
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0)
dis "C3O_RC=" _rc
dis "RESULT C3O N=" e(N)

dis "DONE_R3"
