*==========================================================================
* r3b_perm.do -- repairs and extends r3_perm.do.
*  (i)  the plain (no-permute) control runs, which r3 mis-specified with
*       dots(0) -- dots() requires bootstrap()/permute().
*  (ii) a permonsetcohorts fixture that actually satisfies the bijection
*       (the never-treated must be ONE level of i).
*  (iii) THE QUESTION r3 COULD NOT ANSWER: does the reported permutation
*       p-value itself move between the two designs, or only the reported
*       design size?  r3's p sat on the add-one floor in both runs, which
*       is uninformative.  Here the effect is weak so p is interior.
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
  syntax , nst(integer) ntrt(integer) [nT(integer 10) t0(integer 2000) seed(integer 4242) tau(real 0.6)]
  clear
  set seed `seed'
  qui set obs `=`nst' * 2 * `nT''
  qui gen int state = ceil(_n / (2 * `nT'))
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
  qui gen double y = 1 + a_i + g_t + `tau'*post + rnormal()*.6
end

*==========================================================================
* P: does the reported p move between the 20-cluster and 16-cluster design?
*==========================================================================
dis "===== P: p-VALUE SENSITIVITY ====="
foreach s in 11 22 33 {
  mkpanel , nst(20) ntrt(8) seed(`s') tau(0.12)
  qui replace y = . if state <= 4
  tempfile P
  qui save `P'
  cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(600) rseed(90210) dots(0)
  dis "P_FULL_`s'_RC=" _rc
  dis "RESULT P_FULL_`s' N=" e(N)
  dis "RESULT P_FULL_`s' npossible=" "`e(perm_npossible_str)'"
  dis "RESULT P_FULL_`s' p=" %12.8f e(perm_p_att)
  dis "RESULT P_FULL_`s' sd=" %12.8f e(perm_sd_att)
  use `P', clear
  qui drop if y >= .
  cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(600) rseed(90210) dots(0)
  dis "P_ORAC_`s'_RC=" _rc
  dis "RESULT P_ORAC_`s' N=" e(N)
  dis "RESULT P_ORAC_`s' npossible=" "`e(perm_npossible_str)'"
  dis "RESULT P_ORAC_`s' p=" %12.8f e(perm_p_att)
  dis "RESULT P_ORAC_`s' sd=" %12.8f e(perm_sd_att)
}

*==========================================================================
* C: the three permctl screens, with WORKING plain-run controls
*==========================================================================
dis "===== C: PERMCTL SCREENS ====="
mkpanel , nst(8) ntrt(4) nT(10) seed(1717)
tempfile Cbase
qui save `Cbase'

* C1: one stray record for a unit in state 3 carrying a DIFFERENT onset, y missing
use `Cbase', clear
qui set obs `=_N+1'
qui replace state = 3    in L
qui replace unit  = 5    in L
qui replace year  = 2006 in L
qui replace onset = 2007 in L
qui replace y     = .    in L
tempfile Cnest
qui save `Cnest'
dis "--- C0. ORDINARY run on the stray file (no permute) ---"
cap noi wooldid y unit year onset , cluster(state)
dis "C0_RC=" _rc
dis "RESULT C0 N=" e(N)
dis "--- C0o. ORDINARY run with the stray row deleted ---"
use `Cnest', clear
qui drop if y >= .
cap noi wooldid y unit year onset , cluster(state)
dis "C0O_RC=" _rc
dis "RESULT C0O N=" e(N)

dis "--- C1. permute() on the stray file ---"
use `Cnest', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0)
dis "C1_RC=" _rc
dis "--- C1s. safety(off) ---"
use `Cnest', clear
cap noi wooldid y unit year onset , cluster(state) permute(shuffle) permreps(50) rseed(11) dots(0) safety(off)
dis "C1S_RC=" _rc

* C2: extended missing
use `Cbase', clear
qui set obs `=_N+1'
qui replace state = 3    in L
qui replace unit  = 5    in L
qui replace year  = 2006 in L
qui replace onset = .a   in L
qui replace y     = .    in L
tempfile Cext
qui save `Cext'
dis "--- C2ord. ORDINARY run on the .a file ---"
cap noi wooldid y unit year onset , cluster(state)
dis "C2ORD_RC=" _rc
dis "RESULT C2ORD N=" e(N)

*==========================================================================
* C3: permonsetcohorts -- i IS the onset grouping (never-treated = ONE level)
*==========================================================================
dis "===== C3: PERMONSETCOHORTS ====="
clear
qui set obs `=3*6*12'
* coh 1 (onset 2003), coh 2 (onset 2005), coh 3 (never treated) -- 6 hh each
qui gen int coh = ceil(_n / 72)
qui gen int hh  = ceil(_n / 12)
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
dis "--- C3clean ---"
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0)
dis "C3CLEAN_RC=" _rc
dis "RESULT C3CLEAN N=" e(N)

* one dead row giving a never-treated cohort a stray onset
use `C3base', clear
qui set obs `=_N+1'
qui replace coh   = 3    in L
qui replace hh    = 15   in L
qui replace year  = 2006 in L
qui replace onset = 2008 in L
qui replace y     = .    in L
tempfile C3s
qui save `C3s'
dis "--- C3ord. ORDINARY run on the stray file ---"
cap noi wooldid y coh year onset , cluster(hh)
dis "C3ORD_RC=" _rc
dis "RESULT C3ORD N=" e(N)
dis "--- C3. permonsetcohorts on the stray file ---"
use `C3s', clear
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0)
dis "C3_RC=" _rc
dis "--- C3s. safety(off) ---"
use `C3s', clear
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0) safety(off)
dis "C3S_RC=" _rc
dis "--- C3o. stray row deleted ---"
use `C3s', clear
qui drop if y >= .
cap noi wooldid y coh year onset , cluster(hh) permute(shuffle) permonsetcohorts permreps(50) rseed(11) dots(0)
dis "C3O_RC=" _rc
dis "RESULT C3O N=" e(N)

dis "DONE_R3B"
