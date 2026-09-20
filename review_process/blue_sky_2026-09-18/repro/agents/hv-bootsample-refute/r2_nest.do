*==========================================================================
* r2_nest.do -- SKEPTIC's own repro of HV-BOOT-06.
*
* CLAIM UNDER TEST: wooldid_bootrun decides the T2/T3 nesting regime on
* `smp' (user-if AND nonmissing cluster ONLY -- no markout), so a unit that
* straddles two clusters ONLY through a row the estimator drops is
* classified NON-nesting, the physical re-identification leg takes the
* keep-i rule, and att_i's composition variance collapses.
*
* ISOLATION (this is the part that makes it a test rather than a story):
*  - deleting the one stray row does NOT change G (state 8 exists anyway),
*    so with a fixed rseed() the bsample draw stream must be IDENTICAL.
*    bootsave() writes one w_<value> multiplicity column per cluster per
*    replication; we compare those matrices row-for-row.  If they match,
*    any SE difference is the regime flip, not Monte-Carlo noise.
*  - CONTROL: the same stray row carrying a REAL outcome.  Then the
*    straddle is genuine inside e(sample) and keep-i is the CORRECT rule,
*    so the control must sit with the full-data run, not with the oracle.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

*--------------------------------------------------------------- generator
* 16 states x 3 units x 10 years; i = unit, cluster = state (COARSER than i)
cap prog drop mknest
prog def mknest
  syntax , seed(integer)
  clear
  set seed `seed'
  qui set obs `=16*3*10'
  qui gen int unit  = ceil(_n / 10)
  qui gen int state = ceil(unit / 3)
  qui bysort unit: gen int year = 1999 + _n
  qui gen int onset = .
  qui replace onset = 2003 + mod(state - 1, 4) if state <= 8
  qui gen double a_i = 0
  qui bysort unit: replace a_i = rnormal() if _n == 1
  qui bysort unit: replace a_i = a_i[1]
  * treatment effect VARIES across units -> composition variance is real
  qui gen double tau_i = 0
  qui bysort unit: replace tau_i = .6 + rnormal()*.5 if _n == 1
  qui bysort unit: replace tau_i = tau_i[1]
  qui gen double g_t = 0
  qui bysort year: replace g_t = rnormal()*.4 if _n == 1
  qui bysort year (unit): replace g_t = g_t[1]
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 1 + a_i + g_t + tau_i*post + rnormal()*.3
end

cap prog drop rep2
prog def rep2
  syntax , tag(string)
  tempname B V
  cap mat `B' = e(b)
  cap mat `V' = e(V)
  dis "RESULT `tag' N="      e(N)
  dis "RESULT `tag' Nclust=" e(N_clust)
  dis "RESULT `tag' names="  "`: colfullnames `B''"
  local cn : colfullnames `B'
  local j 0
  foreach c of local cn {
    local ++j
    if "`c'" == "att_i" {
      local wb = el("`B'", 1, `j')
      local wv = el("`V'", `j', `j')
      dis "RESULT `tag' b_att_i="  %21x `wb'
      dis "RESULT `tag' se_att_i=" %21x sqrt(`wv')
      dis "RESULT `tag' se_dec="   %12.9f sqrt(`wv')
    }
  }
end

*==========================================================================
foreach ds in 101 202 303 {
foreach rs in 5551212 8675309 {

dis "===== DATASEED `ds'  RSEED `rs' ====="

*---- build the three files -------------------------------------------
mknest , seed(`ds')
tempfile base
qui save `base'

* STRAY: one extra record for unit 19 coded to state 8 (it belongs to 7),
* with the outcome MISSING -- a crosswalk / boundary-change artefact.
qui set obs `=_N+1'
qui replace unit  = 19   in L
qui replace state = 8    in L
qui replace year  = 2005 in L
qui replace onset = onset[19*10] in L
qui replace y     = .    in L
tempfile stray
qui save `stray'

* CONTROL: same stray record but with a REAL outcome (genuine straddle)
use `stray', clear
qui replace y = 1.4 in L
tempfile ctrl
qui save `ctrl'

*---- FULL (stray row, y missing) --------------------------------------
use `stray', clear
cap noi wooldid y unit year onset , att_i cluster(state) bootstrap(standard) ///
    reps(400) rseed(`rs') dots(0) bootsave(wfull_`ds'_`rs') replace
dis "FULL_RC=" _rc
rep2 , tag(FULL_`ds'_`rs')

*---- ORACLE (stray row physically deleted) -----------------------------
use `base', clear
cap noi wooldid y unit year onset , att_i cluster(state) bootstrap(standard) ///
    reps(400) rseed(`rs') dots(0) bootsave(worac_`ds'_`rs') replace
dis "ORAC_RC=" _rc
rep2 , tag(ORAC_`ds'_`rs')

*---- CONTROL (stray row with real y: genuine straddle, keep-i CORRECT) --
use `ctrl', clear
cap noi wooldid y unit year onset , att_i cluster(state) bootstrap(standard) ///
    reps(400) rseed(`rs') dots(0)
dis "CTRL_RC=" _rc
rep2 , tag(CTRL_`ds'_`rs')

*---- draw-stream identity check ---------------------------------------
* compare the w_ multiplicity columns of the two saved files row-for-row.
cap noi {
  use "wfull_`ds'_`rs'.dta", clear
  qui keep _rep w_*
  qui rename w_* f_*
  qui sort _rep
  tempfile wf
  qui save `wf'
  use "worac_`ds'_`rs'.dta", clear
  qui keep _rep w_*
  qui sort _rep
  qui merge 1:1 _rep using `wf', nogen
  qui gen double maxdiff = 0
  foreach v of varlist w_* {
    local sfx = substr("`v'", 3, .)
    cap qui replace maxdiff = max(maxdiff, abs(`v' - f_`sfx'))
  }
  qui sum maxdiff
  local wmx = r(max)
  dis "RESULT DRAWSTREAM_`ds'_`rs' maxabsdiff=" %21x `wmx'
  qui count
  dis "RESULT DRAWSTREAM_`ds'_`rs' nrows=" r(N)
}
dis "DRAW_RC=" _rc

}
}

dis "DONE_R2"
