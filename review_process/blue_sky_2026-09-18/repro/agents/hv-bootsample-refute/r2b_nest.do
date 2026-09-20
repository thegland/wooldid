*==========================================================================
* r2b_nest.do -- SKEPTIC's repro of HV-BOOT-06, repaired and sharpened.
* Fixes r2_nest.do (the stray row took its onset from the wrong physical
* row because the generator leaves the data sorted by year-then-unit, and
* the estimand column is named main_att_i, not att_i).
*
* THREE INDEPENDENT INSTRUMENTS, so this is not a story about one SE:
*  (1) DIRECT REGIME READ-OUT.  ado:406-407: the A10 caution ("a covariate
*      specification names the i variable ...") is printed IF AND ONLY IF
*      wdb_reid is on AND wdb_nest == 1.  Naming i inside fe() therefore
*      turns the caution into a literal print-out of the regime bit.
*  (2) DRAW-STREAM IDENTITY.  Deleting the stray row does not change G, so
*      with a fixed rseed() the bsample stream is unchanged; bootsave()'s
*      w_ multiplicity columns are compared row-for-row (r2 already
*      measured max|diff| = 0 exactly over 401 rows).
*  (3) The CONTROL: the same stray row carrying a real outcome, where the
*      straddle is genuine inside e(sample) and keep-i is the CORRECT rule.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

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
  qui gen double tau_i = 0
  qui bysort unit: replace tau_i = .6 + rnormal()*.5 if _n == 1
  qui bysort unit: replace tau_i = tau_i[1]
  qui gen double g_t = 0
  qui bysort year: replace g_t = rnormal()*.4 if _n == 1
  qui bysort year (unit): replace g_t = g_t[1]
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 1 + a_i + g_t + tau_i*post + rnormal()*.3
  qui sort unit year
end

cap prog drop rep2
prog def rep2
  syntax , tag(string)
  tempname B V
  cap mat `B' = e(b)
  cap mat `V' = e(V)
  dis "RESULT `tag' N="      e(N)
  dis "RESULT `tag' Nclust=" e(N_clust)
  local cn : colfullnames `B'
  local j 0
  foreach c of local cn {
    local ++j
    if "`c'" == "main_att_i" {
      local wb = el("`B'", 1, `j')
      local wv = el("`V'", `j', `j')
      dis "RESULT `tag' b="  %21x `wb'
      dis "RESULT `tag' se=" %21x sqrt(`wv')
      dis "RESULT `tag' sedec=" %14.10f sqrt(`wv')
    }
  }
end

*==========================================================================
* STEP 1: the regime read-out (reps kept tiny; only the caution matters)
*==========================================================================
dis "===== STEP 1: REGIME READ-OUT via the A10 caution ====="
mknest , seed(101)
tempfile base
qui save `base'
qui sum onset if unit == 19, meanonly
local o19 = r(mean)
dis "RESULT o19=" `o19'
qui set obs `=_N+1'
qui replace unit  = 19     in L
qui replace state = 8      in L
qui replace year  = 2005   in L
qui replace onset = `o19'  in L
qui replace y     = .      in L
tempfile stray
qui save `stray'
use `stray', clear
qui replace y = 1.4 in L
tempfile ctrl
qui save `ctrl'

dis "--- REGIME.oracle (no stray row): expect the A10 caution (nesting) ---"
use `base', clear
cap noi wooldid y unit year onset , att_i cluster(state) fe(unit) bootstrap(standard) reps(12) rseed(1) dots(0)
dis "REG_ORAC_RC=" _rc
dis "--- REGIME.full (stray row, y MISSING) ---"
use `stray', clear
cap noi wooldid y unit year onset , att_i cluster(state) fe(unit) bootstrap(standard) reps(12) rseed(1) dots(0)
dis "REG_FULL_RC=" _rc
dis "--- REGIME.control (stray row, y PRESENT) ---"
use `ctrl', clear
cap noi wooldid y unit year onset , att_i cluster(state) fe(unit) bootstrap(standard) reps(12) rseed(1) dots(0)
dis "REG_CTRL_RC=" _rc

*==========================================================================
* STEP 2: the standard errors, three seeds x two rseeds
*==========================================================================
dis "===== STEP 2: STANDARD ERRORS ====="
foreach ds in 101 202 303 {
foreach rs in 5551212 8675309 {

dis "===== DATASEED `ds'  RSEED `rs' ====="
mknest , seed(`ds')
tempfile base
qui save `base'
qui sum onset if unit == 19, meanonly
local o19 = r(mean)
qui set obs `=_N+1'
qui replace unit  = 19     in L
qui replace state = 8      in L
qui replace year  = 2005   in L
qui replace onset = `o19'  in L
qui replace y     = .      in L
tempfile stray
qui save `stray'
use `stray', clear
qui replace y = 1.4 in L
tempfile ctrl
qui save `ctrl'

use `stray', clear
cap noi wooldid y unit year onset , att_i cluster(state) bootstrap(standard) ///
    reps(400) rseed(`rs') dots(0) bootsave(bf_`ds'_`rs') replace
dis "FULL_RC=" _rc
rep2 , tag(FULL_`ds'_`rs')

use `base', clear
cap noi wooldid y unit year onset , att_i cluster(state) bootstrap(standard) ///
    reps(400) rseed(`rs') dots(0) bootsave(bo_`ds'_`rs') replace
dis "ORAC_RC=" _rc
rep2 , tag(ORAC_`ds'_`rs')

use `ctrl', clear
cap noi wooldid y unit year onset , att_i cluster(state) bootstrap(standard) ///
    reps(400) rseed(`rs') dots(0)
dis "CTRL_RC=" _rc
rep2 , tag(CTRL_`ds'_`rs')

cap noi {
  use "bf_`ds'_`rs'.dta", clear
  qui keep _rep w_*
  qui rename w_* f_*
  qui sort _rep
  tempfile wf
  qui save `wf'
  use "bo_`ds'_`rs'.dta", clear
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
  qui count
  local wnr = r(N)
  dis "RESULT DRAW_`ds'_`rs' maxabsdiff=" %14.10f `wmx' "  nrows=`wnr'"
}
dis "DRAW_RC=" _rc

}
}

dis "DONE_R2B"
