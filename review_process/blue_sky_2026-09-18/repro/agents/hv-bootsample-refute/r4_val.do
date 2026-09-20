*==========================================================================
* r4_val.do -- SKEPTIC's own repro of HV-VAL-01.
* CLAIM: wooldid_validate's parse-time asserts run over the RAW dataset
* (ifassert is empty when no if() is given), so rows wooldid would itself
* have dropped can stop the run with r(411).
* Also probes: does safety(off) lift it?  does the message name a remedy?
* and the sibling asserts (poisson y>=0, weights>0).
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
  syntax , nst(integer) ntrt(integer) [nT(integer 10) t0(integer 2000) seed(integer 4242)]
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
  qui gen byte post = (onset < . & year >= onset)
  qui gen double y = 4 + a_i + .6*post + rnormal()*.3
  qui gen double w = 1 + runiform()
end

*==========================================================================
* A: 20 unusable records carrying the administrative sentinel unit == -9
*    with no outcome.  ifvar zeroes every y>=. row, so they could never
*    have entered the estimation sample.
*==========================================================================
dis "===== VAL A ====="
mkpanel , nst(12) ntrt(6) seed(31)
tempfile base
qui save `base'
qui set obs `=_N+20'
qui replace state = 99 if state >= .
qui replace unit  = -9 if unit  >= .
qui replace year  = 2000 + mod(_n, 10) if year >= .
qui replace onset = .  if unit == -9
qui replace y     = .  if unit == -9
qui replace w     = 1  if w >= .
qui gen byte junk = (unit != -9)
tempfile sent
qui save `sent'

dis "--- A1. full data, no if ---"
cap noi wooldid y unit year onset , cluster(state)
dis "A1_RC=" _rc
dis "RESULT A1 N=" e(N)

dis "--- A2. same data, user supplies if junk == 1 ---"
use `sent', clear
cap noi wooldid y unit year onset if junk == 1 , cluster(state)
dis "A2_RC=" _rc
dis "RESULT A2 N=" e(N)

dis "--- A3. same data, safety(off) ---"
use `sent', clear
cap noi wooldid y unit year onset , cluster(state) safety(off)
dis "A3_RC=" _rc
dis "RESULT A3 N=" e(N)

dis "--- A4. sentinel rows physically deleted ---"
use `sent', clear
qui drop if unit == -9
cap noi wooldid y unit year onset , cluster(state)
dis "A4_RC=" _rc
dis "RESULT A4 N=" e(N)

*==========================================================================
* B: sibling assert 1 -- poisson y >= 0 with a negative y on a row whose
*    UNIT id is missing (so the row can never enter the sample)
*==========================================================================
dis "===== VAL B ====="
use `base', clear
qui set obs `=_N+5'
qui replace unit  = .  if unit >= .
qui replace state = 99 if state >= .
qui replace year  = 2001 if year >= .
qui replace onset = .
qui replace y     = -3 if y >= .
tempfile negy
qui save `negy'
dis "--- B1. poisson, negative y on rows with missing unit ---"
cap noi wooldid y unit year onset , cluster(state) poisson
dis "B1_RC=" _rc
dis "RESULT B1 N=" e(N)
dis "--- B2. those rows deleted ---"
use `negy', clear
qui drop if unit >= .
cap noi wooldid y unit year onset , cluster(state) poisson
dis "B2_RC=" _rc
dis "RESULT B2 N=" e(N)

*==========================================================================
* C: sibling assert 2 -- a negative weight on a row with no outcome
*==========================================================================
dis "===== VAL C ====="
use `base', clear
qui set obs `=_N+5'
qui replace unit  = -9 if unit >= .
qui replace state = 99 if state >= .
qui replace year  = 2001 if year >= .
qui replace onset = .
qui replace y     = .  if y >= .
qui replace w     = -1 if w >= .
qui replace unit  = 500 if unit == -9
tempfile negw
qui save `negw'
dis "--- C1. aweight with a negative weight on a y-missing row ---"
cap noi wooldid y unit year onset [aw=w] , cluster(state)
dis "C1_RC=" _rc
dis "RESULT C1 N=" e(N)
dis "--- C2. those rows deleted ---"
use `negw', clear
qui drop if w < 0
cap noi wooldid y unit year onset [aw=w] , cluster(state)
dis "C2_RC=" _rc
dis "RESULT C2 N=" e(N)

dis "DONE_R4"
