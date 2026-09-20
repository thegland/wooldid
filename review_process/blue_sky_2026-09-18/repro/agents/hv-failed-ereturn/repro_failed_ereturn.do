*=====================================================================
* repro_failed_ereturn.do  -- SELF-CONTAINED
*
* CLAIM: when a wooldid run fails AFTER its internal reghdfe/ppmlhdfe
* regression, nothing clears e().  The user's previous good wooldid
* results are replaced by the internal saturated model, silently.
*
* Run:  stata-mp -e do repro_failed_ereturn.do
* Grep the log for lines beginning "R|".
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

*--- an utterly ordinary staggered panel: 60 units, 10 years,
*    3 treated cohorts + 12 never-treated, one covariate, one dose ----
clear
set obs 60
gen int unit  = _n
gen byte grp  = mod(_n,5)
gen double onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i   = rnormal()
gen byte  state  = 1 + mod(_n,10)
gen double dose0 = 0.5 + runiform()
expand 10
bys unit: gen int year = 1999 + _n
gen byte   post = !missing(onset) & year >= onset
gen double x1   = rnormal()
gen double tau  = 0
replace    tau  = 0.8*(year-onset+1) if post
gen double y    = a_i + 0.10*(year-2000) + tau + 0.4*x1 + 0.30*rnormal()
gen double dosep = dose0*post            // dose, zero outside treatment
egen long sy = group(state year)         // region-by-year FE id

*=====================================================================
* 1. A GOOD RUN.  This is the result the user wants to keep.
*=====================================================================
wooldid y unit year onset, att
dis "R|good_cmd="     "`e(cmd)'"
dis "R|good_cmdline=" "`e(cmdline)'"
dis "R|good_att="     _b[main_att]
dis "R|good_N="       e(N)
dis "R|good_kb="      colsof(e(b))

*=====================================================================
* 2. A SECOND SPEC THAT FAILS AFTER THE INTERNAL REGRESSION.
*    Three independent realistic triggers; each clobbers e().
*=====================================================================
* (a) continuous treatment at unit-level i  -> r(322) from the
*     dropped-cell-dummy screen, which runs AFTER reghdfe
cap noisily wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
dis "R|failA_rc="      _rc
dis "R|failA_cmd="     "`e(cmd)'"
dis "R|failA_cmdline=" substr(`"`e(cmdline)'"',1,90)
dis "R|failA_kb="      colsof(e(b))
dis "R|failA_N="       e(N)
local nm : colnames e(b)
dis "R|failA_name2="   `"`: word 2 of `nm''"'

* restore a good run, then trigger (b)
qui wooldid y unit year onset, att
* (b) saturating region-by-year fixed effects -> r(322)
cap noisily wooldid y unit year onset, att fe(sy)
dis "R|failB_rc="  _rc
dis "R|failB_cmd=" "`e(cmd)'"
dis "R|failB_kb="  colsof(e(b))

qui wooldid y unit year onset, att
* (c) ccc_absorb() at unit-level i -> r(461) raised inside margins
cap noisily wooldid y unit year onset, att ccc_absorb(x1)
dis "R|failC_rc="  _rc
dis "R|failC_cmd=" "`e(cmd)'"
dis "R|failC_kb="  colsof(e(b))
dis "R|failC_N="   e(N)
dis "R|failC_b11=" e(b)[1,1]

*=====================================================================
* 3. WHAT THE USER LOSES.  Start from a good run, fail once, look.
*=====================================================================
qui wooldid y unit year onset, att
scalar truth = _b[main_att]
cap noisily wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)

* 3a. bare `wooldid' (the documented replay) no longer reaches its own
*     staleness message; it reports a command line the user never typed
cap noisily wooldid
dis "R|replay_rc=" _rc

* 3b. a scripted harvest silently takes reghdfe's numbers
qui wooldid y unit year onset, att
matrix RES = J(2,2,.)
matrix RES[1,1] = _b[main_att]
matrix RES[1,2] = e(N)
cap wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
matrix RES[2,1] = e(b)[1,1]
matrix RES[2,2] = e(N)
dis "R|harvest_row2_b=" RES[2,1] "  harvest_row2_N=" RES[2,2]
matrix list RES

* 3c. lincom/test still work, on the internal saturated model
local nm : colnames e(b)
local n1 : word 2 of `nm'
cap noisily lincom _b[`n1']
dis "R|lincom_rc=" _rc "  b=" r(estimate)

* 3d. predict dies on an internal tempvar name
cap noisily predict double xbhat
dis "R|predict_rc=" _rc

* 3e. esttab/estimates table on the "current" model
cap noisily estimates table, keep(main_att)
dis "R|esttab_keep_rc=" _rc

* 3f. the only recovery is a pre-emptive estimates store
qui wooldid y unit year onset, att
est store keepme
cap wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
est restore keepme
dis "R|recovered_cmd=" "`e(cmd)'"
dis "R|recovered_att=" _b[main_att] "  truth=" truth

*=====================================================================
* 4. CONTRAST: a PARSE-stage failure leaves e() alone, so the designed
*    staleness message is reachable there.  Same rc(198), different
*    diagnosis -- the user cannot tell the two classes apart.
*=====================================================================
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, attt
dis "R|parsefail_rc="  _rc
dis "R|parsefail_cmd=" "`e(cmd)'"
cap noisily wooldid
dis "R|replay_after_parse_rc=" _rc

*=====================================================================
* 5. BASELINE: this is not what Stata does, and the fix is one line.
*=====================================================================
sysuse auto, clear
prog def zz_noHold, eclass
  regress price weight foreign
  error 322
end
prog def zz_withHold, eclass
  tempname h
  _estimates hold `h', restore nullok
  regress price weight foreign
  error 322
end
regress price mpg
cap noisily regress price nosuchvar
dis "R|builtin_fail_cmdline=" "`e(cmdline)'"
regress price mpg
cap noisily zz_noHold
dis "R|ado_nohold_cmdline="  "`e(cmdline)'"
regress price mpg
cap noisily zz_withHold
dis "R|ado_withhold_cmdline=" "`e(cmdline)'"

dis "REPRO DONE"
