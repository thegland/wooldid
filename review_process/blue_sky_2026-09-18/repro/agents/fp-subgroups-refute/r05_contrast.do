* ==========================================================================
* r05_contrast.do -- does the missing-subgroup share weighting reach the
* number people actually REPORT in a heterogeneity analysis: the pairwise
* subgroup CONTRAST in e(pwsgs)?
* (The across-subgroup equality test is NOT evidence here: the e() name probed
*  below, e(joint_sg_att_i), does not exist -- the real one is
*  e(joint_att_i_allsgs) -- and in any case both arms have a true nonzero
*  contrast, so that test is not discriminating.  Ignore the SGEQ_P line.)
*
* Realistic story: the subgroup marker (industry code / insurance status /
* self-reported race) is recorded for everyone at baseline but stops being
* collected for the *high-churn* units, who are also the ones with the
* biggest treatment responses.  Classic non-ignorable item nonresponse.
*
* TRUTH BY CONSTRUCTION: sg1 unit effects are 1..30 scaled; sg0 unit effects
* are a constant 1.  The equal-unit contrast att_i(sg1) - att_i(sg0) has a
* known value.  Arm A observes sg everywhere; Arm B loses sg on the later
* post years of the TEN LARGEST-effect sg1 units.  Nothing about the outcome
* data changes between the arms -- only which rows carry a subgroup label.
* ==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 90
gen long id = _n
gen byte sg = cond(id <= 30, 1, 0)
gen double ttre = .
replace ttre = 5 + 2*mod(id-1,3) if id <= 60
* sg1 effects spread 0.5 .. 5.5 ; sg0 effects all 1
gen double eff = 0
replace eff = 0.5 + 5*(id-1)/29 if id <= 30
replace eff = 1                 if id > 30 & id <= 60
* the 10 LARGEST-effect sg1 units are the high-churn ones
gen byte churn = (id > 20 & id <= 30)
gen double u_i = rnormal(0,1)
expand 14
bysort id: gen int t = _n
gen double y = u_i + 0.3*t + eff*(t >= ttre & ttre < .) + rnormal(0,0.5)

* truth: equal-unit att_i in sg1 = mean(0.5..5.5) = 3 ; in sg0 = 1 ; contrast = 2
scalar TRUE_SG1 = 3
scalar TRUE_SG0 = 1
scalar TRUE_CON = 2

gen byte sgA = sg
gen byte sgB = sg
replace sgB = . if churn == 1 & ttre < . & t >= ttre + 2

qui count if sgB >= . 
dis "ROWS_WITH_MISSING_SG=" r(N)
qui count
dis "TOTAL_ROWS=" r(N)

foreach arm in sgA sgB {
  dis "############### ARM `arm' ###############"
  dis "CAUTIONSCAN_`arm'_START"
  capture noisily wooldid y id t ttre, att att_i subgroups(`arm') jointtests
  dis "ARM_`arm'_RC=" _rc
  dis "CAUTIONSCAN_`arm'_END"
  matrix R1_`arm' = e(wdidsgresults_1)
  matrix R0_`arm' = e(wdidsgresults_0)
  matrix PW_`arm' = e(pwsgs)
  dis "--- e(pwsgs) for `arm' ---"
  matrix list PW_`arm', format(%12.8f)
  local pwrn : rowfullnames PW_`arm'
  dis "PWROWNAMES_`arm': `pwrn'"
  local r1rn : rowfullnames R1_`arm'
  dis "R1ROWNAMES_`arm': `r1rn'"
  dis "SGEQ_P_`arm' = " e(joint_sg_att_i)
  ereturn list
}

dis "==================== SUMMARY ===================="
dis "TRUE equal-unit att_i(sg1) = " %12.6f TRUE_SG1
dis "TRUE equal-unit att_i(sg0) = " %12.6f TRUE_SG0
dis "TRUE contrast              = " %12.6f TRUE_CON
dis "A att(sg1)     = " %12.6f R1_sgA[1,1] "   B att(sg1)     = " %12.6f R1_sgB[1,1]
dis "A att_i(sg1)   = " %12.6f R1_sgA[2,1] "   B att_i(sg1)   = " %12.6f R1_sgB[2,1]
dis "A att_i(sg0)   = " %12.6f R0_sgA[2,1] "   B att_i(sg0)   = " %12.6f R0_sgB[2,1]
dis "A contrast rows:"
matrix list PW_sgA, format(%12.8f)
dis "B contrast rows:"
matrix list PW_sgB, format(%12.8f)
dis "ALLDONE"
