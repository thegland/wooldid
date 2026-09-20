* r09_attrit.do -- SELF-CONTAINED.
* (A) BOOTX-4 second attempt: one unit per cohort (every treated state a distinct
*     onset) so each treated (cohort, relyr) cell is a singleton -- the geometry the
*     help says leaves the analytic covariance all-zero under timetrends + cluster(i).
* (B) NEW: how far below reps() can the number of draws BEHIND e(V_es_att) fall on a
*     realistic messy panel, and does anything in e() record it?  wooldid records
*     e(boot_completed) and e(boot_N_V) (listwise over ALL columns) but posts NO count
*     for the per-family V_es listwise pass that wooldid_ptp then consumes.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

dis "======== (A) singleton treated cells: timetrends x cluster(i)"
clear
set obs 26
gen state = _n
gen onset = .
replace onset = 5 + _n if _n <= 10          /* 10 treated states, ALL DISTINCT onsets */
expand 18
bysort state: gen t = _n
bysort state (t): gen double ueff = rnormal(0,.4) if _n==1
bysort state (t): replace ueff = ueff[_n-1] if _n>1
bysort state (t): gen double utr = rnormal(0,.05) if _n==1
bysort state (t): replace utr = utr[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .1*rel, 0)
gen double y = 2 + .05*t + ueff + utr*t + te + e
drop rel te
tempfile SNG
save `SNG'
cap noi wooldid y state t onset, cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) timetrends
dis "SNG_A_RC=" _rc
cap confirm matrix e(V)
dis "SNG_A|vcvstate=`e(wdid_vcvstate)'|confirmV=" _rc
cap noi wooldid y state t onset, cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) timetrends bootstrap(standard) reps(30) rseed(4242)
dis "SNG_B_RC=" _rc
cap confirm matrix e(V)
local r1 = _rc
cap confirm matrix e(V_boot)
dis "SNG_B|vcvstate=`e(wdid_vcvstate)'|confirmV=`r1'|confirmVboot=" _rc "|vposted=" e(boot_v_posted)
cap noi test ES_att_post0 = 0
dis "SNG_B_handtest_rc=" _rc

dis "======== (B) listwise attrition behind e(V_es_att)"
clear
set obs 40
gen state = _n
gen onset = .
replace onset = 12 + mod(_n,4) if inlist(_n, 3, 9, 14, 21, 27, 33)   /* 6 treated states, late */
expand 22
bysort state: gen t = _n
bysort state (t): gen double ueff = rnormal(0,.4) if _n==1
bysort state (t): replace ueff = ueff[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .6 + .08*rel, 0)
gen double y = 2 + .04*t + ueff + te + e
drop rel te
tab onset, m

local OPT "cluster(state) esfixedbaseperiod esprelength(9) espostlength(6) jointtests pretrendpower(80)"
cap noi wooldid y state t onset, `OPT'
dis "B_A_RC=" _rc
dis "BRES|ANALYTIC|keff=" e(roth_keff_att) "|npre=" e(roth_npre_att) "|prep=" %12.8g e(roth_prep_att) "|slope=" %12.8g e(roth_slope_att) "|se=" %12.8g e(roth_slope_se_att) "|mdt=" %12.8g e(roth_power_att)[1,2]
foreach R in 200 {
  cap noi wooldid y state t onset, `OPT' bootstrap(standard) reps(`R') rseed(313)
  dis "B_`R'_RC=" _rc
  dis "BRES|BOOT`R'|keff=" e(roth_keff_att) "|npre=" e(roth_npre_att) "|prep=" %12.8g e(roth_prep_att) "|slope=" %12.8g e(roth_slope_att) "|se=" %12.8g e(roth_slope_se_att) "|mdt=" %12.8g e(roth_power_att)[1,2] "|completed=" e(boot_completed) "|NV=" e(boot_N_V) "|rankjoint=" %12.8g e(joint_horiz)
  * the count wooldid never posts: listwise-complete replications over the ES_att columns
  cap confirm matrix e(bootdist)
  if _rc==0 {
    matrix BD = e(bootdist)
    mata:
      BD = st_matrix("BD")
      nm = st_matrixcolstripe("BD")[.,2]
      keep = J(0,1,.)
      for (b=1; b<=rows(nm); b++) if (substr(nm[b],1,7)=="ES_att_") keep = keep \ b
      C = BD[., keep]
      ok = rowmissing(C) :== 0
      printf("VESN| reps rows = %f ; ES_att columns = %f ; LISTWISE COMPLETE = %f\n", rows(BD), rows(keep), sum(ok))
      // pre columns only
      keep2 = J(0,1,.)
      for (b=1; b<=rows(nm); b++) if (substr(nm[b],1,11)=="ES_att_pre") keep2 = keep2 \ b
      if (rows(keep2)>0) {
        C2 = BD[., keep2]
        printf("VESN| pre columns = %f ; listwise complete over PRE only = %f\n", rows(keep2), sum(rowmissing(C2):==0))
      }
      // per-column availability
      for (b=1; b<=rows(keep); b++) printf("   col %s available %f\n", nm[keep[b]], sum(BD[.,keep[b]] :< .))
    end
  }
}
dis "DONE"
* NOTE: part (B) of this file aborts (--Break--) because its `mata: ... end' block sits
* inside a `foreach { }' loop.  Part (A) - the BOOTX-4 singleton-cell attempt - runs and
* is the part cited in FINDINGS.md.  Part (B) was rerun at top level as r10b_attrit.do.
