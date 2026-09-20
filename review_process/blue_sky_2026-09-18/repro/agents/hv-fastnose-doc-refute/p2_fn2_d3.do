* =========================================================================
* P2 -- tests HV-FN2 / HV-FN3 on MY OWN data, three arms on the same panel:
*   A1  analytic  robust
*   A2  robust fastnose
*   A3  robust safety(off)
* Design: 40 states x 16 years x 20 individuals; dose = state generosity,
* CONSTANT WITHIN STATE (= within i), so every treated cohort-period
* dose-slope term drops for multicollinearity (the help's own named case,
* sthlp:902).  Structural dose effect in the DGP is 0.40, not zero.
* =========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- inline data generation (self-contained) ----
set obs 12800
gen long row = _n
gen int  st  = ceil(row/320)
gen int  yr  = 2005 + mod(ceil(row/20)-1, 16)
gen int adopt = .
replace adopt = 2010 if inrange(st, 1, 10)
replace adopt = 2013 if inrange(st, 11, 20)
replace adopt = 2016 if inrange(st, 21, 28)
bysort st: gen double u_st = runiform() if _n == 1
bysort st (u_st): replace u_st = u_st[1]
sort row
gen double dose = round(0.5 + 2.5*u_st, 0.01)
gen byte D = (yr >= adopt) & !missing(adopt)
bysort st: gen double afe = rnormal() if _n == 1
bysort st (afe): replace afe = afe[1]
sort row
gen double y = 2 + afe + 0.03*(yr-2005) + 0.40*dose*D + rnormal(0,1)
dis "DOSE_RANGE:"
sum dose
* confirm the geometry: dose has zero within-state variance
bysort st: egen double sdw = sd(dose)
sum sdw
dis "MAX_WITHIN_STATE_SD_OF_DOSE=" r(max)
drop sdw

prog def rep2
  args tag
  cap confirm matrix e(b_es_ame)
  dis "`tag'_BESAME_RC=" _rc
  cap confirm matrix e(V_es_ame)
  dis "`tag'_VESAME_RC=" _rc
  cap confirm matrix e(b_es_att)
  dis "`tag'_BESATT_RC=" _rc
  local npa = e(numpre_ame)
  dis "`tag'_NUMPRE_AME=[" "`npa'" "]"
  local npt = e(numpre_att)
  dis "`tag'_NUMPRE_ATT=[" "`npt'" "]"
  dis "`tag'_VCVSTATE=[" e(wdid_vcvstate) "]"
  cap noi mat li e(b_es_ame), format(%21x)
  dis "`tag'_LISTBESAME_RC=" _rc
  dis "--- `tag' wdidmainresults (col 1 = estimate) ---"
  cap noi mat li e(wdidmainresults), format(%12.8f)
  dis "--- `tag' raw e(b) ---"
  cap noi mat li e(b), format(%12.8f)
end

dis "=== ARM 1: analytic robust ==="
cap noi wooldid y st yr adopt, att robust contreat(dose) contreatcontrols(basic) esprelength(3) espostlength(4) esfixedbaseperiod
dis "A1_RC=" _rc
rep2 A1

dis "=== ARM 2: robust fastnose ==="
cap noi wooldid y st yr adopt, att robust contreat(dose) contreatcontrols(basic) esprelength(3) espostlength(4) esfixedbaseperiod fastnose
dis "A2_RC=" _rc
rep2 A2

dis "=== ARM 3: robust safety(off) ==="
cap noi wooldid y st yr adopt, att robust contreat(dose) contreatcontrols(basic) esprelength(3) espostlength(4) esfixedbaseperiod safety(off)
dis "A3_RC=" _rc
rep2 A3

dis "P2_DONE"
