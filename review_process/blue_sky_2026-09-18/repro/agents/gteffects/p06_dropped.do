* p06_dropped.do -- (a) an omitted (collinear) treated cell and what it does to the
* documented weighted-mean identity; (b) fastnose e() surface; (c) permute blanking.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime

*======================================================================
* An administratively plausible panel: the CONTROL sample stops being
* collected after 2008, while the treated cohort keeps being followed.
* Years 2009-2012 then contain treated observations only, so each of
* those calendar-year fixed effects is collinear with the cohort-by-
* relative-period dummy that covers exactly those rows.
*======================================================================
clear
set obs 55
gen long unit = _n
gen byte cohp = cond(_n<=40, 1, 2)
gen double ttre = cond(cohp==2, 2004, .)
gen double ai = rnormal()
qui expand 13
bys unit: gen int year = 1999 + _n
qui drop if cohp == 1 & year > 2008
tempvar tsh
bys year: gen double `tsh' = rnormal() if _n==1
bys year: replace `tsh' = `tsh'[1]
gen double te = 0
qui replace te = 1.0 + 0.2*(year - ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`tsh' + te + rnormal()*0.25
tempfile drp
qui save `drp'

dis "--- G1: default safety ---"
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(3) espostlength(3)
dis "G1_RC=" _rc

dis "--- G2: safety(off) ---"
use `drp', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(3) espostlength(3) safety(off)
dis "G2_RC=" _rc
gtcheck , tag(G2_safetyoff) pre es
cap noi matrix list e(gteffects), format(%12.6g)
* what the estimation sample actually contains
qui gen byte esx = e(sample)
qui levelsof year if esx & cohp==2, local(yy)
dis "G2 treated years in sample: `yy'"
dis "G2 rows in table vs realized treated cell count:"
qui gegen byte tg = tag(cohp year) if esx & cohp==2 & year != 2003
qui count if tg == 1
dis "G2_EXPECTEDCELLS=" r(N)

dis "--- G3: fastnose e() surface (known_bydesign S7 item 11) ---"
use `drp', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(3) espostlength(3) safety(off) fastnose
dis "G3_RC=" _rc
foreach k in b_es_att numpre_att V_es_att b_es_att_it numpre_att_it {
  cap confirm matrix e(`k')
  local m1 = _rc
  cap confirm scalar e(`k')
  local m2 = _rc
  dis "G3_HAS_`k'=" cond(`m1'==0 | `m2'==0, 1, 0)
}
gtcheck , tag(G3_fastnose) pre es

dis "--- G4: permute() blanking ---"
mkpanel , n(60) t(9) onsets("2004 2006") ysd(0.4) seed(808)
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) permute(shuffle) permreps(20) cluster(unit)
dis "G4_RC=" _rc
dis "G4_perm_gt_blank=" "`e(perm_gt_blank)'"
gtcheck , tag(G4_permute) pre es
cap noi matrix list e(gteffects), format(%12.6g)

dis "DONE_P06"
