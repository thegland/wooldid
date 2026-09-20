*==========================================================================
* p4_ame_and_silent.do -- (1) ame-family names under subgroups() on an
*   IDENTIFIED continuous treatment; (2) which downstream consumers fail
*   SILENTLY when fed the help's "_s#" names.  SELF-CONTAINED.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop dumpn
prog def dumpn
  args tag
  local nn : colnames e(b)
  local k = 0
  foreach n of local nn {
    local k = `k' + 1
    di "NAME|`tag'|`k'|`n'"
  }
end
cap prog drop dumpinv
prog def dumpinv
  args tag
  local M : e(matrices)
  foreach m of local M {
    di "INV|`tag'|`m'"
  }
end

* ---------------- DGP: i = COHORT-LEVEL id, many units per cell ----------
* (a coarser i is exactly the design wooldid's help recommends; it also
*  gives each treated cohort-period cell real within-cell dose variation)
set obs 400
gen long unit = _n
gen byte grp  = mod(_n, 20) + 1          // 20 cohort-level ids, 20 units each
gen byte sgA  = 1 + mod(grp, 2)          // subgroup constant within grp
gen int ttre = .
replace ttre = 2003 if inrange(grp, 1, 5)
replace ttre = 2005 if inrange(grp, 6, 10)
replace ttre = 2007 if inrange(grp, 11, 15)
gen double g_fe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double dose = 1 + runiform()          // varies by unit AND year: identified
gen double post = (year >= ttre) & ttre < .
gen double y = g_fe + 0.2*(year-2000) + post*(0.8*dose) + 0.3*rnormal()
save "p4_data.dta", replace

* ============ P: contreat + subgroups, identified dose ============
di "=== RUN P_ame ==="
cap noi wooldid y grp year ttre, att att_it contreat(dose) contreatcontrols(basic) ///
    subgroups(sgA) esprelength(2) espostlength(2) robust
di "P_RC=" _rc
if _rc == 0 {
  dumpn "P_ame"
  dumpinv "P_ame"
}

* ============ Q: which consumers fail SILENTLY on the help's names? =======
di "=== RUN Q_silent ==="
use "p4_data.dta", clear
cap noi wooldid y grp year ttre, att subgroups(sgA) esprelength(2) espostlength(2) robust
di "Q_RC=" _rc

di "--- Q1: esttab keep() with the HELP name (expect loud r(111)) ---"
cap noi esttab ., keep(*_s1*) nostar
di "Q1_RC=" _rc

di "--- Q2: esttab DROP() with the HELP name: does it silently drop nothing? ---"
cap noi esttab ., drop(*_s1*) nostar
di "Q2_RC=" _rc

di "--- Q3: coefplot keep() with the HELP name (silent empty plot?) ---"
cap noi coefplot, keep(*_s1*) nodraw
di "Q3_RC=" _rc

di "--- Q4: coefplot drop() with the HELP name ---"
cap noi coefplot, drop(*_s1*) nodraw
di "Q4_RC=" _rc

di "--- Q5: event_plot with the help-derived subgroup stub ---"
cap noi event_plot , default_look stub_lag(ES_att_post#_s1) stub_lead(ES_att_pre#_s1) graph_opt(nodraw)
di "Q5_RC=" _rc
di "--- Q6: event_plot with the TRUE subgroup stub ---"
cap noi event_plot , default_look stub_lag(ES_att_SG1_post#) stub_lead(ES_att_SG1_pre#) graph_opt(nodraw)
di "Q6_RC=" _rc

di "--- Q7: the help's headline recipe, no subgroups (control) ---"
cap noi event_plot , default_look stub_lag(ES_att_post#) stub_lead(ES_att_pre#) graph_opt(nodraw)
di "Q7_RC=" _rc

di "--- Q8: _b[] / _se[] by the help name ---"
cap noi di _b[main_att_s1]
di "Q8_RC=" _rc
cap noi di _b[main_SG1_att]
di "Q8b_RC=" _rc

di "PROBE4_DONE"
