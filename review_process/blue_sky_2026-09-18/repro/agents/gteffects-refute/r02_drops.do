* r02: (P1) independent repro of GTE-3 -- silent drop of treated units whose first observed
*           period is at/after their own onset, on the POOLED reference branch.
*      (P2) OWN LEAD: under esfixedbaseperiod a treated unit is dropped IN FULL unless it
*           OBSERVES the reference relative period exactly. One missing outcome year at
*           relyr == -1 deletes the whole unit, silently, although it has 5 usable pre-periods
*           and the pooled branch keeps it.
*      (P3) boundary: what happens when EVERY treated unit is removed by each rule.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

cap prog drop attrep
prog def attrep
  syntax , tag(string)
  tempname M
  cap matrix `M' = e(wdidmainresults)
  if _rc != 0 {
    dis "ATT|`tag'|NOMAT"
    exit
  }
  local rn : rownames `M'
  local eqs : roweq `M'
  local nr = rowsof(`M')
  forval r = 1/`nr' {
    local nm : word `r' of `rn'
    local eqn : word `r' of `eqs'
    if "`eqn'" == "main" dis "ATT|`tag'|`eqn':`nm'|b=" %12.6f `M'[`r',1] "|se=" %12.6f `M'[`r',2]
  }
  dis "NOBS|`tag'|N=" e(N) "|N_dropped=" e(N_dropped)
end

cap prog drop gti
prog def gti
  syntax , tag(string)
  tempname G
  cap matrix `G' = e(gteffects)
  if _rc != 0 {
    dis "GTI|`tag'|NOMAT"
    exit
  }
  local cn : colnames `G'
  local ji : list posof "i" in cn
  mata: v = st_matrix("`G'")[.,`ji']; u = uniqrows(v); st_local("nu", strofreal(rows(u))); st_local("nr2", strofreal(rows(v)))
  dis "GTI|`tag'|rows=`nr2'|distinct_i=`nu'"
end

* ===================== P1 =====================
* 30 firms, 2000-2010. 12 treated in 2004, 18 never treated.
* Firms 1-4 are "born" in 2004: their register rows start at 2004 (their own onset year).
* Heterogeneous effect: firms 1-4 have effect 4, firms 5-12 have effect 1.
* True ATT over all 12 treated = 2.0;  ATT over the 8 with pre-coverage = 1.0.
clear
set obs 30
gen long firm = _n
gen double ttre = cond(firm <= 12, 2004, .)
gen double eff  = cond(firm <= 4, 4, cond(firm <= 12, 1, 0))
gen byte lateborn = (firm <= 4)
gen double ufe = rnormal()
expand 11
bysort firm: gen int year = 1999 + _n
gen double tfe = 0.10*(year-1999)
gen byte trt = (ttre < .) & (year >= ttre)
gen double y = 2 + ufe + tfe + trt*eff + rnormal()*0.25
compress
save "r02_p1_full.dta", replace

dis "===== P1A: balanced register (control arm) ====="
use "r02_p1_full.dta", clear
cap noi wooldid y firm year ttre, att gteffects
dis "P1A_RC=" _rc
attrep, tag(P1A_balanced)
gti,    tag(P1A_balanced)

dis "===== P1B: firms 1-4 enter the register only in 2004 (== their onset) ====="
use "r02_p1_full.dta", clear
drop if lateborn == 1 & year < 2004
cap noi wooldid y firm year ttre, att gteffects
dis "P1B_RC=" _rc
attrep, tag(P1B_lateentry)
gti,    tag(P1B_lateentry)

dis "===== P1C: same, i = onset cohort rather than firm ====="
use "r02_p1_full.dta", clear
drop if lateborn == 1 & year < 2004
gen long cohid = cond(ttre < ., 1, 2)
cap noi wooldid y cohid year ttre, att gteffects
dis "P1C_RC=" _rc
attrep, tag(P1C_lateentry_cohort_i)
gti,    tag(P1C_lateentry_cohort_i)

* ===================== P2 =====================
* 30 firms, 2000-2010. 12 treated in 2006, 18 never treated. FULLY observed panel except:
* firms 1-4 have a MISSING OUTCOME in 2005 only  (relyr == -1). They still have
* 2000-2004 (relyr -6..-2) and 2006-2010 (relyr 0..4).
* Heterogeneous effect: firms 1-4 effect 4, firms 5-12 effect 1.  True ATT = 2.0.
clear
set obs 30
gen long firm = _n
gen double ttre = cond(firm <= 12, 2006, .)
gen double eff  = cond(firm <= 4, 4, cond(firm <= 12, 1, 0))
gen byte holey  = (firm <= 4)
gen double ufe = rnormal()
expand 11
bysort firm: gen int year = 1999 + _n
gen double tfe = 0.10*(year-1999)
gen byte trt = (ttre < .) & (year >= ttre)
gen double y = 2 + ufe + tfe + trt*eff + rnormal()*0.25
compress
save "r02_p2_full.dta", replace

dis "===== P2A: complete data, esfixedbaseperiod esrelativeto(-1)  [control] ====="
use "r02_p2_full.dta", clear
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3)
dis "P2A_RC=" _rc
attrep, tag(P2A_complete_esf)

dis "===== P2B: y missing in 2005 (relyr -1) for firms 1-4, esfixedbaseperiod ====="
use "r02_p2_full.dta", clear
replace y = . if holey == 1 & year == 2005
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3) gteffects
dis "P2B_RC=" _rc
attrep, tag(P2B_hole_esf)
gti,    tag(P2B_hole_esf)

dis "===== P2C: same hole, POOLED reference branch ====="
use "r02_p2_full.dta", clear
replace y = . if holey == 1 & year == 2005
cap noi wooldid y firm year ttre, att esprelength(3) espostlength(3) gteffects
dis "P2C_RC=" _rc
attrep, tag(P2C_hole_pooled)
gti,    tag(P2C_hole_pooled)

dis "===== P2D: same hole, esfixedbaseperiod esrelativeto(-2) ====="
use "r02_p2_full.dta", clear
replace y = . if holey == 1 & year == 2005
cap noi wooldid y firm year ttre, att esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(3)
dis "P2D_RC=" _rc
attrep, tag(P2D_hole_esf_rel2)

dis "===== P2E: hole at relyr -2 (2004) instead, esfixedbaseperiod esrelativeto(-1) ====="
use "r02_p2_full.dta", clear
replace y = . if holey == 1 & year == 2004
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3)
dis "P2E_RC=" _rc
attrep, tag(P2E_hole_at_rel2)

dis "===== P2F: hole at relyr -1, esfixedbaseperiod + timetrends ====="
use "r02_p2_full.dta", clear
replace y = . if holey == 1 & year == 2005
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3) timetrends
dis "P2F_RC=" _rc
attrep, tag(P2F_hole_esf_tt)

dis "===== P2G: hole at relyr -2 only, esfixedbaseperiod + timetrends (needs -1 AND -2) ====="
use "r02_p2_full.dta", clear
replace y = . if holey == 1 & year == 2004
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3) timetrends
dis "P2G_RC=" _rc
attrep, tag(P2G_hole_rel2_esf_tt)

* ===================== P3: boundary, every treated unit removed =====================
dis "===== P3A: pooled, EVERY treated firm enters at its own onset ====="
use "r02_p1_full.dta", clear
drop if ttre < . & year < ttre
cap noi wooldid y firm year ttre, att
dis "P3A_RC=" _rc

dis "===== P3B: esfixedbaseperiod, EVERY treated firm missing y at relyr -1 ====="
use "r02_p2_full.dta", clear
replace y = . if ttre < . & year == 2005
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3)
dis "P3B_RC=" _rc

dis "===== R02 DONE ====="
