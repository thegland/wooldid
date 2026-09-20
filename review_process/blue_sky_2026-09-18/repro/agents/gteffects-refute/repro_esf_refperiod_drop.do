* ============================================================================
* SELF-CONTAINED REPRO -- under esfixedbaseperiod, a treated unit is deleted from
* the estimation sample ENTIRELY, with no message, unless it OBSERVES the exact
* relative period named by esrelativeto() (default -1).
*
* One missing outcome year, at relyr == -1 only, removes a unit that still has
* five usable pre-periods and five post-periods. The pooled-reference branch keeps
* the same unit on the same data. Nothing in wooldid.sthlp states this rule
* (the esfixedbaseperiod and esrelativeto() paragraphs, sthlp:722-727, are silent;
* sthlp:73 says only "cohorts treated before the sample begins are dropped").
* Mechanism: wooldid.ado:855-861 -- refperiodobserved = max(_t == refperiod) by i,
* then `replace ifvar = 0 if refperiodobserved == 0 & ttre != .`  (unit-level).
*
* DGP: 12 treated firms. Firms 1-4 have a true effect of 4, firms 5-12 of 1.
*      TRUE ATT over all 12 treated = 2.0.  Firms 1-4 are the ones with the hole.
* Run:  stata-mp -e do repro_esf_refperiod_drop.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop rep
prog def rep
  syntax , tag(string)
  tempname M
  cap matrix `M' = e(wdidmainresults)
  if _rc != 0 {
    dis "RES|`tag'|NOMAT"
    exit
  }
  local rn : rownames `M'
  local eqs : roweq `M'
  forval r = 1/`=rowsof(`M')' {
    local nm : word `r' of `rn'
    local eqn : word `r' of `eqs'
    if "`eqn'" == "main" dis "RES|`tag'|att=" %10.6f `M'[`r',1] "  N=" e(N) "  N_dropped=" e(N_dropped)
  }
  tempname G
  cap matrix `G' = e(gteffects)
  if _rc == 0 {
    local cn : colnames `G'
    local ji : list posof "i" in cn
    mata: st_local("nu", strofreal(rows(uniqrows(st_matrix("`G'")[.,`ji']))))
    dis "RES|`tag'|treated_units_in_e(gteffects)=`nu'"
  }
end

clear
set obs 30
gen long firm = _n
gen double ttre = cond(firm <= 12, 2006, .)
gen double eff  = cond(firm <= 4, 4, cond(firm <= 12, 1, 0))
gen byte holey  = (firm <= 4)
gen double ufe  = rnormal()
expand 11
bysort firm: gen int year = 1999 + _n
gen double tfe = 0.10*(year-1999)
gen byte trt = (ttre < .) & (year >= ttre)
gen double y = 2 + ufe + tfe + trt*eff + rnormal()*0.25
* firms 1-4: outcome missing in 2005 == relyr -1, and in NO other year
gen double ycomplete = y
replace y = . if holey == 1 & year == 2005
tempfile dat
save `dat', replace

dis "########## A: complete outcome (ycomplete), esfixedbaseperiod, esrelativeto(-1) ##########"
use `dat', clear
cap noi wooldid ycomplete firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3) gteffects
dis "A_RC=" _rc
rep, tag(A_complete)

dis "########## B: ONE missing year at relyr -1 for 4 of 12 treated firms ##########"
use `dat', clear
cap noi wooldid y firm year ttre, att esfixedbaseperiod esprelength(3) espostlength(3) gteffects
dis "B_RC=" _rc
rep, tag(B_hole_at_ref)

dis "########## C: identical data, POOLED reference branch ##########"
use `dat', clear
cap noi wooldid y firm year ttre, att esprelength(3) espostlength(3) gteffects
dis "C_RC=" _rc
rep, tag(C_hole_pooled)

dis "########## D: identical data, esfixedbaseperiod esrelativeto(-2) ##########"
use `dat', clear
cap noi wooldid y firm year ttre, att esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(3) gteffects
dis "D_RC=" _rc
rep, tag(D_hole_rel2)

dis "########## VERDICT ##########"
dis "TRUE ATT over all 12 treated firms = 2.0."
dis "A (complete) and C (pooled) recover it. B drops firms 1-4 in full and reports"
dis "the ATT over the surviving 8 firms (~1.0), with no caution of any kind -- the"
dis "only trace is a larger 'obs dropped from initial sample' count (44, not 4)."
dis "===== REPRO DONE ====="
