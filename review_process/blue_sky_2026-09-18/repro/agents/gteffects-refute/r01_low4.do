* r01: independent check of GTE-1 (customweightmultiplier), GTE-2 (fastnose e(b_es_*)),
*      GTE-4 (weight columns for unrequested estimands), GTE-5 (timetrends x se column).
* Self-contained: generates its own data, asserts the pinned build.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------- reporting helpers ----------
cap prog drop gtrep
prog def gtrep
  syntax , tag(string)
  tempname G
  cap matrix `G' = e(gteffects)
  if _rc != 0 {
    dis "GT|`tag'|NOMATRIX"
    exit
  }
  local cn : colnames `G'
  local nr = rowsof(`G')
  local out "GT|`tag'|rows=`nr'"
  foreach v in se p w_att w_att_it w_att_i w_att_itime {
    local j : list posof "`v'" in cn
    if `j' == 0 {
      local out "`out'|`v'=NOCOL"
    }
    else {
      mata: st_local("nm", strofreal(sum(st_matrix("`G'")[.,`j'] :< .)))
      local out "`out'|`v'nm=`nm'"
    }
  }
  dis "`out'"
end

cap prog drop cellmass
prog def cellmass
  syntax , tag(string) col(string)
  tempname G
  cap matrix `G' = e(gteffects)
  if _rc != 0 {
    dis "MASS|`tag'|`col'|NOMATRIX"
    exit
  }
  local cn : colnames `G'
  local j : list posof "`col'" in cn
  local jr : list posof "relyr" in cn
  if `j' == 0 {
    dis "MASS|`tag'|`col'|NOCOL"
    exit
  }
  mata: GG = st_matrix("`G'"); sel = select(GG, GG[.,`jr'] :>= 0); v = select(sel[.,`j'], sel[.,`j'] :< .); ///
        st_local("mn", strofreal(min(v),"%12.8f")); st_local("mx", strofreal(max(v),"%12.8f")); ///
        st_local("nn", strofreal(rows(v)))
  dis "MASS|`tag'|`col'|n=`nn'|min=`mn'|max=`mx'"
end

cap prog drop hasret
prog def hasret
  syntax , tag(string)
  local out "HAS|`tag'"
  foreach m in b_es_att V_es_att b_es_att_it V_es_att_it {
    cap matrix drop __zz
    cap matrix __zz = e(`m')
    local ok = (_rc == 0)
    local out "`out'|`m'=`ok'"
  }
  foreach s in numpre_att numpre_att_it {
    local ok2 = (e(`s') < .)
    local out "`out'|`s'=`ok2'"
  }
  cap matrix drop __zz
  dis "`out'"
end

* ---------- data: 90 units, 3 cohorts, 10 years ----------
clear
set obs 90
gen long unit = _n
gen byte coh = 1 + mod(_n-1,3)
gen double ttre = .
replace ttre = 2004 if coh==1
replace ttre = 2007 if coh==2
gen double ufe = rnormal()
expand 10
bysort unit: gen int t = 1999 + _n
gen double tfe = 0.10*(t-1999)
gen byte trt = (ttre<.) & (t>=ttre)
gen double rel = t - ttre
gen double te = cond(trt==1, 1 + 0.30*rel, 0)
gen double y = 2 + ufe + tfe + te + rnormal()*0.5
* multiplier CONSTANT WITHIN each (coh,t) cell -- so within-cell variation cannot be blamed
gen double pop = 1 + 0.5*coh + 0.20*(t-1999)
compress
save "r01_dat.dta", replace

dis "===== A: GTE-1 customweightmultiplier ====="
use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att att_it att_i att_itime gteffects
dis "A1_RC=" _rc
gtrep, tag(A1_nocwm)
cellmass, tag(A1_nocwm) col(w_att_it)
cellmass, tag(A1_nocwm) col(w_att_i)
cellmass, tag(A1_nocwm) col(w_att_itime)

use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att att_it att_i att_itime gteffects customweightmultiplier(pop)
dis "A2_RC=" _rc
gtrep, tag(A2_cwm)
cellmass, tag(A2_cwm) col(w_att_it)
cellmass, tag(A2_cwm) col(w_att_i)
cellmass, tag(A2_cwm) col(w_att_itime)

* mechanism prediction: cell weight mass of w_att_it == the cell's value of pop exactly
* (pop constant within (coh,t) => sum_cell w*v/W_it = v_cell)
tempname G
matrix `G' = e(gteffects)
local cn : colnames `G'
local ji : list posof "i" in cn
local jc : list posof "caltime" in cn
local jr : list posof "relyr" in cn
local jw : list posof "w_att_it" in cn
mata: GG = st_matrix("`G'"); ///
      sel = select(GG, GG[.,`jr'] :>= 0); ///
      maxdev = 0; ///
      for (r=1; r<=rows(sel); r++) { ///
        pv = 1 + 0.5*sel[r,`ji'] + 0.20*(sel[r,`jc'] - 1999); ///
        d = abs(sel[r,`jw'] - pv); ///
        if (d > maxdev) maxdev = d ; ///
      } ; ///
      st_local("md", strofreal(maxdev,"%12.4e"))
dis "A2_MECH|max|w_att_it_cellmass - pop_cell| = `md'"

dis "===== B: GTE-4 weight columns for unrequested estimands ====="
foreach est in att att_it att_i att_itime {
  use "r01_dat.dta", clear
  cap noi wooldid y coh t ttre, `est' gteffects
  dis "B_`est'_RC=" _rc
  gtrep, tag(B_only_`est')
}

dis "===== C: GTE-2 fastnose e(b_es_*) ====="
use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att att_it esprelength(3) espostlength(3) esfixedbaseperiod
dis "C1_RC=" _rc
hasret, tag(C1_analytic)

use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att att_it esprelength(3) espostlength(3) esfixedbaseperiod fastnose
dis "C2_RC=" _rc
hasret, tag(C2_fastnose)

dis "===== D: GTE-5 timetrends x VCE ====="
use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) timetrends gteffects
dis "D1_RC=" _rc
gtrep, tag(D1_tt_iid)

use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) timetrends gteffects robust
dis "D2_RC=" _rc
gtrep, tag(D2_tt_robust)

use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) timetrends gteffects cluster(unit)
dis "D3_RC=" _rc
gtrep, tag(D3_tt_clusterbelow)

use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) timetrends gteffects cluster(coh)
dis "D4_RC=" _rc
gtrep, tag(D4_tt_clusterati)

use "r01_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) gteffects
dis "D5_RC=" _rc
gtrep, tag(D5_notimetrends_iid)

dis "===== R01 DONE ====="
