*=======================================================================
* r2_esnames.do -- SKEPTIC's independent repro of hv-sgname-02.
*
* CLAIM UNDER TEST: sthlp:1305-1306 promises that for the other estimand
* types you get "the same names with att replaced by att_it, att_i,
* att_itime or ame", i.e. ES_att_i_pre#/ES_att_itime_post#.  The claim
* is that no such column is ever posted: att_it/att_i/att_itime share
* ONE event-study block named ES_att_it, as sthlp:1315 says.
*
* Also checks (beyond the original prober): does the PRINTED table
* mislabel the event study when only att_i was requested?
*
* SELF-CONTAINED.  Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r2_esnames.do
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

cap prog drop dumpcols
prog def dumpcols
  args mname tag
  local cn : colnames `mname'
  local ce : coleq   `mname'
  dis "DUMPCOLS_`tag'_NCOL=" colsof(`mname')
  forval j = 1/`=colsof(`mname')' {
    local nj : word `j' of `cn'
    local ej : word `j' of `ce'
    dis "DC_`tag' `j' eq=<`ej'> nm=<`nj'>"
  }
end

*=================== OWN DATA ==========================================
* Deliberately UNEQUAL cohort sizes with heterogeneous effects, so the
* cell-weighted att and the equal-weighted att_i/att_it genuinely differ.
* onset 2004: 3 cohort-groups, 40 units each (120 units), true effect 0.4
* onset 2007: 2 cohort-groups, 6 units each  ( 12 units), true effect 2.5
* never-treated: 60 units in 2 groups.
clear
set obs 192
gen long unit = _n
gen int grp = .
replace grp = ceil(_n/40)            if _n <= 120        // groups 1,2,3
replace grp = 3 + ceil((_n-120)/6)   if _n > 120 & _n <= 132   // 4,5
replace grp = 5 + ceil((_n-132)/30)  if _n > 132        // 6,7 never
gen int ttre = .
replace ttre = 2004 if inrange(grp,1,3)
replace ttre = 2007 if inrange(grp,4,5)
gen double truetau = cond(inrange(grp,1,3), 0.4, cond(inrange(grp,4,5), 2.5, 0))
gen double ai = rnormal()
expand 10
bysort unit: gen int year = 2000 + _n - 1
gen byte post = (ttre < .) & (year >= ttre)
gen double y = ai + 0.05*(year-2000) + truetau*post + 0.20*rnormal()
compress
dis "DATA_N=" _N
tab grp if year==2000

*=================== RUN 1: att_i ALONE ================================
dis "===== RUN1 att_i alone ====="
cap noi wooldid y grp year ttre, att_i esprelength(2) espostlength(2) robust
dis "RUN1_RC=" _rc
tempname B1
matrix `B1' = e(b)
dumpcols `B1' R1
foreach nm in ES_att_i_post0 ES_att_i_pre1 ES_att_it_post0 ES_att_post0 {
  cap di _b[`nm']
  dis "R1_BREF `nm' RC=" _rc
}
foreach en in b_es_att_i b_es_att_it b_es_att numpre_att_i numpre_att_it {
  cap confirm matrix e(`en')
  local r1 = _rc
  cap confirm scalar e(`en')
  local r2 = _rc
  dis "R1_EOBJ `en' matRC=`r1' sclRC=`r2'"
}

*=================== RUN 2: att_itime ALONE ============================
dis "===== RUN2 att_itime alone ====="
cap noi wooldid y grp year ttre, att_itime esprelength(2) espostlength(2) robust
dis "RUN2_RC=" _rc
tempname B2
matrix `B2' = e(b)
dumpcols `B2' R2
foreach nm in ES_att_itime_post0 ES_att_it_post0 {
  cap di _b[`nm']
  dis "R2_BREF `nm' RC=" _rc
}

*=================== RUN 3: att + att_i together, NUMBERS ==============
dis "===== RUN3 att att_i together ====="
cap noi wooldid y grp year ttre, att att_i esprelength(2) espostlength(2) robust
dis "RUN3_RC=" _rc
tempname B3
matrix `B3' = e(b)
dumpcols `B3' R3
cap dis "R3_ES_att_post0   = " %21x _b[ES_att_post0]
cap dis "R3_ES_att_it_post0= " %21x _b[ES_att_it_post0]
cap dis "R3_ES_att_post0_d   = " _b[ES_att_post0]
cap dis "R3_ES_att_it_post0_d= " _b[ES_att_it_post0]
cap dis "R3_main_att   = " _b[main_att]
cap dis "R3_main_att_i = " _b[main_att_i]
* design truths for post0 (onset year): cell-weighted over units
*   (120*0.4 + 12*2.5)/132 = 0.5909...   ; equal-weighted over 5 groups
*   (3*0.4 + 2*2.5)/5 = 1.24
dis "R3_TRUTH_cellweighted=" (120*0.4 + 12*2.5)/132
dis "R3_TRUTH_equalweighted=" (3*0.4 + 2*2.5)/5

*=================== RUN 4: contreat -> is there ES_ame_i? =============
dis "===== RUN4 contreat ame_i ====="
gen double dose = 0
replace dose = 0.5 + 0.5*mod(unit,3) if ttre < .
cap noi wooldid y grp year ttre, att_i contreat(dose) contreatcontrols(simple) ///
    esprelength(2) espostlength(2) robust
dis "RUN4_RC=" _rc
cap dumpcols e(b) R4
tempname B4
cap matrix `B4' = e(b)
cap dumpcols `B4' R4

dis "R2_DONE"
