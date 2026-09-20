*=======================================================================
* r3_sharpen.do -- third leg.  Closes the loops r1/r2 left open:
*   (A) exact hex: is the shared ES block bit-identical across
*       att_it / att_i / att_itime, and does it equal the equal-weighted
*       design truth rather than the cell-weighted one?
*   (B) contreat(): is there an ES_ame_i / ES_ame_itime?  (sthlp:1306's
*       "or ame" promise, and sthlp:239's "both families are displayed,
*       and both are stored under the same names with att replaced by
*       ame")
*   (C) subgroups() + att_i only: what does the SG event-study block get
*       called in e(b) and in the PRINTED table?
*   (D) bootstrap table COLUMN names: does the _SG infix hold there too?
* SELF-CONTAINED.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

prog def dumpcolsp
  args mname tag
  local cn : colnames `mname'
  dis "DUMPCOLS_`tag'_NCOL=" colsof(`mname')
  forval j = 1/`=colsof(`mname')' {
    local nj : word `j' of `cn'
    dis "DC_`tag' `j' nm=<`nj'>"
  }
end

*--- same design as r2: 5 treated groups, grossly unequal sizes --------
clear
set obs 192
gen long unit = _n
gen int grp = .
replace grp = ceil(_n/40)            if _n <= 120
replace grp = 3 + ceil((_n-120)/6)   if _n > 120 & _n <= 132
replace grp = 5 + ceil((_n-132)/30)  if _n > 132
gen int ttre = .
replace ttre = 2004 if inrange(grp,1,3)
replace ttre = 2007 if inrange(grp,4,5)
gen double truetau = cond(inrange(grp,1,3), 0.4, cond(inrange(grp,4,5), 2.5, 0))
gen double ai = rnormal()
expand 10
bysort unit: gen int year = 2000 + _n - 1
gen byte post = (ttre < .) & (year >= ttre)
gen double y = ai + 0.05*(year-2000) + truetau*post + 0.20*rnormal()
gen byte sg = 1 + mod(unit,2)
gen double dose = cond(ttre<., 0.5 + 0.5*mod(unit,3), 0)
compress

*===================== (A) bit-identical shared ES block ===============
foreach est in att_it att_i att_itime {
  qui wooldid y grp year ttre, `est' esprelength(2) espostlength(2) robust
  noi dis "A_`est'_ES_att_it_post0 = " %21x _b[ES_att_it_post0]
  noi dis "A_`est'_ES_att_it_post0_dec = " _b[ES_att_it_post0]
  noi dis "A_`est'_main = " _b[main_`est']
}
qui wooldid y grp year ttre, att esprelength(2) espostlength(2) robust
noi dis "A_att_ES_att_post0 = " %21x _b[ES_att_post0]
noi dis "A_att_ES_att_post0_dec = " _b[ES_att_post0]
noi dis "A_TRUTH_cellweighted = " (120*0.4 + 12*2.5)/132
noi dis "A_TRUTH_equalweighted = " (3*0.4 + 2*2.5)/5

*===================== (B) contreat: ES_ame_i? =========================
cap noi wooldid y grp year ttre, att_i ame_i contreat(dose) ///
    contreatcontrols(basic) esprelength(2) espostlength(2) robust
noi dis "B_RC=" _rc
if _rc == 0 {
  tempname Bb
  matrix `Bb' = e(b)
  dumpcolsp `Bb' B
  foreach nm in ES_ame_i_post0 ES_ame_itime_post0 ES_ame_it_post0 ES_ame_post0 {
    cap di _b[`nm']
    noi dis "B_BREF `nm' RC=" _rc
  }
  foreach en in b_es_ame_i b_es_ame_it b_es_ame {
    cap confirm matrix e(`en')
    noi dis "B_EMAT `en' RC=" _rc
  }
}

*===================== (C) subgroups + att_i only ======================
cap noi wooldid y grp year ttre, att_i subgroups(sg) esprelength(2) ///
    espostlength(2) robust
noi dis "C_RC=" _rc
tempname Bc
cap matrix `Bc' = e(b)
cap noi dumpcolsp `Bc' C
foreach nm in ES_att_i_SG1_post0 ES_att_it_SG1_post0 ES_att_i_post0_s1 {
  cap di _b[`nm']
  noi dis "C_BREF `nm' RC=" _rc
}

*===================== (D) bootstrap column names ======================
cap noi wooldid y grp year ttre, att subgroups(sg) esprelength(2) ///
    espostlength(2) robust bootstrap(standard) reps(8) rseed(4321)
noi dis "D_RC=" _rc
foreach mm in boot_wald_mainresults bootdist V_boot {
  cap matrix BT = e(`mm')
  if _rc == 0 {
    noi dumpcolsp BT D_`mm'
    matrix drop BT
  }
  else noi dis "D_ABSENT `mm'"
}

dis "R3_DONE"
