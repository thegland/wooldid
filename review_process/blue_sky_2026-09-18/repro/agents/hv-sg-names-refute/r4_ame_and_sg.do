*=======================================================================
* r4_ame_and_sg.do -- fixes r3's three failed legs.
*   (B) contreat(): sthlp:239 says "Under contreat() every requested
*       estimand is reported twice ... (ame, ame_it, ame_i, ame_itime).
*       Both families are displayed, and both are stored under the same
*       names with att replaced by ame."  Test it: request att_i with
*       contreat() and look for main_ame_i and ES_ame_i_post0.
*   (C) subgroups() + att_i only, with sg constant within i.
*   (D) bootstrap table COLUMN names: does the _SG infix hold there?
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
gen double dose = cond(ttre<., 0.5 + 0.5*mod(unit,3), 0)
gen double y = ai + 0.05*(year-2000) + truetau*post + 0.3*dose*post + 0.20*rnormal()
* subgroup CONSTANT within i = grp (wooldid requires this)
gen byte sg = 1 + mod(grp,2)
compress

*===================== (B) contreat + att_i ============================
dis "===== B: contreat(dose), att_i requested ====="
cap noi wooldid y grp year ttre, att_i contreat(dose) ///
    contreatcontrols(basic) esprelength(2) espostlength(2) robust
noi dis "B_RC=" _rc
if _rc == 0 {
  tempname Bb
  matrix `Bb' = e(b)
  dumpcolsp `Bb' B
  foreach nm in main_ame_i main_ame main_ame_it main_ame_itime ///
                ES_ame_i_post0 ES_ame_itime_post0 ES_ame_it_post0 ES_ame_post0 {
    cap di _b[`nm']
    noi dis "B_BREF `nm' RC=" _rc
  }
  foreach en in b_es_ame_i b_es_ame_it b_es_ame b_es_att_it {
    cap confirm matrix e(`en')
    noi dis "B_EMAT `en' RC=" _rc
  }
}

*===================== (B2) contreat + att (baseline) ==================
dis "===== B2: contreat(dose), att requested ====="
cap noi wooldid y grp year ttre, att contreat(dose) ///
    contreatcontrols(basic) esprelength(2) espostlength(2) robust
noi dis "B2_RC=" _rc
tempname Bb2
cap matrix `Bb2' = e(b)
cap noi dumpcolsp `Bb2' B2

*===================== (C) subgroups + att_i only ======================
dis "===== C: subgroups(sg), att_i only ====="
cap noi wooldid y grp year ttre, att_i subgroups(sg) esprelength(2) ///
    espostlength(2) robust
noi dis "C_RC=" _rc
if _rc == 0 {
  tempname Bc
  matrix `Bc' = e(b)
  dumpcolsp `Bc' C
  foreach nm in ES_att_i_SG1_post0 ES_att_it_SG1_post0 ES_att_i_post0_s1 ///
                main_SG1_att_i main_att_i_s1 {
    cap di _b[`nm']
    noi dis "C_BREF `nm' RC=" _rc
  }
}

*===================== (D) bootstrap column names ======================
dis "===== D: bootstrap(standard) + subgroups ====="
cap noi wooldid y grp year ttre, att subgroups(sg) esprelength(2) ///
    espostlength(2) robust bootstrap(standard) reps(8) rseed(4321)
noi dis "D_RC=" _rc
if _rc == 0 {
  foreach mm in boot_wald_mainresults bootdist V_boot {
    cap matrix BT = e(`mm')
    if _rc == 0 {
      noi dumpcolsp BT D_`mm'
      matrix drop BT
    }
    else noi dis "D_ABSENT `mm'"
  }
}

dis "R4_DONE"
