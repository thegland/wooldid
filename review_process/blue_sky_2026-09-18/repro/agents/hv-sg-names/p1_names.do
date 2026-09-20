*==========================================================================
* p1_names.do -- SELF-CONTAINED enumeration of wooldid's e(b) column names
*   across estimands x subgroups x families.  Generates its own data.
*   Verifies (a) whether subgroup rows take the help's "_s#" SUFFIX,
*            (b) whether ES_att_i_* / ES_att_itime_* exist at all.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------- name dumper ----------
cap prog drop dumpn
prog def dumpn
  args tag
  local nf : colfullnames e(b)
  local nn : colnames e(b)
  local k = 0
  foreach n of local nn {
    local k = `k' + 1
    di "NAME|`tag'|`k'|`n'"
  }
  di "NCOL|`tag'|`k'"
  local j = 0
  foreach n of local nf {
    local j = `j' + 1
    di "FULL|`tag'|`j'|`n'"
  }
end

cap prog drop dumpmats
prog def dumpmats
  args tag
  local M : e(matrices)
  foreach m of local M {
    di "MAT|`tag'|`m'"
  }
  local S : e(macros)
  foreach m of local S {
    di "EMAC|`tag'|`m'"
  }
end

* ---------- DGP: plain staggered panel, unit + time FE, 3 treated cohorts ----------
set obs 80
gen long unit = _n
gen double u_fe = rnormal()
gen byte sgA = 1 + mod(_n, 2)          // subgroup values 1,2
gen byte sgB = cond(mod(_n,2)==0, 3, 17)  // subgroup values 3,17 (non-1..k)
gen byte cohcode = mod(_n, 4)
gen int ttre = .
replace ttre = 2003 if cohcode == 1
replace ttre = 2005 if cohcode == 2
replace ttre = 2007 if cohcode == 3
* cohcode == 0 -> never treated
expand 10
bysort unit: gen int year = 1999 + _n
gen double t_fe = 0.2*(year - 2000)
gen double post = (year >= ttre) & ttre < .
gen double tau  = 1.0 + 0.3*(year - ttre) if post
replace tau = 0 if !post
gen double y = u_fe + t_fe + tau + 0.30*rnormal()

save "p1_data.dta", replace

* ==================== RUN A: one estimand at a time, no subgroups ====
foreach est in att att_i att_it att_itime {
  di "=== RUN A_`est' ==="
  cap noi wooldid y unit year ttre, `est' esprelength(2) espostlength(2) robust
  di "A_`est'_RC=" _rc
  if _rc == 0 {
    dumpn "A_`est'"
    dumpmats "A_`est'"
  }
}

* ==================== RUN B: all four estimands at once ==============
di "=== RUN B_all ==="
cap noi wooldid y unit year ttre, att att_i att_it att_itime esprelength(2) espostlength(2) robust
di "B_all_RC=" _rc
if _rc == 0 {
  dumpn "B_all"
  dumpmats "B_all"
}

* ==================== RUN C: all four + subgroups(sgA) 1/2 ===========
di "=== RUN C_sgA ==="
cap noi wooldid y unit year ttre, att att_i att_it att_itime subgroups(sgA) esprelength(2) espostlength(2) robust
di "C_sgA_RC=" _rc
if _rc == 0 {
  dumpn "C_sgA"
  dumpmats "C_sgA"
}

* ==================== RUN D: subgroups(sgB) values 3/17 ==============
di "=== RUN D_sgB ==="
cap noi wooldid y unit year ttre, att att_i att_it att_itime subgroups(sgB) esprelength(2) espostlength(2) robust
di "D_sgB_RC=" _rc
if _rc == 0 {
  dumpn "D_sgB"
  dumpmats "D_sgB"
}

* ==================== RUN E: continuous treatment -> ame family ======
di "=== RUN E_ame ==="
gen double dose = abs(u_fe) + 0.5
replace dose = 0 if ttre >= .
cap noi wooldid y unit year ttre, att att_it contreat(dose) esprelength(2) espostlength(2) robust
di "E_ame_RC=" _rc
if _rc == 0 {
  dumpn "E_ame"
  dumpmats "E_ame"
}

* ==================== RUN F: ame + subgroups =========================
di "=== RUN F_amesg ==="
cap noi wooldid y unit year ttre, att att_it contreat(dose) subgroups(sgA) esprelength(2) espostlength(2) robust
di "F_amesg_RC=" _rc
if _rc == 0 {
  dumpn "F_amesg"
  dumpmats "F_amesg"
}

* ==================== RUN C2: explicit confirm of the THREE naming schemes ====
di "=== RUN C2 confirm ==="
cap noi wooldid y unit year ttre, att att_i att_it att_itime subgroups(sgA) esprelength(2) espostlength(2) robust
di "C2_RC=" _rc
foreach cand in b_es_att b_es_att_s1 b_es_att_SG1 b_es_att_s2 b_es_att_it b_es_att_it_s1 ///
                V_es_att_s1 V_es_att_SG1 wdidsgresults_1 wdidsgresults_s1 wdidsgresults_SG1 ///
                wdidmainresults pwsgs {
  cap confirm matrix e(`cand')
  di "CONFIRM|`cand'|rc=" _rc
}
* rownames of the subgroup results matrix (help: equations NOT suffixed there)
cap noi mat list e(wdidsgresults_1), noheader format(%9.4f)
di "SGRESRC=" _rc

* ---- CONSEQUENCE LEG: address a coefficient by the name the help gives ----
di "--- help name (suffix _s#) ---"
cap noi test main_att_s1 = main_att_s2
di "HELPNAME_TEST_RC=" _rc
cap noi lincom main_att_s1
di "HELPNAME_LINCOM_RC=" _rc
di "--- actual name (infix _SG#_) ---"
cap noi test main_SG1_att = main_SG2_att
di "TRUENAME_TEST_RC=" _rc
cap noi lincom main_SG1_att
di "TRUENAME_LINCOM_RC=" _rc
di "--- ES names the help implies for att_i / att_itime ---"
cap noi lincom ES_att_i_post0
di "ES_ATT_I_RC=" _rc
cap noi lincom ES_att_itime_post0
di "ES_ATT_ITIME_RC=" _rc
cap noi lincom ES_att_it_post0
di "ES_ATT_IT_RC=" _rc
di "--- esttab / coefplot with keep() written from the help ---"
cap noi esttab ., keep(*_s1*) nostar
di "ESTTAB_HELPKEEP_RC=" _rc
cap noi esttab ., keep(*SG1*) nostar
di "ESTTAB_TRUEKEEP_RC=" _rc
cap noi coefplot, keep(*_s1*) nodraw
di "COEFPLOT_HELPKEEP_RC=" _rc

di "PROBE1_DONE"
