*========================================================================
* fpsg-01 THIRD LENS: independent reproduction, built from the claim alone.
* Question: with subgroups() partially missing on a unit's rows, is the
* reported subgroup att_i / att_itime the equal-unit average it is
* documented to be, the share-of-own-post-weight average the claim
* predicts, or the row-weighted average (innocent alternative)?
*========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

*---------------- helper: pull one row of a wooldid result matrix -------
cap prog drop l0get
prog def l0get, rclass
  syntax , MAT(string) EQ(string) NM(string)
  tempname M
  cap matrix `M' = `mat'
  if _rc != 0 {
    return scalar val = .
    exit
  }
  local rn : rownames `M'
  local re : roweq `M'
  local k = rowsof(`M')
  local got = .
  forval j = 1/`k' {
    local tn : word `j' of `rn'
    local te : word `j' of `re'
    if "`tn'" == "`nm'" & "`te'" == "`eq'" local got = `M'[`j',1]
  }
  return scalar val = `got'
end

*================= DATA =================================================
* 45 firms x 12 years (2000-2011), administrative shape.
* 6 treated cohorts, onsets 2004..2009, 5 firms each, cohort-constant
* true effects 6,5,4,3,2,1.  15 never-treated firms.
* Exactly ONE firm per cohort carries sector==1, so each sector-1 firm's
* own effect is separately identified and the six are all different.
clear
set obs 45
gen long id = _n
gen byte cohortno = ceil(id/5) if id <= 30
gen int ttre = 2003 + cohortno if cohortno < .
gen double effect = 7 - cohortno if cohortno < .
replace effect = 0 if effect >= .
gen double unitfe = 0.37*id - 3
expand 12
bysort id: gen int year = 1999 + _n
gen double yearfe = 0.21*(year-2000)
gen byte post = (ttre < .) & (year >= ttre)
gen double ynoiseless = unitfe + yearfe + effect*post

* sector: 1 for ids 1,6,11,16,21,26 (one per cohort); 0 for the other
* treated firms; alternating for the never-treated (they never enter the
* post subpopulation, so this is cosmetic realism only).
gen byte sector = 0
replace sector = 1 if inlist(id,1,6,11,16,21,26)
replace sector = mod(id,2) if id > 30

* ITEM NONRESPONSE: the sector field stops being recorded for two firms.
* firm 1  (onset 2004, P=8 post rows) keeps it for 2004 only   -> L=1
* firm 6  (onset 2005, P=7 post rows) keeps it for 2005-2006   -> L=2
* firms 11,16,21,26 keep it throughout                         -> L=P
gen byte blank = (id==1 & year>=2005) | (id==6 & year>=2007)
gen byte sector_miss = sector
replace sector_miss = . if blank
* the adjacent case wooldid REFUSES: same rows, value CHANGES 1->0
gen byte sector_chg = sector
replace sector_chg = 0 if blank

* per-firm subgroup ids, for the noisy-arm oracle (single-firm subgroups,
* where the share-vs-equal-unit distinction vanishes by construction)
gen byte sgu_full = .
gen byte sgu_lab  = .
local k = 0
foreach u of numlist 1 6 11 16 21 26 {
  local ++k
  replace sgu_full = `k' if id==`u'
  replace sgu_lab  = `k' if id==`u' & blank==0
}

* noisy twin of the same design
gen double ynoisy = ynoiseless + rnormal(0,0.3)

* sanity: post-row counts and labelled post-row counts for the 6 sector-1 firms
preserve
  keep if sector==1 & post==1
  gen byte one = 1
  gen byte lab = 1 - blank
  collapse (sum) P=one Llab=lab (mean) eff=effect, by(id)
  gen double share = Llab/P
  list, noobs
restore

compress
tempfile base
save "`base'", replace

*================= PREDICTIONS (arithmetic, stated before the runs) =====
* sector==1 firms: effects 6,5,4,3,2,1 ; P = 8,7,6,5,4,3 ; L = 1,2,6,5,4,3
scalar EQUALUNIT = (6+5+4+3+2+1)/6
scalar SHAREWTD  = (6/8 + 5*2/7 + 4 + 3 + 2 + 1) / (1/8 + 2/7 + 1 + 1 + 1 + 1)
scalar ROWWTD    = (1*6 + 2*5 + 6*4 + 5*3 + 4*2 + 3*1) / (1+2+6+5+4+3)
scalar ATTCLEAN  = (8*6 + 7*5 + 6*4 + 5*3 + 4*2 + 3*1) / (8+7+6+5+4+3)
dis "PRED_EQUALUNIT = " %20.10f EQUALUNIT
dis "PRED_SHAREWTD  = " %20.10f SHAREWTD
dis "PRED_ROWWTD    = " %20.10f ROWWTD
dis "PRED_ATTCLEAN  = " %20.10f ATTCLEAN

*================= ARM A: clean sector, noiseless =======================
dis "===ARM_A_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoiseless id year ttre, att att_i att_it att_itime subgroups(sector)
dis "ARM_A_RC=" _rc
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_i)
scalar A_sg1_att_i = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_itime)
scalar A_sg1_att_itime = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att)
scalar A_sg1_att = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_it)
scalar A_sg1_att_it = r(val)
l0get, mat(e(wdidsgresults_0)) eq(main) nm(att_i)
scalar A_sg0_att_i = r(val)
l0get, mat(e(wdidmainresults)) eq(main) nm(att_i)
scalar A_all_att_i = r(val)
scalar A_N = e(N)
scalar A_Ndrop = e(N_dropped)
cap matrix GTA = e(gteffects)
cap matrix PWA = e(pwsgs)
dis "A rownames of wdidsgresults_1:"
cap matrix list e(wdidsgresults_1), format(%12.7f)
dis "A pwsgs:"
cap matrix list e(pwsgs), format(%12.7f)
dis "===ARM_A_END==="

*================= ARM B: sector partially MISSING, noiseless ===========
dis "===ARM_B_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoiseless id year ttre, att att_i att_it att_itime subgroups(sector_miss)
dis "ARM_B_RC=" _rc
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_i)
scalar B_sg1_att_i = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_itime)
scalar B_sg1_att_itime = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att)
scalar B_sg1_att = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_it)
scalar B_sg1_att_it = r(val)
l0get, mat(e(wdidsgresults_0)) eq(main) nm(att_i)
scalar B_sg0_att_i = r(val)
l0get, mat(e(wdidmainresults)) eq(main) nm(att_i)
scalar B_all_att_i = r(val)
scalar B_N = e(N)
scalar B_Ndrop = e(N_dropped)
cap matrix GTB = e(gteffects)
dis "B pwsgs:"
cap matrix list e(pwsgs), format(%12.7f)
dis "B ereturn scalars:"
cap ereturn list
dis "===ARM_B_END==="

*================= ARM B2: same, explicit safety(off) ===================
dis "===ARM_B2_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoiseless id year ttre, att att_i att_it att_itime subgroups(sector_miss) safety(off)
dis "ARM_B2_RC=" _rc
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_i)
scalar B2_sg1_att_i = r(val)
dis "===ARM_B2_END==="

*================= ARM C: sector CHANGES 1->0, default safety ===========
dis "===ARM_C_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoiseless id year ttre, att att_i att_it att_itime subgroups(sector_chg)
dis "ARM_C_RC=" _rc
scalar C_rc = _rc
dis "===ARM_C_END==="

*================= ARM C2: sector CHANGES 1->0, safety(off) =============
dis "===ARM_C2_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoiseless id year ttre, att att_i att_it att_itime subgroups(sector_chg) safety(off)
dis "ARM_C2_RC=" _rc
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_i)
scalar C2_sg1_att_i = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_itime)
scalar C2_sg1_att_itime = r(val)
dis "===ARM_C2_END==="

*================= ARM D0: clean sector, NOISY ==========================
dis "===ARM_D0_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoisy id year ttre, att att_i att_it att_itime subgroups(sector)
dis "ARM_D0_RC=" _rc
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_i)
scalar D0_sg1_att_i = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_itime)
scalar D0_sg1_att_itime = r(val)
dis "D0 pwsgs:"
cap matrix list e(pwsgs), format(%12.7f)
dis "===ARM_D0_END==="

*================= ARM D1: sector partially MISSING, NOISY ==============
dis "===ARM_D1_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoisy id year ttre, att att_i att_it att_itime subgroups(sector_miss)
dis "ARM_D1_RC=" _rc
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_i)
scalar D1_sg1_att_i = r(val)
l0get, mat(e(wdidsgresults_1)) eq(main) nm(att_itime)
scalar D1_sg1_att_itime = r(val)
dis "D1 pwsgs:"
cap matrix list e(pwsgs), format(%12.7f)
dis "D1 sgresults_1:"
cap matrix list e(wdidsgresults_1), format(%12.7f)
dis "===ARM_D1_END==="

*========= ARM E0/E1: single-firm subgroups = independent noisy oracle ==
* With ONE firm in a subgroup, share weighting and equal-unit weighting
* coincide identically, so these runs measure each firm's own post-cell
* mean without assuming which theory is right.
dis "===ARM_E0_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoisy id year ttre, att att_i att_itime subgroups(sgu_full)
dis "ARM_E0_RC=" _rc
forval g = 1/6 {
  l0get, mat(e(wdidsgresults_`g')) eq(main) nm(att_i)
  scalar E0_`g' = r(val)
  dis "E0_`g' = " %20.10f E0_`g'
}
dis "===ARM_E0_END==="

dis "===ARM_E1_BEGIN==="
use "`base'", clear
cap noisily wooldid ynoisy id year ttre, att att_i att_itime subgroups(sgu_lab)
dis "ARM_E1_RC=" _rc
forval g = 1/6 {
  l0get, mat(e(wdidsgresults_`g')) eq(main) nm(att_i)
  scalar E1_`g' = r(val)
  dis "E1_`g' = " %20.10f E1_`g'
}
dis "===ARM_E1_END==="

*================= VERDICT BLOCK ========================================
dis "===VERDICT_BEGIN==="
dis "PRED_EQUALUNIT   = " %20.10f EQUALUNIT
dis "PRED_SHAREWTD    = " %20.10f SHAREWTD
dis "PRED_ROWWTD      = " %20.10f ROWWTD
dis "PRED_ATTCLEAN    = " %20.10f ATTCLEAN
dis ""
dis "A_sg1_att_i      = " %20.10f A_sg1_att_i
dis "A_sg1_att_itime  = " %20.10f A_sg1_att_itime
dis "A_sg1_att        = " %20.10f A_sg1_att
dis "A_sg1_att_it     = " %20.10f A_sg1_att_it
dis "A_sg0_att_i      = " %20.10f A_sg0_att_i
dis "A_all_att_i      = " %20.10f A_all_att_i
dis "A_N              = " A_N "   A_Ndrop = " A_Ndrop
dis ""
dis "B_sg1_att_i      = " %20.10f B_sg1_att_i
dis "B_sg1_att_itime  = " %20.10f B_sg1_att_itime
dis "B_sg1_att        = " %20.10f B_sg1_att
dis "B_sg1_att_it     = " %20.10f B_sg1_att_it
dis "B_sg0_att_i      = " %20.10f B_sg0_att_i
dis "B_all_att_i      = " %20.10f B_all_att_i
dis "B_N              = " B_N "   B_Ndrop = " B_Ndrop
dis ""
dis "DEV_B_vs_EQUALUNIT = " %20.12e (B_sg1_att_i - EQUALUNIT)
dis "DEV_B_vs_SHAREWTD  = " %20.12e (B_sg1_att_i - SHAREWTD)
dis "DEV_B_vs_ROWWTD    = " %20.12e (B_sg1_att_i - ROWWTD)
dis "DEV_Bitime_vs_SHARE= " %20.12e (B_sg1_att_itime - SHAREWTD)
dis "DEV_A_vs_EQUALUNIT = " %20.12e (A_sg1_att_i - EQUALUNIT)
dis "DEV_Batt_vs_ROWWTD = " %20.12e (B_sg1_att - ROWWTD)
dis ""
dis "B2_sg1_att_i     = " %20.10f B2_sg1_att_i "  (safety(off), should equal B)"
dis "C_rc             = " C_rc "   (expect 459 if the changing case is refused)"
dis "C2_sg1_att_i     = " %20.10f C2_sg1_att_i
dis "C2_sg1_att_itime = " %20.10f C2_sg1_att_itime
dis "HEX B_sg1_att_i  = " %21x B_sg1_att_i
dis "HEX C2_sg1_att_i = " %21x C2_sg1_att_i
dis "B_MINUS_C2       = " %20.12e (B_sg1_att_i - C2_sg1_att_i)
dis ""
dis "gteffects A vs B mreldif = " %20.12e mreldif(GTA,GTB)
dis ""
dis "D0_sg1_att_i     = " %20.10f D0_sg1_att_i
dis "D1_sg1_att_i     = " %20.10f D1_sg1_att_i
dis "D0_sg1_att_itime = " %20.10f D0_sg1_att_itime
dis "D1_sg1_att_itime = " %20.10f D1_sg1_att_itime
scalar ORACLE_EQ = (E0_1+E0_2+E0_3+E0_4+E0_5+E0_6)/6
scalar SW = (1/8)*E1_1 + (2/7)*E1_2 + E1_3 + E1_4 + E1_5 + E1_6
scalar SD = (1/8) + (2/7) + 1 + 1 + 1 + 1
scalar ORACLE_SHARE = SW/SD
dis "NOISY ORACLE equal-unit  = " %20.10f ORACLE_EQ
dis "NOISY ORACLE share-wtd   = " %20.10f ORACLE_SHARE
dis "DEV_D1_vs_ORACLE_EQ      = " %20.12e (D1_sg1_att_i - ORACLE_EQ)
dis "DEV_D1_vs_ORACLE_SHARE   = " %20.12e (D1_sg1_att_i - ORACLE_SHARE)
dis "DEV_D0_vs_ORACLE_EQ      = " %20.12e (D0_sg1_att_i - ORACLE_EQ)
dis "===VERDICT_END==="
