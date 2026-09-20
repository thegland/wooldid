* ==========================================================================
* r01_sgmissing.do  --  SKEPTIC's independent test of claim fpsg-01
*
* CLAIM UNDER TEST: when subgroups() is missing on SOME of a unit's post rows,
* the subgroup att_i / att_itime become SHARE-weighted, not equal-unit averages.
*
* Independent discriminating design (NOT the original prober's):
* FIVE units in subgroup 1 with FIVE DISTINCT effects and FIVE DISTINCT
* observed shares, so equal-unit and share-weighted differ by a lot and the
* share-weighted prediction is one exact number that cannot be hit by accident.
*
*   unit  effect  post rows kept in sg1   share s
*     1      1          1 of 8            0.125
*     2      2          2 of 8            0.250
*     3      3          4 of 8            0.500
*     4      4          8 of 8            1.000
*     5      5          8 of 8            1.000
*
*   equal-unit average  = (1+2+3+4+5)/5   = 3.0000000000
*   share-weighted mean = 11.125/2.875    = 3.8695652174
* ==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid
dis "BUILDCHECK done"

* ---------------- panel ----------------------------------------------------
clear
set obs 15
gen long id = _n
gen double ttre = .
replace ttre = 9 if id <= 10
gen byte sg = .
replace sg = 1 if id <= 5
replace sg = 0 if id > 5
gen double eff = 0
replace eff = id       if id <= 5
replace eff = 9        if id > 5 & id <= 10
expand 16
bysort id: gen int t = _n
gen double y = 3*id + 0.5*t + eff*(t >= ttre & ttre < .)

gen byte sgA = sg
gen byte sgB = sg
replace sgB = . if id == 1 & t >= 10
replace sgB = . if id == 2 & t >= 11
replace sgB = . if id == 3 & t >= 13
gen byte sgC = sg
replace sgC = 0 if id == 1 & t >= 10
replace sgC = 0 if id == 2 & t >= 11
replace sgC = 0 if id == 3 & t >= 13

preserve
  keep if t >= 9 & ttre < . & sgB == 1
  bysort id: gen n = _N
  bysort id: keep if _n == 1
  list id n, noobs
restore

* ---------------- RUN C FIRST (the refusal), so nothing later can lose it ---
dis "=============== RUN C (sg changes 1->0, same rows) ==============="
cap noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgC)
dis "RUNC_RC=" _rc

* ---------------- RUN A : clean subgroup -----------------------------------
dis "=============== RUN A (clean sg) ==============="
dis "CAUTIONSCAN_A_START"
cap noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgA) histogramcohorteffects
dis "RUNA_RC=" _rc
dis "CAUTIONSCAN_A_END"
capture noisily mat A = e(wdidsgresults_1)
dis "GOT_A_RC=" _rc
capture noisily mat list A, format(%12.9f)
capture noisily mat HA = e(histogramestimates)
capture noisily mat list HA, format(%12.9f)
local rnA : rowfullnames A
dis "ROWNAMES_A: `rnA'"

* ---------------- RUN B : partially MISSING subgroup -----------------------
dis "=============== RUN B (sg missing on some post rows) ==============="
dis "CAUTIONSCAN_B_START"
cap noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgB) histogramcohorteffects
dis "RUNB_RC=" _rc
dis "CAUTIONSCAN_B_END"
capture noisily mat B = e(wdidsgresults_1)
dis "GOT_B_RC=" _rc
capture noisily mat list B, format(%12.9f)
capture noisily mat HB = e(histogramestimates)
capture noisily mat list HB, format(%12.9f)
local rnB : rowfullnames B
dis "ROWNAMES_B: `rnB'"
capture noisily mat GT = e(gteffects)
dis "GOT_GT_RC=" _rc
capture noisily mat list GT, format(%12.6f)

* ---------------- the arithmetic -------------------------------------------
capture noisily scalar bA_i     = A["main:att_i", 1]
capture noisily scalar bA_itime = A["main:att_itime", 1]
capture noisily scalar bB_i     = B["main:att_i", 1]
capture noisily scalar bB_itime = B["main:att_itime", 1]
capture noisily scalar bB_att   = B["main:att", 1]
capture noisily scalar bB_it    = B["main:att_it", 1]
dis "EXTRACT_RC=" _rc

scalar EQUAL  = 3
scalar SHARE  = (1*0.125 + 2*0.25 + 3*0.5 + 4*1 + 5*1)/(0.125+0.25+0.5+1+1)

dis "----------------------------------------------------------"
capture noisily dis "RUNA_att_i     = " %18.10f bA_i
capture noisily dis "RUNA_att_itime = " %18.10f bA_itime
capture noisily dis "RUNB_att_i     = " %18.10f bB_i
capture noisily dis "RUNB_att_itime = " %18.10f bB_itime
capture noisily dis "RUNB_att       = " %18.10f bB_att
capture noisily dis "RUNB_att_it    = " %18.10f bB_it
dis "TARGET_EQUAL   = " %18.10f EQUAL
dis "TARGET_SHARE   = " %18.10f SHARE
capture noisily dis "DEV_Bi_EQUAL   = " %18.10e (bB_i - EQUAL)
capture noisily dis "DEV_Bi_SHARE   = " %18.10e (bB_i - SHARE)
capture noisily dis "DEV_Bitime_SHR = " %18.10e (bB_itime - SHARE)
capture noisily dis "DEV_Ai_EQUAL   = " %18.10e (bA_i - EQUAL)
dis "----------------------------------------------------------"


dis "--- positional fallback (rows: 1 att, 2 att_it, 3 att_i, 4 att_itime) ---"
dis "POS_A_att   = " %18.10f A[1,1]
dis "POS_A_att_it= " %18.10f A[2,1]
dis "POS_A_att_i = " %18.10f A[3,1]
dis "POS_A_atti_m= " %18.10f A[4,1]
dis "POS_B_att   = " %18.10f B[1,1]
dis "POS_B_att_it= " %18.10f B[2,1]
dis "POS_B_att_i = " %18.10f B[3,1]
dis "POS_B_atti_m= " %18.10f B[4,1]
dis "POSDEV_Bi_EQUAL = " %18.10e (B[3,1] - EQUAL)
dis "POSDEV_Bi_SHARE = " %18.10e (B[3,1] - SHARE)
dis "POSDEV_Bm_SHARE = " %18.10e (B[4,1] - SHARE)
dis "POSDEV_Ai_EQUAL = " %18.10e (A[3,1] - EQUAL)
dis "POSDEV_Am_EQUAL = " %18.10e (A[4,1] - EQUAL)
dis "ALLDONE"
