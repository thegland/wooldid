* ==========================================================================
* r06_attit.do -- NEW probe (mine).  Both sthlp:887 and wooldid's own r(459)
* error text assert that "by-subgroup att and att_it are unaffected" by the
* share-weighting problem.  That is TRUE when every (i,t) cell holds exactly
* one row.  It is FALSE as soon as a cell holds several rows and only some of
* them carry the subgroup label -- the normal shape when i is a COHORT and the
* rows are individuals, which is the geometry wooldid's own help recommends
* (a coarser i with unit fixed effects through fe()).
*
* att_it is documented as "the simple average of cohort-period effects" -- each
* treatment CELL equally weighted.  A cell labelled on 1 of its 5 rows should
* still be one cell.
*
* DESIGN (noiseless, so every target is exact):
*   i=1  treated t>=5, sg=1, effect 10, 4 post cells of 5 rows, ALL labelled
*   i=2  treated t>=5, sg=1, effect  2, post cells of 5,5,2,2 rows,
*                                      labelled 1,1,2,2 rows respectively
*   i=3  treated t>=5, sg=0, effect  1
*   i=4,5,6 never treated
*
*   cell shares for i=2: .2 .2 1 1  (sum 2.4)      for i=1: 1 1 1 1 (sum 4)
*
*   TARGETS for subgroup 1
*     att_it  equal-cell (documented) = 48/8            = 6.0000000000
*     att_it  cell-share-weighted     = 44.8/6.4        = 7.0000000000
*     att_i   equal-unit (documented) = 12/2            = 6.0000000000
*     att_i   share-weighted          = 10.857143/1.428571 = 7.6000000000
*     att_itime share-weighted        = 11.2/1.6        = 7.0000000000
*     att     row-weighted            = 212/26          = 8.1538461538
* ==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- (i,t) skeleton ----
clear
set obs 6
gen long id = _n
expand 8
bysort id: gen int t = _n
* rows per (i,t)
gen int nrow = 5
replace nrow = 2 if id == 2 & t >= 7
expand nrow
bysort id t: gen int k = _n

gen double ttre = .
replace ttre = 5 if id <= 3
gen double eff = 0
replace eff = 10 if id == 1
replace eff = 2  if id == 2
replace eff = 1  if id == 3
gen double y = 7*id + 0.5*t + eff*(t >= ttre & ttre < .)

* clean subgroup: constant within (i,t) and within i
gen byte sgA = .
replace sgA = 1 if id <= 2
replace sgA = 0 if id == 3
replace sgA = 0 if id >= 4

* partially-labelled cells: i=2's post cells keep the label on 1,1,2,2 rows
gen byte sgB = sgA
replace sgB = . if id == 2 & t == 5 & k > 1
replace sgB = . if id == 2 & t == 6 & k > 1
* t==7 and t==8 cells have only 2 rows and both stay labelled

* audit the shares the design intends
preserve
  keep if id == 2 & t >= 5
  gen byte lab = (sgB == 1)
  collapse (sum) lab (count) tot = k, by(t)
  gen double share = lab/tot
  list t lab tot share, noobs
restore

foreach arm in sgA sgB {
  dis "############### ARM `arm' ###############"
  dis "CAUTIONSCAN_`arm'_START"
  capture noisily wooldid y id t ttre, att att_it att_i att_itime subgroups(`arm')
  dis "ARM_`arm'_RC=" _rc
  dis "CAUTIONSCAN_`arm'_END"
  matrix M_`arm' = e(wdidsgresults_1)
  local rn : rowfullnames M_`arm'
  dis "ROWNAMES_`arm': `rn'"
  matrix list M_`arm', format(%14.10f)
}

scalar T_attit_equal = 48/8
scalar T_attit_share = 44.8/6.4
scalar T_atti_equal  = 6
scalar T_atti_share  = (10 + 2*(6/14))/(1 + 6/14)
scalar T_attim_share = 11.2/1.6
scalar T_att_row     = 212/26

dis "================= RESULTS ================="
dis "CLEAN   att       = " %16.10f M_sgA[1,1]
dis "CLEAN   att_it    = " %16.10f M_sgA[2,1]
dis "CLEAN   att_i     = " %16.10f M_sgA[3,1]
dis "CLEAN   att_itime = " %16.10f M_sgA[4,1]
dis "PARTIAL att       = " %16.10f M_sgB[1,1]
dis "PARTIAL att_it    = " %16.10f M_sgB[2,1]
dis "PARTIAL att_i     = " %16.10f M_sgB[3,1]
dis "PARTIAL att_itime = " %16.10f M_sgB[4,1]
dis "-------------------------------------------"
dis "TARGET att_it equal-cell   = " %16.10f T_attit_equal
dis "TARGET att_it cell-share   = " %16.10f T_attit_share
dis "TARGET att_i  equal-unit   = " %16.10f T_atti_equal
dis "TARGET att_i  share        = " %16.10f T_atti_share
dis "TARGET att_itime share     = " %16.10f T_attim_share
dis "TARGET att    row-weighted = " %16.10f T_att_row
dis "-------------------------------------------"
dis "DEV attit_vs_EQUALCELL = " %16.8e (M_sgB[2,1]-T_attit_equal)
dis "DEV attit_vs_SHARECELL = " %16.8e (M_sgB[2,1]-T_attit_share)
dis "DEV atti_vs_EQUALUNIT  = " %16.8e (M_sgB[3,1]-T_atti_equal)
dis "DEV atti_vs_SHARE      = " %16.8e (M_sgB[3,1]-T_atti_share)
dis "DEV attim_vs_SHARE     = " %16.8e (M_sgB[4,1]-T_attim_share)
dis "DEV att_vs_ROW         = " %16.8e (M_sgB[1,1]-T_att_row)
dis "CLEANDEV attit_EQUAL   = " %16.8e (M_sgA[2,1]-6)
dis "CLEANDEV atti_EQUAL    = " %16.8e (M_sgA[3,1]-6)
dis "ALLDONE"
