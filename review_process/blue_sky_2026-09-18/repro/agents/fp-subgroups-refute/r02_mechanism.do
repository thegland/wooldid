* ==========================================================================
* r02_mechanism.do -- separates the THREE candidate weightings, and tests the
* PRE-period analogue nobody has probed.
*
* In r01 every unit had the same number of post rows, so "share of own post
* weight" and "row count" coincide.  Here they are pulled apart with an
* UNBALANCED panel (a real administrative panel shape: units enter and leave).
*
*   sg==1 units      eff   post rows IN PANEL   post rows WITH sg==1   share
*     unit 1          1            8                    8             1.00
*     unit 2          5            2 (panel ends t=10)  2             1.00
*     unit 3          9            8                    2 (sg misses) 0.25
*
*   equal-unit  (att_i as documented) = (1+5+9)/3          = 5.0000000000
*   SHARE-weighted                    = 8.25/2.25          = 3.6666666667
*   ROW-weighted (what att/att_it do) = 36/12              = 3.0000000000
*
* Three distinct numbers -> the mechanism is identified, not just the error.
*
* PART 2 repeats the exercise with the missingness confined to PRE-treatment
* rows, to see whether pretreatment:att_i is distorted the same way.
* ==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ------------------------------ PART 1 -------------------------------------
clear
set obs 10
gen long id = _n
gen double ttre = .
replace ttre = 9 if id <= 6
gen byte sg = .
replace sg = 1 if id <= 3
replace sg = 0 if id > 3
gen double eff = 0
replace eff = 1 if id == 1
replace eff = 5 if id == 2
replace eff = 9 if id == 3
replace eff = 2 if id >= 4 & id <= 6
expand 16
bysort id: gen int t = _n
* unit 2's panel ENDS at t=10 (leaves the data) -> only 2 post rows exist
drop if id == 2 & t > 10
gen double y = 3*id + 0.5*t + eff*(t >= ttre & ttre < .)

gen byte sgA = sg
gen byte sgB = sg
replace sgB = . if id == 3 & t >= 11     /* unit 3 keeps 2 of its 8 post rows */

dis "=== PART1 RUN A (clean) ==="
capture noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgA)
dis "P1A_RC=" _rc
matrix A = e(wdidsgresults_1)
matrix list A, format(%12.9f)

dis "=== PART1 RUN B (sg missing on 6 of unit 3's 8 post rows) ==="
dis "CAUTIONSCAN_P1B_START"
capture noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgB)
dis "P1B_RC=" _rc
dis "CAUTIONSCAN_P1B_END"
matrix B = e(wdidsgresults_1)
matrix list B, format(%12.9f)

scalar EQUAL = 5
scalar SHARE = 8.25/2.25
scalar ROWWT = 36/12
dis "TARGET_EQUAL = " %18.10f EQUAL
dis "TARGET_SHARE = " %18.10f SHARE
dis "TARGET_ROWWT = " %18.10f ROWWT
local rn : rowfullnames B
dis "ROWNAMES: `rn'"
dis "P1A_att    = " %18.10f A[1,1]
dis "P1A_att_it = " %18.10f A[2,1]
dis "P1A_att_i  = " %18.10f A[3,1]
dis "P1A_atti_m = " %18.10f A[4,1]
dis "P1B_att    = " %18.10f B[1,1]
dis "P1B_att_it = " %18.10f B[2,1]
dis "P1B_att_i  = " %18.10f B[3,1]
dis "P1B_atti_m = " %18.10f B[4,1]
dis "D_att_i_vs_EQUAL = " %18.10e (B[3,1]-EQUAL)
dis "D_att_i_vs_SHARE = " %18.10e (B[3,1]-SHARE)
dis "D_att_i_vs_ROWWT = " %18.10e (B[3,1]-ROWWT)
dis "D_atim_vs_SHARE  = " %18.10e (B[4,1]-SHARE)
dis "D_A_att_i_EQUAL  = " %18.10e (A[3,1]-EQUAL)

* ------------------------------ PART 2 -------------------------------------
* Missingness confined to PRE-treatment rows.  Does pretreatment:att_i move?
* Give the units differential PRE-trend effects so a distortion is visible.
clear
set obs 10
gen long id = _n
gen double ttre = .
replace ttre = 9 if id <= 6
gen byte sg = .
replace sg = 1 if id <= 3
replace sg = 0 if id > 3
gen double eff  = cond(id==1,1,cond(id==2,5,cond(id==3,9,2)))
gen double pre  = cond(id==1,0.2,cond(id==2,1.0,cond(id==3,3.0,0)))
expand 16
bysort id: gen int t = _n
* a genuine pre-trend deviation at relative year -2 and -3 only (t = 6,7)
gen double y = 3*id + 0.5*t + eff*(t >= ttre & ttre < .) + pre*inrange(t,6,7)*(ttre<.)
gen byte sgA = sg
gen byte sgB = sg
replace sgB = . if id == 3 & inrange(t,1,7)   /* sg unknown before t=8 */

dis "=== PART2 RUN A (clean, pre rows) ==="
capture noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgA) esprelength(4)
dis "P2A_RC=" _rc
matrix PA = e(wdidsgresults_1)
matrix list PA, format(%12.9f)
local rnPA : rowfullnames PA
dis "P2A_ROWNAMES: `rnPA'"

dis "=== PART2 RUN B (sg missing on PRE rows only) ==="
dis "CAUTIONSCAN_P2B_START"
capture noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(sgB) esprelength(4)
dis "P2B_RC=" _rc
dis "CAUTIONSCAN_P2B_END"
matrix PB = e(wdidsgresults_1)
matrix list PB, format(%12.9f)
local rnPB : rowfullnames PB
dis "P2B_ROWNAMES: `rnPB'"

dis "ALLDONE"
