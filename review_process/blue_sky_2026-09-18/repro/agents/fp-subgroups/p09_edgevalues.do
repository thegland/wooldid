* ===========================================================================
* p09 -- SUBGROUP VALUE EDGE CASES.
*   N1  a subgroup level that exists ONLY among never-treated units
*       (realistic: a region with no treated unit).  Does margins' over()
*       still line up with sglist, and are the contrasts still right?
*   N2  a NEGATIVE subgroup value on treated rows -- informative refusal, or
*       a mangled name "main_SG-1_att" and a bare r(198)?
*   N3  a NON-INTEGER subgroup value on treated post rows.
*   N4  8-digit subgroup values WITH att_itime: the ado comment at :12283 says
*       "pretreatment_SG<value>_att_itime rownames all crash at 8-digit values".
*       Caution first, then a clean error?  Or a silently truncated / colliding
*       stored name?
*   N5  7-digit subgroup values with att_itime (the name is exactly 32 chars).
*   N6  a subgroup with NO estimable pre-treatment period.
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

local gN 240
local gT 9
local gt0 2000
local gonsets "2003 2004 2006"
local gysd 0.60
local gtsd 0.50
local gaisd 1
include "`GEN'/gen_core.do"

egen long grp = group(cohort region)
qui save base, replace

* ============ N1: a level only among never-treated units ============
dis "@@@@@@@@@@ N1 @@@@@@@@@@"
use base, clear
gen byte sgn1 = cond(ttre >= ., 9, cond(region <= 2, 0, 1))
tab sgn1 cohort, m
cap noi wooldid y grp t ttre , att att_it cluster(unit) fe(i.unit) ///
    subgroups(sgn1) esprelength(2) espostlength(2)
dis "N1_RC=" _rc
cap dis "N1 sglist = " e(wdid_sglist)
cap mat list e(pwsgs)
cap mat list e(wdidsgresults_0)
cap mat list e(wdidsgresults_1)
cap mat list e(wdidsgresults_9)
matrix N1_S0 = e(wdidsgresults_0)
matrix N1_S1 = e(wdidsgresults_1)
matrix N1_PW = e(pwsgs)
dis "N1 b0-b1 check: b1-b0 = " %20.16f (N1_S1[1,1]-N1_S0[1,1]) "  pw = " %20.16f N1_PW[1,1]

* ============ N2: negative subgroup value ============
dis "@@@@@@@@@@ N2 @@@@@@@@@@"
use base, clear
gen int sgn2 = cond(region <= 2, -1, 1)
cap noi wooldid y grp t ttre , att att_it cluster(unit) fe(i.unit) ///
    subgroups(sgn2) esprelength(2) espostlength(2)
dis "N2_RC=" _rc

* ============ N3: non-integer subgroup value on post rows ============
dis "@@@@@@@@@@ N3 @@@@@@@@@@"
use base, clear
gen double sgn3 = cond(region <= 2, 1.5, 2)
cap noi wooldid y grp t ttre , att att_it cluster(unit) fe(i.unit) ///
    subgroups(sgn3) esprelength(2) espostlength(2)
dis "N3_RC=" _rc

* ============ N4: 8-digit values + att_itime ============
dis "@@@@@@@@@@ N4 @@@@@@@@@@"
use base, clear
gen double sgn4 = cond(region <= 2, 48000001, 48000002)
cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) fe(i.unit) ///
    subgroups(sgn4) esprelength(2) espostlength(2)
dis "N4_RC=" _rc
cap dis "N4 sglist = " e(wdid_sglist)
cap matrix N4B = e(b)
cap local n4n : colnames N4B
cap dis "N4 colnames: `n4n'"

* ============ N4b: 8-digit values, att only (short names) ============
dis "@@@@@@@@@@ N4b @@@@@@@@@@"
use base, clear
gen double sgn4 = cond(region <= 2, 48000001, 48000002)
cap noi wooldid y grp t ttre , att cluster(unit) fe(i.unit) ///
    subgroups(sgn4) esprelength(2) espostlength(2)
dis "N4b_RC=" _rc
cap matrix N4bB = e(b)
cap local n4bn : colnames N4bB
cap dis "N4b colnames: `n4bn'"
cap mat list e(pwsgs)

* ============ N5: 7-digit values + att_itime ============
dis "@@@@@@@@@@ N5 @@@@@@@@@@"
use base, clear
gen double sgn5 = cond(region <= 2, 4800001, 4800002)
cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) fe(i.unit) ///
    subgroups(sgn5) esprelength(2) espostlength(2)
dis "N5_RC=" _rc
cap matrix N5B = e(b)
cap local n5n : colnames N5B
cap dis "N5 colnames: `n5n'"
cap dis "N5 name length check: " strlen("pretreatment_SG4800001_att_itime")

* ============ N6: subgroup with no estimable pre period ============
dis "@@@@@@@@@@ N6 @@@@@@@@@@"
use base, clear
* subgroup 1 == the cohort whose onset is the FIRST usable period, so it has
* essentially no pre-treatment window
gen byte sgn6 = (ttre == 2003)
qui replace sgn6 = 0 if ttre >= .
qui drop if ttre == 2003 & t < 2002
tab sgn6 ttre, m
cap noi wooldid y grp t ttre , att att_it cluster(unit) fe(i.unit) ///
    subgroups(sgn6) esprelength(3) espostlength(2)
dis "N6_RC=" _rc
cap mat list e(pwsgs)
cap dis "N6 numpre_att_s0=" e(numpre_att_s0) " numpre_att_s1=" e(numpre_att_s1)

* ============ N7: subgroup value missing on SOME of a unit's rows ============
* Realistic: item nonresponse in a subgroup variable (the unit's sg is known in
* some years and missing in others).  wooldid REFUSES a unit whose sg CHANGES
* over time, and its own error text says why: att_i / att_itime would then enter
* the unit with a SHARE weight rather than the equal unit weight those estimands
* are defined by.  The partially-missing case produces exactly that same share
* weighting -- but passes every screen silently.  Measure which weighting is used.
dis "@@@@@@@@@@ N7 @@@@@@@@@@"
use base, clear
gen byte sgn7 = (region <= 2)
* for region-1 units, blank the subgroup in the second half of the panel
qui replace sgn7 = . if region == 1 & t >= 2005
tab sgn7, m
cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) fe(i.unit) ///
    subgroups(sgn7) esprelength(2) espostlength(2)
dis "N7_RC=" _rc
qui gen byte insmp7 = e(sample)
matrix N7_S1 = e(wdidsgresults_1)
matrix N7_MAIN = e(wdidmainresults)
dis "N7 att_i for sg1 (wooldid) = " %20.16f N7_S1[3,1]

* rebuild att_i for sg==1 two ways
cap drop __*
qui gen byte __post = insmp7 & ttre < . & t >= ttre
qui gen double __w = 1
qui egen double __wip = total(cond(__post,__w,0)), by(grp)          // ALL post rows of i
qui egen double __wipS = total(cond(__post & sgn7==1,__w,0)), by(grp) // post rows of i INSIDE sg1
* (a) SHARE weighting: w / (total post weight of i over ALL rows)
qui gen double __Wa = __w/__wip  if __post & sgn7==1
* (b) EQUAL-UNIT weighting: w / (total post weight of i INSIDE the subgroup)
qui gen double __Wb = __w/__wipS if __post & sgn7==1
* the per-unit effect wooldid is averaging is not observable, so instead compare
* the WEIGHTS themselves: equal-unit weighting gives each unit total weight 1.
qui egen byte __tagi = tag(grp) if __post & sgn7==1
qui summarize __Wa if __post & sgn7==1, meanonly
local SA = r(sum)
qui summarize __Wb if __post & sgn7==1, meanonly
local SB = r(sum)
qui count if __tagi==1
local NU = r(N)
dis "N7 units in sg1 = `NU';  total SHARE weight = `SA';  total EQUAL-UNIT weight = `SB'"
dis "N7 (if SHARE < units, some unit enters sg1 with less than full weight)"
* per-unit share weights
qui egen double __uw = total(__Wa), by(grp)
qui tabstat __uw if __tagi==1, stat(min mean max) save
matrix N7T = r(StatTotal)
mat list N7T
