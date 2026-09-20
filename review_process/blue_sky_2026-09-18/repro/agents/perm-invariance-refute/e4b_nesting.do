*==========================================================================
* e4b_nesting.do -- SELF-CONTAINED.  NEW probe (not seeded).
*
* CLAIM UNDER TEST: permute()'s cluster-nesting screen (wooldid_permctl,
* ado:5969-5981) evaluates onset consistency on `insmp' = (if==1 & cluster
* nonmissing), which INCLUDES rows wooldid's own markout will discard.  A
* panel that estimates cleanly -- because the offending unit contributes no
* estimable row -- is therefore REFUSED r(459) the moment permute() is added,
* and the error's observation count is dominated by rows no fit can use.
*
* Shape: 8 states x 3 counties.  State 3 adopts in 2006, but its county 9 is
* absent from BOTH the policy file and the outcome file (a merge miss kept in
* the panel skeleton): ttre == . and y == . on every one of its rows.  Onset
* is constant within every i, so wooldid's general within-i screen (ado:772)
* does not fire.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do e4b_nesting.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

clear
set seed 4242
set sortseed 4242
qui set obs 24
qui gen long unit  = _n
qui gen long state = ceil(unit/3)
qui gen double ai  = rnormal()
qui gen double ttre = .
qui replace ttre = 2004 if state == 1
qui replace ttre = 2006 if state == 2
qui replace ttre = 2006 if state == 3
qui expand 8
qui bysort unit: gen int t = 1999 + _n
qui gen double eps  = rnormal()*0.4
qui gen byte   post = (t >= ttre) & !missing(ttre)
qui gen double y    = ai + 0.15*(t-2000) + 0.9*post + eps
qui drop eps
tempfile base
save `base'

*--- V0: clean reference, no merge miss
use `base', clear
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permexact dots(0)
dis "RES V0_clean rc=" _rc " N=" e(N) " Nclust=" e(N_clust) " npos=" e(perm_npossible)

*--- V1: county 9 of state 3 is a merge miss: ttre AND y missing on every row.
*        (a) estimation WITHOUT permute() -- expect clean
use `base', clear
qui replace ttre = . if unit == 9
qui replace y    = . if unit == 9
cap noi wooldid y unit t ttre , att cluster(state)
dis "RES V1a_noperm rc=" _rc " N=" e(N) " Nclust=" e(N_clust)
dis "RES V1a_noperm b=" %21x e(wdidmainresults)[1,1]

*        (b) SAME data, permute() added
use `base', clear
qui replace ttre = . if unit == 9
qui replace y    = . if unit == 9
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permexact dots(0)
dis "RES V1b_perm rc=" _rc

*        (c) the estimation sample is unchanged if the dead rows are DROPPED,
*            and permute() then runs.  So the refusal is caused only by rows
*            that cannot enter any fit.
use `base', clear
qui drop if unit == 9
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permexact dots(0)
dis "RES V1c_dropped rc=" _rc " N=" e(N) " Nclust=" e(N_clust) " npos=" e(perm_npossible)
dis "RES V1c_dropped b=" %21x e(wdidmainresults)[1,1]

*        (d) and the point estimate in (a) equals the one in (c): same sample
use `base', clear
qui replace ttre = . if unit == 9
qui replace y    = . if unit == 9
cap noi wooldid y unit t ttre , att cluster(state)
scalar bA = e(wdidmainresults)[1,1]
scalar nA = e(N)
use `base', clear
qui drop if unit == 9
cap noi wooldid y unit t ttre , att cluster(state)
dis "RES V1d_same b_keep=" %21x bA " b_drop=" %21x e(wdidmainresults)[1,1] " identical=" (bA == e(wdidmainresults)[1,1]) " N_keep=" nA " N_drop=" e(N)

*--- V2: same shape but the dead county carries a DIFFERENT ONSET rather than
*        a missing one (a stale policy code on an unmergeable unit)
use `base', clear
qui replace ttre = 2010 if unit == 9
qui replace y    = .    if unit == 9
cap noi wooldid y unit t ttre , att cluster(state)
dis "RES V2a_noperm rc=" _rc " N=" e(N)
use `base', clear
qui replace ttre = 2010 if unit == 9
qui replace y    = .    if unit == 9
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permexact dots(0)
dis "RES V2b_perm rc=" _rc

dis "DONE_E4"
