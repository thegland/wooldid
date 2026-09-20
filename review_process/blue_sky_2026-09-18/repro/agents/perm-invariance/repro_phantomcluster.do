*==========================================================================
* repro_phantomcluster.do  -- SELF-CONTAINED
*
* CLAIM: under permute(), a cluster that contributes NOTHING to the
* estimation sample (every one of its rows has a missing outcome) is still
* included in the randomization design.  Two datasets with a BIT-IDENTICAL
* estimation sample and identical point estimates therefore report
* different permutation p-values and different confidence intervals, with
* no caution.  Excluding the same cluster through a missing CLUSTER value,
* or through the user's `if', does remove it from the design -- so the
* three exclusion routes are inconsistent with one another.
*
* Run:
*   cd <this directory>
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_phantomcluster.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*--- 40 counties in 10 states, 10 annual periods, state-level adoption
clear
set seed 90210
set sortseed 90210
qui set obs 40
qui gen long unit  = _n
qui gen long state = ceil(unit/4)
qui gen double ai  = rnormal()
qui gen double ttre = .
qui replace ttre = 2004 if state <= 3
qui replace ttre = 2006 if state > 3 & state <= 6
qui replace ttre = 2008 if state == 7
qui expand 10
qui bysort unit: gen int t = 1999 + _n
qui gen double eps  = rnormal()*0.5
qui gen byte   post = (t >= ttre) & !missing(ttre)
qui gen double y    = ai + 0.2*(t-2000) + 0.8*post + eps
qui drop eps
tempfile base
save `base'

local OPT att att_i cluster(state) permute(shuffle) permreps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0)
scalar A_b = .
scalar B_b = .

*--- (A) state 5 (a 2006 adopter) has a MISSING OUTCOME on every row; rows kept
use `base', clear
qui replace y = . if state == 5
qui wooldid y unit t ttre , `OPT'
scalar A_b = e(wdidmainresults)[1,1]
local A_p   = e(perm_p_att)
local A_lo  = e(wdidmainresults)[1,5]
local A_hi  = e(wdidmainresults)[1,6]
local A_n   = e(N)
local A_g   = e(N_clust)
local A_npos= e(perm_npossible)
matrix PA = e(perm_p)

*--- (B) the same rows physically dropped
use `base', clear
qui drop if state == 5
qui wooldid y unit t ttre , `OPT'
scalar B_b = e(wdidmainresults)[1,1]
local B_p   = e(perm_p_att)
local B_lo  = e(wdidmainresults)[1,5]
local B_hi  = e(wdidmainresults)[1,6]
local B_n   = e(N)
local B_g   = e(N_clust)
local B_npos= e(perm_npossible)
matrix PB = e(perm_p)

*--- (C) the same cluster removed by setting its CLUSTER value to missing
use `base', clear
qui replace state = . if state == 5
qui wooldid y unit t ttre , `OPT'
local C_p   = e(perm_p_att)
local C_npos= e(perm_npossible)

*--- (D) the same cluster removed by the user's if
use `base', clear
qui wooldid y unit t ttre if state != 5 , `OPT'
local D_p   = e(perm_p_att)
local D_npos= e(perm_npossible)

dis _n "=================== RESULT ==================="
dis "estimation sample   A: N=`A_n' G=`A_g'   B: N=`B_n' G=`B_g'"
dis "att point estimate  A: " %21x A_b "   B: " %21x B_b
dis "                    identical = " (A_b == B_b)
dis "perm_npossible      A: `A_npos'   B: `B_npos'   C: `C_npos'   D: `D_npos'"
dis "permutation p(att)  A: " %9.6f `A_p' "  B: " %9.6f `B_p' "  C: " %9.6f `C_p' "  D: " %9.6f `D_p'
dis "95% interval        A: [" %7.4f `A_lo' ", " %7.4f `A_hi' "]"
dis "                    B: [" %7.4f `B_lo' ", " %7.4f `B_hi' "]"
dis _n "per-coefficient permutation p-values, A (rows kept) then B (rows dropped):"
matlist PA, format(%9.6f)
matlist PB, format(%9.6f)
dis "=============================================="
dis "EXPECTED: A and B agree -- they have the same estimation sample and the"
dis "same estimates, and differ only in rows the estimator discards."
dis "OBSERVED: A keeps 10 clusters in the randomization design, B keeps 9;"
dis "the p-values and intervals differ.  C and D agree with B, not with A."

*--- the same asymmetry in the bootstrap frame, which shares the design population
use `base', clear
qui replace y = . if state == 5
qui wooldid y unit t ttre , att att_i cluster(state) bootstrap(standard) reps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0)
local AB_se = e(V)[1,1]
local AB_g  = e(boot_nclus)
use `base', clear
qui drop if state == 5
qui wooldid y unit t ttre , att att_i cluster(state) bootstrap(standard) reps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0)
local BB_se = e(V)[1,1]
local BB_g  = e(boot_nclus)
dis "bootstrap V[1,1]    A: " %12.8f `AB_se' "   B: " %12.8f `BB_se'
dis "PI1_DONE" 
