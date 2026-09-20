*=============================================================================
* p15  SELF-CONTAINED REPRO.
*
* CLAIM: under bootstrap(bayes) + safety(off), two estimands that are THE SAME
* NUMBER on the data - att and att_it are bit-identical here, with bit-identical
* analytic standard errors - are reported with bootstrap standard errors that
* differ by a factor of ~5. bootstrap(standard) reports them bit-identically, as
* it must. So the bayes standard error for the cell-weighted estimands is not the
* standard error of the statistic that is printed beside it.
*
* Run:  stata-mp -e do p15_attvsattit.do    (needs the pinned build on the adopath)
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
which wooldid

*--- a FULLY BALANCED staggered panel: every treated cell holds exactly one row,
*    so the observation-weighted att and the equal-cell-weighted att_it coincide
set seed 20260917
set sortseed 20260917
clear
set obs 24
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2005 if st <= 6
replace onset = 2008 if st > 6 & st <= 12
replace onset = 2011 if st > 12 & st <= 16
* st 17-24 never treated
* treatment effect varies ACROSS clusters - the regime the package's own caution names
gen double shk = rnormal()
expand 12
bysort st: gen int year = 2000 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.03*(year-2000) + (1.0 + 1.2*shk)*post + rnormal()*0.5

*--- 1. analytic: are att and att_it the same number? -----------------------
qui wooldid y st year onset, att att_it cluster(st)
matrix A = e(wdidmainresults)
dis "ANALYTIC  att    b = " %21x A[1,1] "   se = " %21x A[1,2]
dis "ANALYTIC  att_it b = " %21x A[2,1] "   se = " %21x A[2,2]
dis "ANALYTIC  identical point estimates: " (A[1,1]==A[2,1]) "   identical ses: " (A[1,2]==A[2,2])

*--- 2. bootstrap(standard): must report them identically -------------------
qui wooldid y st year onset, att att_it cluster(st) bootstrap(standard) reps(400) rseed(31337) dots(0)
matrix VS = e(V)
matrix BS = e(b)
dis "STANDARD  att    b = " %21x BS[1,1] "  se = " %12.8f sqrt(VS[1,1])
dis "STANDARD  att_it b = " %21x BS[1,2] "  se = " %12.8f sqrt(VS[2,2])
dis "STANDARD  se ratio att_it/att = " %8.4f sqrt(VS[2,2]/VS[1,1])

*--- 3. bootstrap(bayes) + safety(off): the same two numbers ----------------
qui wooldid y st year onset, att att_it cluster(st) bootstrap(bayes) reps(400) rseed(31337) dots(0) safety(off)
matrix VB = e(V)
matrix BB = e(b)
dis "BAYES     att    b = " %21x BB[1,1] "  se = " %12.8f sqrt(VB[1,1])
dis "BAYES     att_it b = " %21x BB[1,2] "  se = " %12.8f sqrt(VB[2,2])
dis "BAYES     se ratio att_it/att = " %8.4f sqrt(VB[2,2]/VB[1,1])
dis "BAYES vs STANDARD: se(att)    ratio = " %8.4f sqrt(VB[1,1]/VS[1,1])
dis "BAYES vs STANDARD: se(att_it) ratio = " %8.4f sqrt(VB[2,2]/VS[2,2])
dis "ANALYTIC se        = " %12.8f A[1,2]
dis "PROBE_DONE"
