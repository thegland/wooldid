*==========================================================================
* repro_SORTORDER.do
*   RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_SORTORDER.do
*   Self-contained.  Pinned build only.  No option is typed anywhere below: both
*   runs are a PLAIN DEFAULT wooldid call.
*
* CLAIM: under the default cell collapse, the standard errors, t statistics,
* p-values and confidence bounds wooldid PRINTS depend on the order in which the
* user's observations happen to be stored.  Sorting the same dataset differently
* and re-running moves the pretreatment standard error by 3.3e-06 relative, into
* the printed digits.  With nofastcollapse the same permutation moves nothing
* (1.3e-15).  The direction is the reverse of the calibration paragraph at
* wooldid.sthlp:357, which tells the reader that a row permutation of the
* uncollapsed path moves a result about 19 times MORE than the collapse does.
*
* Cause (verified separately with a one-line patched copy of the ado): the
* collapsed row keeps the cell's FIRST row's raw controls() value.  Re-sorting
* changes which row that is.  Carrying the cell's weighted MEAN instead - what
* wooldid_fcol already does for the contreat() dose - removes the dependence.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 515151
set sortseed 515151

set obs 150
gen long unit = _n
gen double _u = runiform()
sort _u unit
gen long rk = _n
gen long cohort = 4
gen double ttre = .
replace cohort = 1 in 1/33
replace ttre = 2004 in 1/33
replace cohort = 2 in 34/66
replace ttre = 2006 in 34/66
replace cohort = 3 in 67/98
replace ttre = 2008 in 67/98
gen double a_i = rnormal()
drop _u rk
expand 10
bysort unit: gen int t = 1999 + _n
sort t unit
by t: gen double _f = rnormal() if _n==1
by t: replace _f = _f[1]
gen double ft = 0.5*_f
drop _f
gen double x1 = rnormal()
gen double x2 = 0.5*a_i + 0.5*rnormal()
gen byte post = (ttre < . & t >= ttre)
gen double tau = 0
replace tau = 0.50 if cohort==1 & post
replace tau = 0.80 if cohort==2 & post
replace tau = 0.30 if cohort==3 & post
gen double y = a_i + ft + tau + 0.5*x1 - 0.25*x2 + 0.30*rnormal()
* an ordinary percentage-scale covariate that varies across units within a
* cohort-year cell, e.g. a county's uninsured rate in a state-level design
gen double xpct = 50 + 8*x1
compress
save "_so.dta", replace

dis _n as text "RUN 1: the data as built (sorted by unit, t)"
use "_so.dta", clear
sort unit t
wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xpct)
matrix M1 = e(wdidmainresults)
matrix V1 = e(V)

dis _n as text "RUN 2: the SAME data, sorted differently.  Same command, no options changed."
use "_so.dta", clear
gsort -t -unit
wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xpct)
matrix M2 = e(wdidmainresults)
matrix V2 = e(V)

dis _n as text "RUN 3 / RUN 4: the same two sorts under nofastcollapse"
use "_so.dta", clear
sort unit t
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xpct) nofastcollapse
matrix M3 = e(wdidmainresults)
use "_so.dta", clear
gsort -t -unit
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xpct) nofastcollapse
matrix M4 = e(wdidmainresults)

mata:
function mx(string scalar an, string scalar bn)
{
    real matrix A, B, S, Dm
    A = st_matrix(an)
    B = st_matrix(bn)
    S  = abs(A)
    S  = S :* (S:>=abs(B)) + abs(B) :* (S:<abs(B))
    Dm = abs(A-B) :/ (S :+ (S:==0))
    Dm = editmissing(Dm, 0)
    return(max(Dm))
}
printf("\nmax relative move, DEFAULT collapse, sort 1 vs sort 2 :  e(wdidmainresults) %9.2e   e(V) %9.2e\n", mx("M1","M2"), mx("V1","V2"))
printf("max relative move, NOFASTCOLLAPSE,  sort 1 vs sort 2 :  e(wdidmainresults) %9.2e\n", mx("M3","M4"))
end

dis as text "Compare the two printed pretreatment blocks above: the standard errors and"
dis as text "confidence bounds differ between two runs of the same command on the same data."

erase "_so.dta"
dis "REPRO_SORTORDER_DONE"
