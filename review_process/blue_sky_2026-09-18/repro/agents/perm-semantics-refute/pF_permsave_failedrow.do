*! SKEPTIC probe F -- SELF-CONTAINED repro of a documented-vs-actual divergence I found
*! while refuting PERMSEM-01.
*!
*! The help, describing permsave() (wooldid.sthlp line 493), says:
*!   "(_rep = 1 to permreps(), or 1 to M under permexact; a failed replication keeps its
*!    row, ALL MISSING) plus a _rep = 0 row carrying the observed estimates"
*! It is not all missing.  The same paragraph also adds "one onset_<clustervalue> column
*! per cluster holding the onset that cluster carried in each replication", and THOSE
*! columns are fully populated on a failed replication's row.  So a user who identifies
*! the failed replications the way the sentence invites - "the row where everything is
*! missing" - finds NONE of them, on a file that has several.
*! (Contrast bootsave(), sthlp:557, which says "a replication that failed keeps its row
*!  with every value missing" and whose w_<cluster> resample columns raise the same
*!  question; only the permutation file is tested here.)
*!
*! DGP: 10 districts, 2 adopt in 2017, 3 stop reporting after 2014.  A shuffled assignment
*! that puts BOTH onsets on early leavers leaves no treated observation, so that
*! replication fails: P = C(3,2)/C(10,2) = 3/45.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pF_permsave_failedrow.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

clear
set obs 10
gen long dist = _n
gen double onset = cond(inlist(dist,4,9), 2017, .)
gen byte leaver = inlist(dist,2,6,10)
gen double a_i = rnormal()*2
expand 10
bysort dist: gen int year = 2010 + _n - 1
drop if leaver & year > 2014
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()

cap noi wooldid y dist year onset, cluster(dist) att esprelength(2) espostlength(2) ///
    permute(shuffle) permreps(100) rseed(99) permsave(pF_save) replace
dis "PF_RC=" _rc
dis "PF_reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
local NFAIL = e(perm_failed)

use pF_save.dta, clear
dis "PF_ROWS=" _N

* capture BOTH variable lists BEFORE generating anything (a generated flag would
* otherwise be swept into the second ds)
qui ds _rep, not
local allvars `r(varlist)'
qui ds _rep onset_*, not
local coefvars `r(varlist)'
dis "PF_ALLCOLS = `: word count `allvars''   PF_COEFCOLS = `: word count `coefvars''"

* (1) the reading the help invites: a failed replication's row is ALL MISSING
qui gen byte allmiss_helpreading = 1
foreach v of local allvars {
  qui replace allmiss_helpreading = 0 if !missing(`v')
}
qui count if _rep >= 1 & allmiss_helpreading == 1
dis "PF_FAILED_BY_HELP_READING = " r(N) "   (e(perm_failed) was `NFAIL')"

* (2) the reading that actually works: only the COEFFICIENT columns are missing
qui gen byte allmiss_coefonly = 1
foreach v of local coefvars {
  qui replace allmiss_coefonly = 0 if !missing(`v')
}
qui count if _rep >= 1 & allmiss_coefonly == 1
dis "PF_FAILED_BY_COEF_COLUMNS = " r(N)

* (3) show one failed row in full: coefficients missing, onsets carry the assignment
qui levelsof _rep if _rep >= 1 & allmiss_coefonly == 1, local(fr)
local f1 : word 1 of `fr'
dis "PF_FAILED_REPS = `fr'"
list _rep `coefvars' if _rep == `f1', noobs abbrev(24)
list _rep onset_* if _rep == `f1', noobs abbrev(24)
qui egen byte nonmiss_onsets = rownonmiss(onset_*)
qui sum nonmiss_onsets if _rep == `f1', meanonly
dis "PF_FAILEDROW_NONMISSING_ONSET_CELLS = " r(mean) " of 10"
dis "PF_NOTE: a missing onset_ cell is the never-treated code, not an absent value, so the"
dis "PF_NOTE: assignment on a failed row is complete; at least the treated clusters' cells"
dis "PF_NOTE: carry a value, which is why the row is not 'all missing'."

dis "PF_DONE"
