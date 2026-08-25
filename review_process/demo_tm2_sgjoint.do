* demo_tm2_sgjoint.do — isolates T-M2 / SI-10: the per-subgroup CONTINUOUS joint
* pre-trend p-values are stored as literal missing because the code reads
* r(p)[1,2] where the continuous overjoint r(p) is 1x1 (real p at [1,1]).
* Stata's lenient expression-context subscripting returns missing for the
* out-of-range read with rc 0, so nothing ever errors — e(joint_Choriz_s#) simply
* stores "." where genuinely significant p-values exist. The BINARY twins read
* [1,2] correctly because the binary overjoint r(p) is 1x2 (vacuous base column
* at [1,1], real p at [1,2]).
*
* NOTE: run this with a VCE that does not trip the Group A cell-containment screen
* (robust below) — under cohort-level cluster(i) the whole subgroup ame joint block
* is now deliberately skipped, so the T-M2 sites never execute there.
*
* Run from the repo root:   "<stata>" -e do review_process/demo_tm2_sgjoint.do
* (or pass the wooldid.ado path as argument 1). Deterministic; ~30 seconds.
*
* STATUS NOTE (2026-08-25): the author-ordered fix has landed — the two
* sg-continuous sites now read r(p)[1,1], so running this demo against the
* current build shows REAL ame p-values in Exhibit A. To reproduce the defect,
* pass a pre-fix build as argument 1 (e.g. stage main's wooldid.ado to a file).
* The oldnew suite pins the old-vs-new contrast via trig_tm2_sgjoint (ledger
* rows TM2-sg-p / TM2-ext-p): expected old p=., new real p.

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if `"`ADO'"' == "" local ADO wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'

* ---- a rich clean panel: 15 cohorts x 8 units x 16 periods, dose varying within
* cells, 3 subgroups with genuinely DIFFERENT pre-trend behavior so the subgroup
* joint pre-trend tests have something real to find ----
clear
set seed 41005
qui set obs 15
gen i = _n
gen ttre = cond(_n <= 4, 8, cond(_n <= 8, 10, cond(_n <= 12, 12, .)))
qui expand 8
bysort i: gen unit = _n
gen sg = mod(i, 3) + 1
gen dose = exp(0.5*runiform() + 0.15*mod(unit, 4))
qui expand 16
bysort i unit: gen t = _n
gen rel = t - ttre
* real dose effects post-treatment + subgroup-differentiated PRE-trends in dose
* cond() keeps never-treated rows nonmissing (rel is missing there, and Stata's
* 0*missing is missing, so bare (rel<0)*(rel+8) products would silently drop the
* entire control group)
gen y = 0.05*i + 0.04*t + 0.35*dose*(rel >= 0 & rel < .) ///
      + 0.06*(sg == 2)*dose*cond(rel < 0 & rel < ., rel + 8, 0) ///
      + 0.10*(sg == 3)*dose*cond(rel < 0 & rel < ., rel + 8, 0) ///
      + rnormal()
sort i unit t

* ================= Exhibit A: the symptom on a full wooldid run =================
* Binary att subgroup joint pre-trend tests print REAL p-values; the ame twins for
* the same subgroups print "." — inspect the "p-value for joint test that sg=..."
* console lines and the e() scalars dumped below.
wooldid y i t ttre, att robust contreat(dose) contreatcontrols(basic) esprelength(3) espostlength(3) subgroups(sg) jointtests

dis _n as result "---- Exhibit A: stored subgroup joint pre-trend p-values ----"
foreach sgc in 1 2 3 {
  dis as text "  sg=`sgc'   att: joint_horiz_s`sgc' = [" `"`e(joint_horiz_s`sgc')'"' "]   ame: joint_Choriz_s`sgc' = [" `"`e(joint_Choriz_s`sgc')'"' "]"
}
dis as text "The att column carries real p-values; the ame column is empty/missing"
dis as text "even though subgroups 2 and 3 have PLANTED dose pre-trends the joint"
dis as text "test should detect (and the instrumented evidence in"
dis as text "review_process/tm2_evidence/ shows the mid-run r(p)[1,1] values are"
dis as text "highly significant on exactly this class of data)."

* ================= Exhibit B: the mechanics, reproduced standalone ==============
* The margins overjoint r(p) geometry differs between binary and continuous
* estimands, and Stata's expression-context subscripting silently returns missing
* for an out-of-range read. Reproduced here on a plain regression so the indexing
* trap is visible without wooldid internals:
dis _n as result "---- Exhibit B: the r(p) geometry trap ----"
qui gen post = (rel >= 0 & rel < .)
qui reg y i.post##c.dose i.sg
qui margins, dydx(dose) over(sg) contrast(overjoint)
tempname RP
matrix `RP' = r(p)
dis as text "continuous (dydx) overjoint r(p) is " rowsof(`RP') "x" colsof(`RP') ":"
matrix list `RP'
dis as text "  r(p)[1,1] = " r(p)[1,1] "   <- the real p lives here"
dis as text "  r(p)[1,2] = " r(p)[1,2] "   <- out of range: silently missing, rc 0"
qui margins, dydx(post) over(sg) contrast(overjoint)
matrix `RP' = r(p)
dis as text "binary (dydx of an indicator) overjoint r(p) is " rowsof(`RP') "x" colsof(`RP') ":"
matrix list `RP'
dis as text "  For 1x2 geometry [1,2] is correct (vacuous base at [1,1]) - which is"
dis as text "  why the binary att code path reads [1,2] and works, and why copying"
dis as text "  that read onto the 1x1 continuous geometry stores missing instead of"
dis as text "  the real p. Fix (per the mandated protocol): [1,2] -> [1,1] at the two"
dis as text "  sg-continuous sites, then re-bless cont_sg.csv (sole pinning baseline)."

dis _n as result "DEMO_TM2_DONE"
exit, clear
