* NOTE (post-R1, 2026-08-24): wooldid now DETECTS this design (the cell-containment
* predicate), prints an explanatory caution, sets the ame SEs/tests to missing, and
* withholds e(V_es_ame) while still posting e(b_es_ame). Exhibit A below therefore no
* longer shows the silent identical SEs on current builds — it shows the refusal.
* To reproduce the original pathological OUTPUT, run Exhibit A against a pre-R1 build
* (e.g. git show <pre-R1 sha>:wooldidinstall/wooldid.ado). Exhibit B is unchanged.
*
* demo_vesame_pathology.do — minimal demonstration of NF-2 (wooldid_old_bug_report.md
* §5.3): under cohort-level cluster(i) with contreat(), EVERY reported ame standard
* error (main ame, pretreatment ame, every ES_ame horizon) is the same number — the SE
* of the absorbed baseline dose slope — and e(V_es_ame) is EXACTLY rank 1, so every
* cross-horizon contrast carries zero reported sampling variance even though the point
* estimates differ. The true sampling covariance is full rank with per-horizon
* variances ~5x the reported common value (120-rep MC; mechanism and old-build
* comparison: review_process/inv1_v_es_ame_2026-08-24.md).
*
* Run from the repo root:   "<stata>" -e do review_process/demo_vesame_pathology.do
* (or pass the wooldid.ado path as argument 1).
* Everything is deterministic (fixed seed). Runtime ~15 seconds.

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if `"`ADO'"' == "" local ADO wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'

* ---- a well-behaved continuous-treatment panel: 30 cohorts (i) x 6 units x 12
* periods; dose varies IRREGULARLY across units within every cohort-period cell
* (so nothing here is a balanced-design artifact); half the cohorts treated at t=6;
* y contains NO dose effect and NO treatment effect whatsoever ----
set seed 2010
qui set obs 30
gen i = _n
gen ttre = cond(mod(_n,2)==0, 6, .)
qui expand 6
bysort i: gen unit = _n
gen dose = exp(0.8*runiform() + 0.2*mod(i,5))
qui expand 12
bysort i unit: gen t = _n
gen y = 0.1*i + 0.05*t + rnormal()
sort i unit t
tempfile demo
save `demo'

* ================= Exhibit A: the pathological run (cohort-level clustering) =====
* Look at the printed table: main ame, pretreatment ame, and ALL eight ES_ame rows
* carry the IDENTICAL standard error, while their point estimates differ freely.
use `demo', clear
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod

tempname V B
matrix `B' = e(b_es_ame)
matrix `V' = e(V_es_ame)
dis _n as result "---- Exhibit A diagnostics: e(V_es_ame) under cluster(i) ----"
dis as text "per-horizon reported SEs (sqrt of diag):"
forval j = 1/8 {
  local l : word `j' of `: colnames `B''
  dis as text "  `l'   b = " %10.6f `B'[1,`j'] "   se = " %12.9f sqrt(`V'[`j',`j'])
}
mata: V = st_matrix(st_local("V")); C = corr(V); printf("min pairwise correlation = %16.14g (all exactly 1)\n", min(C))
mata: ev = symeigenvalues(V); printf("eigenvalues: largest = %10.6g, second-largest = %10.3g  (rank 1)\n", max(ev), sort(ev',1)[7])
dis as text "reported variance of the (post0 - pre2) CONTRAST (V11+V22-2V12 form):"
dis as result "  " %12.4e (`V'[5,5] + `V'[4,4] - 2*`V'[5,4]) as text "   while the estimates differ by " %10.6f (`B'[1,5] - `B'[1,4])
dis as text "=> every cross-horizon comparison is reported as having ZERO sampling"
dis as text "   variance; honestdid/lincom/test consumers of e(V_es_ame) ingest this"
dis as text "   silently. (True per-horizon variances are ~5x the common reported"
dis as text "   value and the true covariance is full rank - see the INV-1 report.)"

* ================= Exhibit B: same data, robust VCE — the collapse is gone ======
* (Full-rank V, horizon-specific SEs. NOTE: full-rank does not mean well-sized ame
* inference — see NF-3 in the report — but it isolates the clustering collapse.)
use `demo', clear
qui wooldid y i t ttre, att robust contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod
tempname V2 B2
matrix `B2' = e(b_es_ame)
matrix `V2' = e(V_es_ame)
dis _n as result "---- Exhibit B diagnostics: e(V_es_ame) under robust ----"
forval j = 1/8 {
  local l : word `j' of `: colnames `B2''
  dis as text "  `l'   b = " %10.6f `B2'[1,`j'] "   se = " %12.9f sqrt(`V2'[`j',`j'])
}
mata: V2 = st_matrix(st_local("V2")); ev2 = symeigenvalues(V2); printf("eigenvalues: largest = %10.6g, smallest = %10.6g  (full rank)\n", max(ev2), min(ev2))

dis _n as result "DEMO_DONE"
exit, clear
