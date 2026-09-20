*=======================================================================
* pE_sediverge.do -- the CONSEQUENCE of the listwise collapse that a
* user can actually see: wooldid's PRINTED standard error for a
* coefficient is the sd over that column's own available replications
* (up to 60 here), while sqrt(diag(e(V))) -- what esttab, coefplot,
* honestdid, test and lincom all read -- is the sd over the LISTWISE
* complete-case rows only (e(boot_N_V)).  When those two counts are far
* apart the two standard errors for the SAME coefficient are far apart,
* and nothing in the run says so.
*
* Same 30-state panel and long window as pB2, `att' only so the
* replications use frequency reweighting and the run is quick.
* SELF-CONTAINED: needs gen_admin.do next to it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank/gen_admin.do"

dis "=== E1_BEGIN ==="
cap noisily wooldid y state t ttre, att esfixedbaseperiod ///
    esprelength(6) espostlength(11) cluster(state) subgroups(reg) ///
    bootstrap(standard) reps(60) rseed(4242)
dis "E1_RC=" _rc
dis "E1_COMPLETED=" e(boot_completed)
dis "E1_N_V=" e(boot_N_V)
dis "E1_K=" colsof(e(b))
dis "=== E1_END ==="

* how many cautions did the run print?  (grep the log for "Caution")

mata:
V   = st_matrix("e(V)")
na  = st_matrix("e(boot_navail)")
cs  = st_matrixcolstripe("e(V)")
comp = st_numscalar("e(boot_completed)")
nv   = st_numscalar("e(boot_N_V)")
printf("E1 rank(V)=%f of %f   boot_N_V=%f  completed=%f\n", rank(V), cols(V), nv, comp)
printf("E1 navail min=%f max=%f\n", min(na), max(na))
printf("E1 columns below half of completed: %f\n", sum(na :< 0.5*comp))

// the printed table: e(wdidmainresults) carries b and se as displayed
W  = st_matrix("e(wdidmainresults)")
ws = st_matrixrowstripe("e(wdidmainresults)")
wc = st_matrixcolstripe("e(wdidmainresults)")
printf("wdidmainresults dims %f x %f; columns:\n", rows(W), cols(W))
wc[.,2]'
end

mat list e(wdidmainresults), format(%12.6f)

mata:
// sqrt of the e(V) diagonal, by e(b) name
d  = diagonal(st_matrix("e(V)"))
nm = st_matrixcolstripe("e(V)")
printf("\n%-34s %14s %14s\n", "coefficient", "sqrt(diag eV)", "navail")
for (i = 1; i <= rows(nm); i++) {
    printf("SEDIAG %-30s %14.6f %10.0f\n", nm[i,2], (d[i] > 0 ? sqrt(d[i]) : .), st_matrix("e(boot_navail)")[1,i])
}
end

dis "PROBE_E_DONE"
