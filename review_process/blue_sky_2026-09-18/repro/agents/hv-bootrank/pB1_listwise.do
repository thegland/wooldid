*=======================================================================
* pB1_listwise.do -- MY OWN lead.  e(V) is the LISTWISE complete-case
* covariance over all non-excluded columns (wooldid_bootfin,
* ado:14536-14552), so it rests on e(boot_N_V), NOT on reps().  The
* per-column standard errors use a different, larger row set, the
* low-availability caution (ado:10296) keys on per-column counts at a
* 50% threshold, and the footer (ado:11786) prints the TOKEN
* "e(boot_N_V)" rather than its value.  Measured here on a realistic
* 30-state staggered-adoption panel with a long event-study window.
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
qui levelsof state, local(ss)
dis "NSTATES=" `: word count `ss''

cap prog drop vdump
prog def vdump
  args tag
  dis "`tag'_COMPLETED=" e(boot_completed)
  dis "`tag'_FAILED=" e(boot_failed)
  dis "`tag'_N_V=" e(boot_N_V)
  dis "`tag'_REPS=" e(boot_reps)
  dis "`tag'_VPOSTED=" e(boot_v_posted)
  dis "`tag'_DFR=" e(df_r)
  dis "`tag'_K=" colsof(e(b))
  dis "`tag'_VESFAILED=" e(boot_ves_failed)
  dis "`tag'_JTABSENT=" e(boot_jt_absent)
  dis "`tag'_SGJTABSENT=" e(boot_sgjt_absent)
  dis "`tag'_SGJTREASON=" e(boot_sgjt_reason)
  mata: V = st_matrix("e(V)")
  mata: d = diagonal(V)
  mata: printf("`tag'_RANK=%f\n", rank(V))
  mata: printf("`tag'_NSENT=%f\n", sum(d :< -1e99))
  mata: printf("`tag'_NZERODIAG=%f\n", sum(d :== 0))
  mata: printf("`tag'_NPOSDIAG=%f\n", sum(d :> 0))
  mata: ev = symeigenvalues(makesymmetric(V))
  mata: printf("`tag'_MINEIG=%g MAXEIG=%g\n", min(ev), max(ev))
  cap confirm matrix e(V_es_att)
  if _rc == 0 {
    mata: printf("`tag'_RANK_VESATT=%f of %f\n", rank(st_matrix("e(V_es_att)")), cols(st_matrix("e(V_es_att)")))
  }
  * per-column availability: the ONLY disclosure surface, and its caution
  * threshold is 50% of completed reps
  cap confirm matrix e(boot_navail)
  if _rc == 0 {
    mata: na = st_matrix("e(boot_navail)")
    mata: printf("`tag'_NAVAIL_MIN=%f MAX=%f MEAN=%f\n", min(na), max(na), mean(na'))
    mata: printf("`tag'_NCOL_BELOW_HALF=%f\n", sum(na :< 0.5*st_numscalar("e(boot_completed)")))
  }
  * recompute the listwise count MYSELF from e(bootdist), and the
  * ES-family-only listwise count, so both B's are on the record
  cap confirm matrix e(bootdist)
  if _rc == 0 {
    mata: BD = st_matrix("e(bootdist)")
    mata: cs = st_matrixcolstripe("e(bootdist)")
    mata: printf("`tag'_BD_ROWS=%f COLS=%f\n", rows(BD), cols(BD))
    mata: printf("`tag'_BD_LISTWISE_ALL=%f\n", sum(rowmissing(BD) :== 0))
    mata: cst = st_matrixcolstripe("e(bootdist)")
    mata: nm = cst[., 2]
    mata: esel = selectindex(strpos(nm, "ES_att_") :== 1)
    mata: if (length(esel) > 0) printf("`tag'_BD_LISTWISE_ESATT=%f over %f cols\n", sum(rowmissing(BD[., esel]) :== 0), length(esel))
  }
end

*---------------------------------------------------------------
* SPEC 1: a long-window state-policy event study, NO subgroups
*---------------------------------------------------------------
dis "=== S1_BEGIN ==="
cap noisily wooldid y state t ttre, att att_it esfixedbaseperiod ///
    esprelength(6) espostlength(11) cluster(state) jointtests ///
    pretrendpower(50 80) bootstrap(standard) reps(60) rseed(4242)
dis "S1_RC=" _rc
if _rc == 0 {
  vdump S1
  mata: st_matrixcolstripe("e(b)")
  cap noisily mat list e(boot_navail), format(%6.0f)
}
dis "=== S1_END ==="

dis "=== S1TEST_BEGIN ==="
cap noisily test ES_att_pre2 ES_att_pre3 ES_att_pre4 ES_att_pre5 ES_att_pre6 ES_att_pre7
dis "S1TEST_RC=" _rc
dis "S1TEST_F=" r(F)
dis "S1TEST_chi2=" r(chi2)
dis "S1TEST_df=" r(df)
dis "S1TEST_dfr=" r(df_r)
dis "S1TEST_p=" r(p)
dis "=== S1TEST_END ==="

dis "=== S1PTP_BEGIN ==="
foreach m in wdidpretrendpower pretrendpower_att wdidptp {
  cap noisily mat list e(`m'), format(%12.6f)
}
dis "=== S1PTP_END ==="

*---------------------------------------------------------------
* SPEC 3: narrow window control - should NOT collapse
*---------------------------------------------------------------
dis "=== S3_BEGIN ==="
cap noisily wooldid y state t ttre, att att_it esfixedbaseperiod ///
    esprelength(2) espostlength(2) cluster(state) ///
    bootstrap(standard) reps(60) rseed(4242)
dis "S3_RC=" _rc
if _rc == 0 vdump S3
dis "=== S3_END ==="

dis "PROBE_B1_DONE"
