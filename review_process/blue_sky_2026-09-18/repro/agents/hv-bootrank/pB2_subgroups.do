*=======================================================================
* pB2_subgroups.do -- companion to pB1: the same long-window state
* event study WITH regional subgroups, which is where the listwise
* intersection actually collapses.  SELF-CONTAINED.
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
* SPEC 2: the SAME long window with regional subgroups.
* reg = mod(state-1,4)+1, so the four earliest adopters (states
* 1,2,3,4) sit in four DIFFERENT regions.  Each region's deepest
* post horizon is therefore carried by ONE state, and those four
* dropout events are independent across resamples: the listwise
* complete-case share should be about 0.638^4 ~ 0.17.
*---------------------------------------------------------------
dis "=== S2_BEGIN ==="
cap noisily wooldid y state t ttre, att att_it esfixedbaseperiod ///
    esprelength(6) espostlength(11) cluster(state) subgroups(reg) ///
    bootstrap(standard) reps(60) rseed(4242)
dis "S2_RC=" _rc
if _rc == 0 {
  vdump S2
  cap noisily mat list e(boot_navail), format(%6.0f)
}
dis "=== S2_END ==="

dis "=== S2TEST_BEGIN ==="
cap noisily test ES_att_pre2 ES_att_pre3 ES_att_pre4 ES_att_pre5 ES_att_pre6 ES_att_pre7
dis "S2TEST_RC=" _rc
dis "S2TEST_F=" r(F)
dis "S2TEST_df=" r(df)
dis "S2TEST_dfr=" r(df_r)
dis "S2TEST_p=" r(p)
dis "=== S2TEST_END ==="

dis "PROBE_B2_DONE"
