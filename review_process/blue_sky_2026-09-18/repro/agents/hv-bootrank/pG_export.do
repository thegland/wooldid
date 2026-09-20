*=======================================================================
* pG_export.do -- the consequence, demonstrated rather than argued.
* The run's own footer says "esttab, coefplot and honestdid export the
* standard-deviation family".  This shows what esttab and lincom
* actually print next to what wooldid printed, on the same estimate.
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

cap noisily wooldid y state t ttre, att esfixedbaseperiod ///
    esprelength(6) espostlength(11) cluster(state) subgroups(reg) ///
    bootstrap(standard) reps(60) rseed(4242)
dis "G_RC=" _rc
dis "G_COMPLETED=" e(boot_completed)
dis "G_N_V=" e(boot_N_V)
dis "G_DFR=" e(df_r)

dis "=== WHAT WOOLDID PRINTED (first rows of e(wdidmainresults)) ==="
cap noisily mat list e(wdidmainresults), format(%12.6f)

dis "=== WHAT lincom SEES ==="
cap noisily lincom main_att
dis "G_LINCOM_RC=" _rc
dis "G_LINCOM_se=" r(se)
dis "G_LINCOM_est=" r(estimate)
dis "G_LINCOM_lb=" r(lb)
dis "G_LINCOM_ub=" r(ub)
cap noisily lincom ES_att_post9
dis "G_LINCOM9_se=" r(se)

dis "=== WHAT test SEES (6 pre-period coefficients, a TRUE null here) ==="
cap noisily test ES_att_pre2 ES_att_pre3 ES_att_pre4 ES_att_pre5 ES_att_pre6 ES_att_pre7
dis "G_TEST_RC=" _rc
dis "G_TEST_F=" r(F)
dis "G_TEST_df=" r(df)
dis "G_TEST_dfr=" r(df_r)
dis "G_TEST_p=" r(p)

dis "=== WHAT esttab EXPORTS ==="
cap noisily estimates store wdm
dis "G_STORE_RC=" _rc
cap noisily esttab wdm, se b(%12.6f) se(%12.6f) nostar keep(main_att ES_att_post0 ES_att_post9 ES_att_post10 ES_att_pre2) nonumbers
dis "G_ESTTAB_RC=" _rc

dis "PROBE_G_DONE"
