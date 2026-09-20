*=======================================================================
* r4_collapse.do -- can the listwise rule drive e(boot_N_V) to 0 or 1?
* wooldid_bootfin only fills V when >= 2 listwise-complete rows exist;
* otherwise V stays all-zero and EVERY column is given the -1e100
* sentinel.  e(boot_v_posted) is still 1 and the footer still says
* "e(V) holds the bootstrap covariance".  Realistic trigger: the same
* 36-state panel with SIX subgroups (think census divisions) instead of
* three, at reps(60).
* SELF-CONTAINED: needs gen_skep.do beside it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/gen_skep.do"
gen byte div = mod(state - 1, 6) + 1
label define divl 1 "NewEng" 2 "MidAtl" 3 "SouthAtl" 4 "ENCentral" 5 "WNCentral" 6 "Pacific"
label values div divl

cap noisily wooldid y state t ttre, att esfixedbaseperiod esprelength(5) ///
    espostlength(10) cluster(state) subgroups(div) bootstrap(standard) ///
    reps(60) rseed(7701)
dis "C_RC=" _rc
dis "C_COMPLETED=" e(boot_completed)
dis "C_FAILED=" e(boot_failed)
dis "C_N_V=" e(boot_N_V)
dis "C_VPOSTED=" e(boot_v_posted)
dis "C_DFR=" e(df_r)
tempname VV
matrix `VV' = e(V)
dis "C_K=" colsof(`VV')
mata:
  V = st_matrix("`VV'")
  d = diagonal(V)'
  printf("C_DIAG_SENTINEL=%f of %f\n", sum(d :== -1e100), cols(d))
  printf("C_DIAG_POSITIVE=%f\n", sum(d :> 0 :& d :< .))
  printf("C_ALLZERO_OFFDIAG=%f\n", all((V - diag(diagonal(V))) :== 0))
end
dis "=== what the consumers do with an all-sentinel e(V) ==="
cap noisily lincom main_att
dis "C_LINCOM_RC=" _rc " se=" r(se) " b=" r(estimate) " lb=" r(lb) " ub=" r(ub)
cap noisily test ES_att_pre2 ES_att_pre3
dis "C_TEST_RC=" _rc " df=" r(df) " p=" r(p)
cap noisily estimates store cm
cap noisily esttab cm, se b(%12.6f) se(%12.6f) nostar nonumbers keep(main_att ES_att_post0 ES_att_pre2)
dis "C_ESTTAB_RC=" _rc
dis "R4_COLLAPSE_DONE"
