*=======================================================================
* r3_big.do -- THE DECISIVE RUN.  Same realistic spec as the pilot, but
* reps(300), so the listwise complete-case count is ~60 instead of ~6-14.
* At N_lw = 60 the sampling noise in an SD is ~9%, so any systematic
* downward bias of the exported covariance has to show up as a median
* ratio well below 1.  The block study then reproduces the B=60 regime
* the original finding was measured in, from the same 300 draws.
* SELF-CONTAINED: needs gen_skep.do and anal.mata beside it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/gen_skep.do"
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/anal.mata"

cap noisily wooldid y state t ttre, att esfixedbaseperiod esprelength(5) ///
    espostlength(10) cluster(state) subgroups(reg) bootstrap(standard) ///
    reps(300) rseed(9001)
dis "BIG_RC=" _rc
dis "BIG_COMPLETED=" e(boot_completed)
dis "BIG_FAILED=" e(boot_failed)
dis "BIG_N_V=" e(boot_N_V)
dis "BIG_VPOSTED=" e(boot_v_posted)
dis "BIG_DFR=" e(df_r)
tempname BD VV NAV
matrix `BD' = e(bootdist)
matrix `VV' = e(V)
matrix `NAV' = e(boot_navail)
cap noisily lincom main_att
dis "BIG_LINCOM_se=" r(se)
mata: wd_analyze("`BD'", "`VV'", "`NAV'", 60, 1000, "BIG")
dis "R3_BIG_DONE"
