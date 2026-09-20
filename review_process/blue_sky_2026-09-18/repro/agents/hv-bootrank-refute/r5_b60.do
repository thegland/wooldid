*=======================================================================
* r5_b60.do -- the ORIGINAL REGIME on my own data: the identical spec at
* reps(60).  Two questions: (1) how many cautions does the run print?
* (2) what is the exported/printed SE ratio here, and in which direction?
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

* four independent seeds of the SAME spec - direction is the whole question
foreach s in 5101 5102 5103 5104 {
  cap noisily wooldid y state t ttre, att esfixedbaseperiod esprelength(5) ///
      espostlength(10) cluster(state) subgroups(reg) bootstrap(standard) ///
      reps(60) rseed(`s')
  dis "S`s'_RC=" _rc
  dis "S`s'_COMPLETED=" e(boot_completed)
  dis "S`s'_N_V=" e(boot_N_V)
  tempname BD VV NAV
  matrix `BD' = e(bootdist)
  matrix `VV' = e(V)
  matrix `NAV' = e(boot_navail)
  cap noisily lincom main_att
  dis "S`s'_LINCOM_se=" r(se)
  mata: wd_analyze("`BD'", "`VV'", "`NAV'", 0, 0, "S`s'")
}
dis "R5_B60_DONE"
