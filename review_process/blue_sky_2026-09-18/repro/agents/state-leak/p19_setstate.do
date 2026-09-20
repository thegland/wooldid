* p19: SESSION SETTINGS a real user plausibly has switched on, and whether wooldid
* is invariant to them.  `set dp comma' in particular changes what strofreal()
* produces, and wooldid builds matrix names, cell keys and messages from numbers.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ BASELINE (default settings)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) subgroups(sg) jointtests
dis "BASE_RC=" _rc
wddumpe set_base
dis "BASE_COLNAMES=" "`: colnames e(wdidmainresults)'"

* ============ A: set dp comma  (European decimal separator)
set dp comma
dis "A_DP=" "`c(dp)'"
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) subgroups(sg) jointtests
dis "A_RC=" _rc
wddumpe set_dpcomma
dis "A_COLNAMES=" "`: colnames e(wdidmainresults)'"
dis "A_SGLIST=" "`e(wdid_sglist)'"
* engines under dp comma
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(20) dots(0) rseed(7)
dis "A2_PERM_RC=" _rc
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(20) dots(0) rseed(7)
dis "A3_BOOT_RC=" _rc
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) simultaneousbands pretrendpower(50 80)
dis "A4_PTP_RC=" _rc
set dp period
dis "A_DP_BACK=" "`c(dp)'"

* ============ B: set varabbrev off
set varabbrev off
dis "B_VARABBREV=" "`c(varabbrev)'"
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) subgroups(sg) jointtests
dis "B_RC=" _rc
wddumpe set_varabbrevoff
set varabbrev on

* ============ C: set type double (affects `gen' default type)
set type double
dis "C_TYPE=" "`c(type)'"
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) subgroups(sg) jointtests
dis "C_RC=" _rc
wddumpe set_typedouble
set type float

* ============ D: non-default display formats
set cformat %12.8f
set pformat %8.6f
set sformat %9.5f
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3)
dis "D_RC=" _rc
wddumpe set_cformat
set cformat
set pformat
set sformat

* ============ E: a different sortseed (sort tie-breaking)
set sortseed 999
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) subgroups(sg) jointtests
dis "E_RC=" _rc
wddumpe set_sortseed999
set sortseed 20260916

* ============ F: ACCUMULATION -- 30 runs in a loop, the way a user sweeping
* outcomes or specifications actually writes it.  Does anything pile up?
set sortseed 20260916
wdgen
wdsnap loop_pre
forvalues k = 1/30 {
  capture wooldid y unit year ttre, cluster(state) att
  if _rc != 0 dis "F_ITER`k'_RC=" _rc
}
dis "F_LOOP_DONE_RC=" _rc
wdsnap loop_post
wddumpe set_loop30
* and 10 runs with plots + frames + subgroups, the heaviest shape
wdgen
wdsnap loopB_pre
forvalues k = 1/10 {
  capture wooldid y unit year ttre, cluster(state) att subgroups(sg) esprelength(2) espostlength(2) makeplots
  if _rc != 0 dis "FB_ITER`k'_RC=" _rc
}
dis "FB_LOOP_DONE_RC=" _rc
wdsnap loopB_post
capture noisily graph dir
capture noisily frames dir
capture noisily matrix dir
capture noisily estimates dir

dis "P19_DONE"
