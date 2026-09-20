* p18: (a) does the replication line re-evaluate the user's `if' expression on the
* RESAMPLED frame, where `_n' means something else?  (b) do two DIFFERENT engine
* runs back to back leave each other's e() behind?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

cap prog drop wdse
prog def wdse
  args tag
  capture confirm matrix e(wdidmainresults)
  if _rc != 0 {
    dis "SE|`tag'|ABSENT"
    exit
  }
  tempname M
  matrix `M' = e(wdidmainresults)
  dis "SE|`tag'|N=" e(N) "|est=" %18.12g `M'[1,1] "|se=" %18.12g `M'[1,2] "|lb=" %18.12g `M'[1,5] "|ub=" %18.12g `M'[1,6]
  dis "SE|`tag'|completed=" e(boot_completed) "|failed=" e(boot_failed) "|permn=" e(perm_n)
end

* ============ A: identical SAMPLES, two spellings of the same restriction.
* keepflag is the user's row order; `_n <= 360' is the same 360 rows on the
* observed run (p08 proved that).  Under bootstrap the replications run on a
* RESAMPLED frame, where `_n <= 360' is a different, arbitrary set of rows.
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if keepflag == 1, cluster(state) att bootstrap(standard) reps(60) dots(0) rseed(11)
dis "A1_RC=" _rc
wdse A1_flagvar_boot

wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if _n <= 360, cluster(state) att bootstrap(standard) reps(60) dots(0) rseed(11)
dis "A2_RC=" _rc
wdse A2_underscore_n_boot

* the same pair under permute()
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if keepflag == 1, cluster(state) att permute(shuffle) permreps(60) dots(0) rseed(11)
dis "A3_RC=" _rc
wdse A3_flagvar_perm
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if _n <= 360, cluster(state) att permute(shuffle) permreps(60) dots(0) rseed(11)
dis "A4_RC=" _rc
wdse A4_underscore_n_perm

* ============ B: a NONDETERMINISTIC `if' under an engine
wdgen
set seed 31415
qui gen double uu = runiform()
qui gen byte uflag = (uu < 0.7)
capture noisily wooldid y unit year ttre if uflag == 1, cluster(state) att bootstrap(standard) reps(60) dots(0) rseed(11)
dis "B1_RC=" _rc
wdse B1_uflag_boot
wdgen
set seed 31415
capture noisily wooldid y unit year ttre if runiform() < 0.7, cluster(state) att bootstrap(standard) reps(60) dots(0) rseed(11)
dis "B2_RC=" _rc
wdse B2_runiform_boot

* ============ C: two DIFFERENT engines back to back on the same data
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(40) dots(0) rseed(11)
dis "C1_RC=" _rc
dis "C1_perm_method=" "`e(perm_method)'" "|boot_method=" "`e(boot_method)'"
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(11)
dis "C2_RC=" _rc
dis "C2_perm_method=" "`e(perm_method)'" "|boot_method=" "`e(boot_method)'"
capture confirm matrix e(perm_p_lo)
dis "C2_STALE_perm_p_lo=" (_rc == 0)
capture confirm matrix e(bootdist)
dis "C2_bootdist=" (_rc == 0)
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(40) dots(0) rseed(11)
dis "C3_RC=" _rc
dis "C3_perm_method=" "`e(perm_method)'" "|boot_method=" "`e(boot_method)'"
capture confirm matrix e(bootdist)
dis "C3_STALE_bootdist=" (_rc == 0)
capture confirm matrix e(V_boot)
dis "C3_STALE_V_boot=" (_rc == 0)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "C4_RC=" _rc
dis "C4_perm_method=" "`e(perm_method)'" "|boot_method=" "`e(boot_method)'"
capture confirm matrix e(perm_p_lo)
dis "C4_STALE_perm=" (_rc == 0)
capture confirm matrix e(V_boot)
dis "C4_STALE_V_boot=" (_rc == 0)
wddumpe engines_then_plain

* ============ D: `if e(sample)' -- a very common applied pattern ("estimate on
* the same sample my first stage used").  e() is SESSION state, and wooldid
* overwrites it with its own post partway through the call, then hands the raw
* expression to every replication.
wdgen
qui regress y post if state <= 8
qui count if e(sample)
dis "D0_PRIOR_SAMPLE_ROWS=" r(N)
capture noisily wooldid y unit year ttre if e(sample), cluster(state) att
dis "D1_RC=" _rc
wdse D1_plain_ifesample

wdgen
qui regress y post if state <= 8
qui gen byte esflag = e(sample)
qui count if esflag == 1
dis "D2_FLAG_ROWS=" r(N)
capture noisily wooldid y unit year ttre if esflag == 1, cluster(state) att
dis "D2_RC=" _rc
wdse D2_plain_flag

wdgen
qui regress y post if state <= 8
capture noisily wooldid y unit year ttre if e(sample), cluster(state) att bootstrap(standard) reps(60) dots(0) rseed(11)
dis "D3_RC=" _rc
wdse D3_boot_ifesample

wdgen
qui regress y post if state <= 8
qui gen byte esflag = e(sample)
capture noisily wooldid y unit year ttre if esflag == 1, cluster(state) att bootstrap(standard) reps(60) dots(0) rseed(11)
dis "D4_RC=" _rc
wdse D4_boot_flag

wdgen
qui regress y post if state <= 8
capture noisily wooldid y unit year ttre if e(sample), cluster(state) att permute(shuffle) permreps(60) dots(0) rseed(11)
dis "D5_RC=" _rc
wdse D5_perm_ifesample

wdgen
qui regress y post if state <= 8
qui gen byte esflag = e(sample)
capture noisily wooldid y unit year ttre if esflag == 1, cluster(state) att permute(shuffle) permreps(60) dots(0) rseed(11)
dis "D6_RC=" _rc
wdse D6_perm_flag

dis "P18_DONE"
