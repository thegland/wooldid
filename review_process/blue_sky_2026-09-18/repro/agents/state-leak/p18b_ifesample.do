* p18b: the high-value half of p18, with the wdse type-mismatch fixed.
*  D: `if e(sample)' -- e() is SESSION state, and ado:9057 hands the RAW `if'
*     string to every nested replication, which re-evaluates it in a frame where
*     e() has since been overwritten by the run itself.
*  A: `if _n <= 360' vs the same rows as a flag variable, under permute().
*  C: two DIFFERENT engines back to back -- stale e() between them.
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
  dis "SE|`tag'|N=" e(N)
  dis "SE|`tag'|est=" %18.12g `M'[1,1]
  dis "SE|`tag'|se=" %18.12g `M'[1,2]
  dis "SE|`tag'|p=" %18.12g `M'[1,4]
  dis "SE|`tag'|lb=" %18.12g `M'[1,5]
  dis "SE|`tag'|ub=" %18.12g `M'[1,6]
  capture dis "SE|`tag'|bootdone=" e(boot_completed) "|bootfail=" e(boot_failed)
end

* ============ D: `if e(sample)'
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
capture noisily wooldid y unit year ttre if e(sample), cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(11)
dis "D3_RC=" _rc
wdse D3_boot_ifesample

wdgen
qui regress y post if state <= 8
qui gen byte esflag = e(sample)
capture noisily wooldid y unit year ttre if esflag == 1, cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(11)
dis "D4_RC=" _rc
wdse D4_boot_flag

* ============ A: `_n <= 360' vs a flag, under permute()
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if keepflag == 1, cluster(state) att permute(shuffle) permreps(40) dots(0) rseed(11)
dis "A3_RC=" _rc
wdse A3_flagvar_perm
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if _n <= 360, cluster(state) att permute(shuffle) permreps(40) dots(0) rseed(11)
dis "A4_RC=" _rc
wdse A4_underscore_n_perm

* ============ C: two DIFFERENT engines back to back
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(25) dots(0) rseed(11)
dis "C1_RC=" _rc
dis "C1_perm_method=" "`e(perm_method)'" "|boot_method=" "`e(boot_method)'"
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(25) dots(0) rseed(11)
dis "C2_RC=" _rc
dis "C2_perm_method=" "`e(perm_method)'" "|boot_method=" "`e(boot_method)'"
capture confirm matrix e(perm_p_lo)
dis "C2_STALE_perm_p_lo=" (_rc == 0)
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(25) dots(0) rseed(11)
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

dis "P18B_DONE"
