* p07: bootstrap() and the session RNG -- seeding, reproducibility, leftovers
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

cap prog drop wdboot
prog def wdboot
  args tag
  dis "BT|`tag'|method=" "`e(boot_method)'"
  dis "BT|`tag'|seed=" "`e(boot_seed)'"
  dis "BT|`tag'|rng_head=" substr("`e(boot_rngstate)'",1,24)
  dis "BT|`tag'|completed=" e(boot_completed) "|failed=" e(boot_failed)
  capture confirm matrix e(wdidmainresults)
  if _rc == 0 {
    tempname M
    matrix `M' = e(wdidmainresults)
    local rn : rownames `M'
    local cn : colnames `M'
    forvalues a = 1/`=rowsof(`M')' {
      local rr : word `a' of `rn'
      forvalues b = 1/`=colsof(`M')' {
        local cc2 : word `b' of `cn'
        dis "BT|`tag'|`rr'|`cc2'=" %21x `M'[`a',`b']
      }
    }
  }
end

* ============ ARM A: bootstrap, no rseed, user set seed 777
wdgen
set seed 777
global WDH_R1 "`c(rngstate)'"
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0)
dis "A_RC=" _rc
global WDH_R2 "`c(rngstate)'"
dis "A_RNG_ADVANCED=" ("$WDH_R1" != "$WDH_R2")
global WDH_BOOTRNG "`e(boot_rngstate)'"
dis "A_BOOTRNG_EQ_PRECALL=" ("$WDH_BOOTRNG" == "$WDH_R1")
wdboot A_seed777
wdline A_post

* ============ ARM B: same call after set seed 777 -- must reproduce exactly
wdgen
set seed 777
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0)
dis "B_RC=" _rc
wdboot B_seed777_again

* ============ ARM C: the documented replay handle -- set rngstate to e(boot_rngstate)
wdgen
capture noisily set rngstate "$WDH_BOOTRNG"
dis "C_SETRNG_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0)
dis "C_RC=" _rc
wdboot C_rngstate_replay

* ============ ARM D: rseed() from two different starting streams -- must match
wdgen
set seed 111
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(4242)
dis "D_RC=" _rc
wdboot D_rseed_from111
wdgen
set seed 999
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(4242)
dis "D2_RC=" _rc
wdboot D2_rseed_from999

* ============ ARM E: TWO bootstrap runs back to back with NO reseeding.
* The observed point estimates must be identical between them.
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0)
dis "E1_RC=" _rc
wdboot E1_first
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0)
dis "E2_RC=" _rc
wdboot E2_second_sameData

* ============ ARM F: bayes mode + bootsave file characteristics
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(bayes) reps(40) dots(0) rseed(4242) bootsave(p07boot) replace safety(off)
dis "F_RC=" _rc
wdline F_post
wdboot F_bayes

* ============ ARM G: back to a plain analytic run -- does anything stale survive?
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "G_RC=" _rc
dis "G_boot_method=" "`e(boot_method)'"
dis "G_perm_method=" "`e(perm_method)'"
wddumpe boot_then_plain
wdline G_post

dis "P07_DONE"
