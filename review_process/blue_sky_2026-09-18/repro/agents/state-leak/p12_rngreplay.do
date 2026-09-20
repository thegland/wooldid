* p12: the documented reproducibility handle -- sthlp:542 says a run without
* rseed() "can still be replayed exactly by set rngstate to that value before
* repeating the call".  Measure it, for permute() and for bootstrap().
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

cap prog drop wdmain
prog def wdmain
  args tag
  capture confirm matrix e(wdidmainresults)
  if _rc != 0 {
    dis "MR|`tag'|ABSENT"
    exit
  }
  tempname M
  matrix `M' = e(wdidmainresults)
  local rn : rownames `M'
  local cn : colnames `M'
  forvalues a = 1/`=rowsof(`M')' {
    local rr : word `a' of `rn'
    forvalues b = 1/`=colsof(`M')' {
      local cc2 : word `b' of `cn'
      dis "MR|`tag'|`rr'|`cc2'=" %21x `M'[`a',`b']
    }
  }
end

* ============ PERMUTE: run with no rseed, capture e(perm_rngstate), replay
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(50) dots(0)
dis "P1_RC=" _rc
global WDH_PRNG "`e(perm_rngstate)'"
dis "P1_RNGLEN=" strlen("$WDH_PRNG")
wdmain P1_original

wdgen
* NOTE: set rngstate takes an UNQUOTED state token; quoting it is r(198)
capture noisily set rngstate $WDH_PRNG
dis "P2_SETRNG_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(50) dots(0)
dis "P2_RC=" _rc
wdmain P2_replayed

* and the quoted spelling a user would naturally write
wdgen
capture noisily set rngstate "$WDH_PRNG"
dis "P3_QUOTED_SETRNG_RC=" _rc

* ============ BOOTSTRAP: same claim
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(50) dots(0)
dis "B1_RC=" _rc
global WDH_BRNG "`e(boot_rngstate)'"
dis "B1_RNGLEN=" strlen("$WDH_BRNG")
wdmain B1_original

wdgen
capture noisily set rngstate $WDH_BRNG
dis "B2_SETRNG_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(50) dots(0)
dis "B2_RC=" _rc
wdmain B2_replayed

* ============ and the rseed() path, which must be seed-only reproducible
wdgen
set seed 101
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(50) dots(0) rseed(2024)
dis "R1_RC=" _rc
wdmain R1_rseed_from101
wdgen
set seed 202
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(50) dots(0) rseed(2024)
dis "R2_RC=" _rc
wdmain R2_rseed_from202

dis "P12_DONE"
