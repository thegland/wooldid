* probe02: hand-recompute every enumerated replication and compare bitwise
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* harvest assignments
use "permA.dta", clear
local R = _N - 1
forval r = 0/`R' {
  local rr = `r' + 1
  forval s = 1/5 {
    local a`r'_`s' = onset_`s'[`rr']
  }
  local saved`r' = main_att[`rr']
  local savedx`r' : di %21x main_att[`rr']
}

tempname memhold
postfile `memhold' int rep double sav double extA double extB double extC using "probe02res.dta", replace

forval r = 0/`R' {
  use "dataA.dta", clear
  gen double repons = .
  forval s = 1/5 {
    qui replace repons = `a`r'_`s'' if state == `s'
  }
  * (A) observed-code map, exactly as the ado does: onset4->1, onset6->2, missing->3
  gen double repiA = cond(repons==4,1,cond(repons==6,2,3))
  * (B) fresh codes from egen group
  egen double repiB = group(repons), missing
  * (C) observed-code map, relabelled (10/20/30)
  gen double repiC = cond(repons==4,10,cond(repons==6,20,30))
  local vA = .
  local vB = .
  local vC = .
  cap qui wooldid y repiA t repons, att fastnose cluster(state) esrelativeto(-1)
  if _rc == 0 local vA = _b[main_att]
  cap qui wooldid y repiB t repons, att fastnose cluster(state) esrelativeto(-1)
  if _rc == 0 local vB = _b[main_att]
  cap qui wooldid y repiC t repons, att fastnose cluster(state) esrelativeto(-1)
  if _rc == 0 local vC = _b[main_att]
  post `memhold' (`r') (`saved`r'') (`vA') (`vB') (`vC')
  dis "REP `r' saved=`savedx`r'' A=" %21x `vA' " B=" %21x `vB' " C=" %21x `vC'
}
postclose `memhold'

use "probe02res.dta", clear
gen byte okA = (sav == extA)
gen byte okB = (sav == extB)
gen byte okC = (sav == extC)
gen double dA = abs(sav-extA)
gen double dB = abs(sav-extB)
gen double dC = abs(sav-extC)
sum okA okB okC dA dB dC
list, sep(0)
dis "SUMOK_A=" r(N)
count if okA == 0
dis "NMISMATCH_A=" r(N)
count if okB == 0
dis "NMISMATCH_B=" r(N)
count if okC == 0
dis "NMISMATCH_C=" r(N)
