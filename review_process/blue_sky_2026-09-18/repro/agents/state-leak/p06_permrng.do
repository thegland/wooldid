* p06: permute() and the session RNG -- seeding, reproducibility, and what is left behind
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

cap prog drop wdperm
prog def wdperm
  args tag
  dis "PRM|`tag'|method=" "`e(perm_method)'"
  dis "PRM|`tag'|seed=" "`e(perm_seed)'"
  dis "PRM|`tag'|rng_head=" substr("`e(perm_rngstate)'",1,24)
  dis "PRM|`tag'|rng_len=" strlen("`e(perm_rngstate)'")
  capture confirm matrix e(wdidmainresults)
  if _rc == 0 {
    tempname M
    matrix `M' = e(wdidmainresults)
    local rn : rownames `M'
    local cn : colnames `M'
    local nr = rowsof(`M')
    local nc = colsof(`M')
    forvalues a = 1/`nr' {
      local rr : word `a' of `rn'
      forvalues b = 1/`nc' {
        local cc2 : word `b' of `cn'
        dis "PRM|`tag'|`rr'|`cc2'=" %21x `M'[`a',`b']
      }
    }
  }
end

* ============ ARM 0: plain reference (no permute) -- the observed estimates
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "R0_RC=" _rc
wdperm R0_plain

* ============ ARM A: permute, no rseed, user set seed 777 beforehand
wdgen
set seed 777
global WDH_R1 "`c(rngstate)'"
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(30) dots(0)
dis "A_RC=" _rc
global WDH_R2 "`c(rngstate)'"
dis "A_RNG_ADVANCED=" ("$WDH_R1" != "$WDH_R2")
global WDH_PERMRNG "`e(perm_rngstate)'"
dis "A_PERMRNG_EQ_PRECALL_STATE=" ("$WDH_PERMRNG" == "$WDH_R1")
wdperm A_seed777
wdline A_post

* ============ ARM B: identical call after set seed 777 again -- must reproduce exactly
wdgen
set seed 777
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(30) dots(0)
dis "B_RC=" _rc
wdperm B_seed777_again

* ============ ARM C: replay by e(perm_rngstate), as the help promises
wdgen
capture noisily set rngstate "$WDH_PERMRNG"
dis "C_SETRNGSTATE_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(30) dots(0)
dis "C_RC=" _rc
wdperm C_rngstate_replay

* ============ ARM D: rseed(), twice, from DIFFERENT starting streams
wdgen
set seed 111
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(30) dots(0) rseed(4242)
dis "D_RC=" _rc
wdperm D_rseed4242_from111
wdgen
set seed 999
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(30) dots(0) rseed(4242)
dis "D2_RC=" _rc
wdperm D2_rseed4242_from999

* ============ ARM E: does rseed() leak the seeded stream into the user's session?
wdgen
set seed 555
gen double __u0 = .
qui replace __u0 = runiform() in 1
scalar __u0s = __u0[1]
set seed 555
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(20) dots(0) rseed(4242)
dis "E_RC=" _rc
qui gen double __u1 = .
qui replace __u1 = runiform() in 1
dis "E_U_WITHOUT_WOOLDID=" %21x __u0s
dis "E_U_AFTER_WOOLDID=" %21x __u1[1]
dis "E_STREAM_DISTURBED=" (__u0s != __u1[1])
qui drop __u0 __u1

* ============ ARM F: permsave file + its characteristics
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(20) dots(0) rseed(4242) permsave(p06perm) replace
dis "F_RC=" _rc
wdline F_post
capture noisily use "p06perm.dta", clear
dis "F_USE_RC=" _rc
local chs : char _dta[]
dis "F_CHARS=" "`chs'"
foreach ch of local chs {
  local cv : char _dta[`ch']
  dis "F_CHAR|`ch'|len=" strlen("`cv'") "|head=" substr("`cv'",1,20)
}

dis "P06_DONE"
