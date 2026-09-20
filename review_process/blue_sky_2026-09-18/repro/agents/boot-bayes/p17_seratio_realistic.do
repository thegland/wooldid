* p17: sd(bayes,equal)/sd(standard) on a REALISTIC weight distribution.
* p08b's weighted arm used p04_panel, whose weight has one cluster at 500x the
* others - a deliberately pathological design built for a different probe. Redo it
* with ordinary lognormal population weights, and with no weight.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
local R 600
foreach arm in none pop {
  local wt ""
  if "`arm'" == "pop" local wt "[aw=pop]"
  foreach tag in std eq siz {
    if "`tag'" == "std" {
      local md standard
      local alopt ""
    }
    if "`tag'" == "eq" {
      local md bayes
      local alopt bootalpha(equal)
    }
    if "`tag'" == "siz" {
      local md bayes
      local alopt bootalpha(size)
    }
    use p16_panel, clear
    cap noi wooldid y st year onset `wt', att cluster(st) bootstrap(`md') `alopt' reps(`R') rseed(5150) dots(0)
    matrix VV = e(V)
    local sd_`arm'_`tag' = sqrt(VV[1,1])
    dis "[`arm' `tag'] rc=" _rc " completed=" e(boot_completed) " failed=" e(boot_failed) " sd=" %12.8f sqrt(VV[1,1])
  }
}
use p16_panel, clear
qui count
local NN = r(N)
qui bysort st: gen byte f1 = _n==1
qui count if f1
local CC = r(N)
qui su pop, detail
dis "weight p10=" r(p10) " p50=" r(p50) " p90=" r(p90) "  max/min=" r(max)/r(min)
dis "N=`NN' C=`CC' sqrt(C/N)=" sqrt(`CC'/`NN')
foreach arm in none pop {
  dis "[`arm'] sd(bayes equal)/sd(standard)   = " %6.4f `sd_`arm'_eq'/`sd_`arm'_std'
  dis "[`arm'] sd(bayes size)/sd(bayes equal) = " %6.4f `sd_`arm'_siz'/`sd_`arm'_eq' " vs sqrt(C/N)=" %6.4f sqrt(`CC'/`NN')
}
dis "PROBE_DONE"
