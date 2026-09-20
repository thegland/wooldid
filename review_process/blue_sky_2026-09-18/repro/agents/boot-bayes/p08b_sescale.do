* p08b: two DOCUMENTED quantitative claims about the Dirichlet engine, plus the
* F001 centring re-verification.  (p08 died on an invalid macro name; same design.)
*  (1) sthlp:532 "shares [the pairs cluster bootstrap's] first two moments"
*      => sd(bayes,equal) ~= sd(standard) on the same data.
*  (2) sthlp:548 bootalpha(size) "shrinks the cluster-weight dispersion by roughly
*      the square root of the mean cluster size" => sd(size)/sd(equal) ~= sqrt(C/N).
*  (3) F001 (ado:14289-14312, fixed 2026-09-13): weighted unbalanced replications
*      must stay CENTRED. |mean(bootdist)-b| / sd(bootdist) was 9.36 before the fix.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

local R 400
foreach arm in unw wtd {
  local wt ""
  if "`arm'" == "wtd" local wt "[aw=bigpop]"
  foreach tag in std eq siz {
    if "`tag'" == "std" {
      local md standard
      local alopt ""
    }
    if "`tag'" == "eq"  {
      local md bayes
      local alopt bootalpha(equal)
    }
    if "`tag'" == "siz" {
      local md bayes
      local alopt bootalpha(size)
    }
    use p04_panel, clear
    cap noi wooldid y st year onset `wt', att cluster(st) bootstrap(`md') `alopt' reps(`R') rseed(808) dots(0)
    local rc = _rc
    matrix BB = e(b)
    matrix VV = e(V)
    matrix BD = e(bootdist)
    local b0 = BB[1,1]
    local se = sqrt(VV[1,1])
    preserve
      clear
      qui svmat double BD
      qui drop if missing(BD1)
      qui su BD1
      local nn = r(N)
      local md_ = r(mean)
      local sd_ = r(sd)
      restore
    dis "[`arm' `tag'] rc=`rc' n=`nn' b=" %12.8f `b0' " se(e(V))=" %12.8f `se' " sd(draws)=" %12.8f `sd_' " (mean-b)/sd = " %8.4f (`md_'-`b0')/`sd_'
    local sd_`arm'_`tag' = `sd_'
  }
}
use p04_panel, clear
qui count
local NN = r(N)
qui bysort st: gen byte f1 = _n==1
qui count if f1
local CC = r(N)
dis "N=`NN' C=`CC'  sqrt(C/N)=" sqrt(`CC'/`NN')
foreach arm in unw wtd {
  dis "[`arm'] sd(bayes equal)/sd(standard)   = " %6.4f `sd_`arm'_eq'/`sd_`arm'_std'  "   (documented claim ~1)"
  dis "[`arm'] sd(bayes size)/sd(bayes equal) = " %6.4f `sd_`arm'_siz'/`sd_`arm'_eq' "   (documented claim ~ sqrt(C/N) = " %6.4f sqrt(`CC'/`NN') ")"
}
dis "PROBE_DONE"
