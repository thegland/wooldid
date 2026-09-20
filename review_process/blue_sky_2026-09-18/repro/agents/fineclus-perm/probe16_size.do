* probe16: SIZE of the permutation test under the rebuild, on a clean EXCHANGEABLE design
* 9 clusters, balanced, i = 3 onset cohorts, cluster(state) strictly finer than i,
* sharp null (no treatment effect at all). Rejection rate must be <= nominal.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

local NSIM 80
tempname pf
postfile `pf' int sim double p double b using "probe16res.dta", replace
forval s = 1/`NSIM' {
  clear
  set seed `=20000 + `s''
  set sortseed `=20000 + `s''
  qui set obs 9
  qui gen int state = _n
  qui gen double sfx = rnormal()*1.5
  qui expand 10
  qui bysort state: gen int t = _n
  qui gen double eps = rnormal()
  qui gen int coh = cond(state<=3,1,cond(state<=6,2,3))
  qui gen double onset = cond(coh==1,4,cond(coh==2,7,.))
  * SHARP NULL: no treatment effect anywhere
  qui gen double y = 1 + sfx + 0.25*t + eps
  local pp = .
  local bb = .
  cap qui wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(99) rseed(`=700+`s'') permonsetcohorts dots(0)
  if _rc == 0 {
    local pp = e(perm_p_att)
    local bb = _b[main_att]
  }
  post `pf' (`s') (`pp') (`bb')
}
postclose `pf'
use "probe16res.dta", clear
count
dis "SIZE_NSIM=" r(N)
count if missing(p)
dis "SIZE_NMISSING=" r(N)
gen byte r10 = (p <= 0.10) if !missing(p)
gen byte r05 = (p <= 0.05) if !missing(p)
gen byte r20 = (p <= 0.20) if !missing(p)
sum r05 r10 r20 p
dis "SIZE_R05=" 0
qui sum r05
dis "SIZE_R05=" r(mean) " n=" r(N)
qui sum r10
dis "SIZE_R10=" r(mean) " n=" r(N)
qui sum r20
dis "SIZE_R20=" r(mean) " n=" r(N)
