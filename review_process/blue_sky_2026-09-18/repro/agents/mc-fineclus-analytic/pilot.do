clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
set graphics off
which wooldid

cap program drop mkS
program define mkS
  * args: nsd = noise scale (1 = full, 1e-6 = "truth" draw)
  args nsd
  clear
  local G  = 40
  local M  = 5
  local T  = 12
  qui set obs `G'
  gen double state = _n
  gen double s_g   = rnormal()*0.8*`nsd'
  gen double ttre  = cond(state<=10, 5, cond(state<=20, 8, .))
  qui expand `M'
  bysort state: gen double unit = (state-1)*`M' + _n
  gen double a_u   = rnormal()*0.5*`nsd'
  qui expand `T'
  bysort unit: gen double t = _n
  gen double e     = rnormal()*0.4*`nsd'
  * state-by-time shock: one draw per (state,t)
  tempvar one v
  bysort state t: gen byte `one' = (_n==1)
  gen double `v' = rnormal()*0.4*`nsd' if `one'
  bysort state t (`one'): replace `v' = `v'[_N]
  gen double v_gt = `v'
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_gt + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double g   = cond(ttre<., ttre, 0)
  gen double coh = cond(ttre<., ttre, 99)
  gen double icoh= cond(ttre<., ttre, 1000+unit)
  gen double era = ceil(t/4)
  egen double uera = group(unit era)
  qui compress
end

mkS 1
qui count
dis "NROWS=" r(N)
tab ttre, missing
qui levelsof coh, local(CC)
dis "COHVALS=`CC'"

timer clear
foreach arm in S1 S2 S3 S4 {
  local IV
  local CV
  if "`arm'"=="S1" {
    local IV unit
    local CV state
  }
  if "`arm'"=="S2" {
    local IV state
    local CV state
  }
  if "`arm'"=="S3" {
    local IV coh
    local CV state
  }
  if "`arm'"=="S4" {
    local IV coh
    local CV unit
  }
  timer clear 1
  timer on 1
  cap noisily wooldid y `IV' t ttre, att cluster(`CV') esprelength(4) espostlength(4) esfixedbaseperiod
  local rc = _rc
  dis "ARM_`arm'_RC=`rc'"
  timer off 1
  qui timer list 1
  dis "ARM_`arm'_TIME=" r(t1)
  if `rc'==0 {
    dis "ARM_`arm'_NCLUST=" e(N_clust)
    matrix list e(b)
    matrix list e(V)
  }
}

timer clear 2
timer on 2
cap noisily csdid y, ivar(unit) time(t) gvar(g) cluster(state) long2
dis "CS_RC=" _rc
if _rc==0 {
  cap noisily estat event, post
  dis "CSEV_RC=" _rc
  if _rc==0 matrix list e(b)
}
timer off 2
qui timer list 2
dis "CS_TIME=" r(t2)

timer clear 3
timer on 3
cap noisily did_imputation y unit t ttre, horizons(0/4) pretrends(4) minn(0) cluster(state)
dis "DI_RC=" _rc
if _rc==0 matrix list e(b)
timer off 3
qui timer list 3
dis "DI_TIME=" r(t3)
