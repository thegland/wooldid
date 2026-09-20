* ---------------------------------------------------------------------------------------
* mc_wave3.do — THE LADDER. Everything held fixed except the number of CLUSTERS INSIDE A
* TREATED COHORT-BY-PERIOD CELL, which is the condition wooldid.sthlp:251 names for the
* finer-than-i remedy ("so that a treated cohort-by-period cell contains several clusters").
*
* args: 1 = K (states per treated cohort: 3, 5 or 10)  2 = shard  3 = ndraws  4 = seedbase  5 = outdir
*
* 40 states x 4 units x 12 periods in every rung. Two adoption waves, t=5 and t=8, with K
* states each; the remaining 40-2K states never adopt. Dependence at the STATE level
* (state RE .8, state-by-year shock .4) plus unit RE .5 and AR(1) rho .5 errors, so
* cluster(state) is the correct level in every rung.
* Arms: R3  i=coh cluster(state) esfixedbaseperiod  <- the remedy
*       R3p i=coh cluster(state) esprelength(0)     <- the remedy, pooled
*       CS  csdid long2 cluster(state)              <- external control (same points)
*       DI  did_imputation cluster(state)           <- external control
* ---------------------------------------------------------------------------------------
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local K      = `1'
local SHARD  = `2'
local NDRAWS = `3'
local BASE   = `4'
local OUTDIR `5'

cap program drop mkk
program define mkk
  args K nsd
  clear
  qui set obs 40
  gen double state = _n
  gen double s_g   = rnormal()*0.8*`nsd'
  gen double ttre  = cond(state<=`K', 5, cond(state<=2*`K', 8, .))
  qui expand 4
  bysort state: gen double upos = _n
  gen double unit = (state-1)*4 + upos
  gen double a_u  = rnormal()*0.5*`nsd'
  qui expand 12
  bysort unit: gen double t = _n
  tempvar first vv eps
  bysort state t: gen byte `first' = (_n==1)
  gen double `vv' = rnormal()*0.4*`nsd' if `first'
  bysort state t (`first'): replace `vv' = `vv'[_N]
  gen double v_gt = `vv'
  gen double `eps' = rnormal()*0.4*`nsd'
  bysort unit (t): gen double e = `eps'[1]
  by unit (t): replace e = 0.5*e[_n-1] + sqrt(1-0.25)*`eps'[_n] if _n > 1
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_gt + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double g   = cond(ttre<., ttre, 0)
  gen double coh = cond(ttre<., ttre, 99)
  keep y state unit t ttre rel g coh
  qui compress
end

cap program drop wrow
program define wrow
  args H tag draw
  tempname B V
  matrix `B' = e(b)
  matrix `V' = e(V)
  local nc = e(N_clust)
  if "`nc'" == "" local nc .
  local dr = e(df_r)
  if "`dr'" == "" local dr .
  local names : colnames `B'
  local k = colsof(`B')
  forval j = 1/`k' {
    local nm : word `j' of `names'
    file write `H' "`tag'|`draw'|`nm'|" %21.0g (`B'[1,`j']) "|" %21.0g (sqrt(`V'[`j',`j'])) "|`nc'|`dr'" _n
  }
end

file open TR using "`OUTDIR'/truth_K`K'.txt", write replace
set seed 777
set sortseed 777
mkk `K' 1e-4
tempfile td
qui save `td'
use `td', clear
cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod
di "TRUTH_R3_RC=" _rc
if _rc == 0 wrow TR R3 0
use `td', clear
cap qui wooldid y coh t ttre, att cluster(state) esprelength(0) espostlength(4)
di "TRUTH_R3p_RC=" _rc
if _rc == 0 wrow TR R3p 0
file close TR

file open MC using "`OUTDIR'/mc_K`K'_s`SHARD'.txt", write replace
timer clear 1
timer on 1
tempfile dd
forval d = 1/`NDRAWS' {
  local r = (`SHARD'-1)*100000 + `d'
  local sd = `BASE' + `SHARD'*1000000 + `d'*19 + `K'
  set seed `sd'
  set sortseed `sd'
  mkk `K' 1
  qui save `"`dd'"', replace
  cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod
  if _rc == 0 wrow MC R3 `r'
  else file write MC "err|R3|`r'|" (_rc) _n
  use `dd', clear
  cap qui wooldid y coh t ttre, att cluster(state) esprelength(0) espostlength(4)
  if _rc == 0 wrow MC R3p `r'
  else file write MC "err|R3p|`r'|" (_rc) _n
  use `dd', clear
  cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(state) long2
  if _rc == 0 {
    cap qui estat event, post
    if _rc == 0 wrow MC CS `r'
    else file write MC "err|CSEV|`r'|" (_rc) _n
  }
  else file write MC "err|CS|`r'|" (_rc) _n
  use `dd', clear
  cap qui did_imputation y unit t ttre, horizons(0/3) pretrends(4) minn(0) cluster(state)
  if _rc == 0 wrow MC DI `r'
  else file write MC "err|DI|`r'|" (_rc) _n
  timer off 1
  qui timer list 1
  file write MC "progress|`r'|" %9.1f (r(t1)) _n
  di "MCPROG K=`K' shard=`SHARD' draw=`d' cum=" %8.1f r(t1)
  timer on 1
}
file write MC "MC_DONE K`K' `SHARD'" _n
file close MC
di "MC_SHARD_DONE K`K' `SHARD'"
exit, clear
