* ---------------------------------------------------------------------------------------
* mc_fineclus.do — does the documented remedy ("choose i COARSE, cluster FINER") deliver
* correct analytic-SE coverage?  Monte Carlo.  2026-09-17.
*
* args: 1 = kind (S | N)   2 = shard   3 = ndraws   4 = seedbase   5 = outdir
*
* DGP-S ("state policy"): 40 states x 5 units x 12 periods = 2400 rows.
*   onset at the STATE level: states 1-10 ttre=5, 11-20 ttre=8, 21-40 never.
*   y = s_g + v_gt + a_u + .2 t + .5 cos t + tau + e
*     s_g  ~ N(0,.8^2)  state random effect
*     v_gt ~ N(0,.4^2)  state-by-time shock  <- the component a saturated state-by-relperiod
*                                               dummy ABSORBS, which is the whole mechanism
*     a_u  ~ N(0,.5^2)  unit random effect
*     e    AR(1) rho .5, innovation sd .4, within unit
*   tau = 1.0 + .05*rel (ttre=5 states), 0.4 + .05*rel (ttre=8 states)
*   TRUE cluster level = state.
* DGP-N: the same code with s_g = v_gt = 0, so the ONLY dependence is within unit and
*   cluster(unit) is the correct level.
*
* Every arm is run on the SAME draw. Truth is a zero-noise (1e-6) draw, computed once per
* kind per arm at the top, because the treatment design is fixed across draws.
* ---------------------------------------------------------------------------------------
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local KIND   `1'
local SHARD  = `2'
local NDRAWS = `3'
local BASE   = `4'
local OUTDIR `5'

cap program drop mkd
program define mkd
  * args: kind (S|N)  nsd (noise scale)
  args kind nsd
  clear
  local G = 40
  local M = 5
  local T = 12
  local SS = 0.8
  local SV = 0.4
  if "`kind'" == "N" {
    local SS = 0
    local SV = 0
  }
  qui set obs `G'
  gen double state = _n
  gen double s_g   = rnormal()*`SS'*`nsd'
  gen double ttre  = cond(state<=10, 5, cond(state<=20, 8, .))
  qui expand `M'
  bysort state: gen double unit = (state-1)*`M' + _n
  gen double a_u   = rnormal()*0.5*`nsd'
  qui expand `T'
  bysort unit: gen double t = _n
  * state-by-time shock: one draw per (state,t)
  tempvar first vv
  bysort state t: gen byte `first' = (_n==1)
  gen double `vv' = rnormal()*`SV'*`nsd' if `first'
  bysort state t (`first'): replace `vv' = `vv'[_N]
  gen double v_gt = `vv'
  * AR(1) idiosyncratic error within unit
  tempvar eps
  gen double `eps' = rnormal()*0.4*`nsd'
  bysort unit (t): gen double e = `eps'[1]
  by unit (t): replace e = 0.5*e[_n-1] + sqrt(1-0.25)*`eps'[_n] if _n > 1
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_gt + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double g    = cond(ttre<., ttre, 0)
  gen double coh  = cond(ttre<., ttre, 99)
  gen double icoh = cond(ttre<., ttre, 1000+unit)
  keep y state unit t ttre rel g coh icoh
  qui compress
end

* ---- writers ---------------------------------------------------------------------------
cap program drop wrow
program define wrow
  * args: handle tag draw  (uses e(b), e(V), e(N_clust), e(df_r))
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

* ---- arm runner ------------------------------------------------------------------------
cap program drop runarm
program define runarm
  * args: handle tag draw ivar clusvar esmode
  args H tag draw IV CV ESM
  local esopt esprelength(4) espostlength(4) esfixedbaseperiod
  if "`ESM'" == "pool" local esopt esprelength(0) espostlength(4)
  cap qui wooldid y `IV' t ttre, att cluster(`CV') `esopt'
  local rc = _rc
  if `rc' == 0 wrow `H' `tag' `draw'
  else file write `H' "err|`tag'|`draw'|`rc'" _n
  c_local armrc `rc'
end

* ---- arm table -------------------------------------------------------------------------
* tag  ivar  clustervar  esmode
if "`KIND'" == "S" {
  local ARMS `" "S2 state state fb" "S2p state state pool" "S3 coh state fb" "S3p coh state pool" "S4 coh unit fb" "S4p coh unit pool" "'
  local EXTCL state
}
else {
  local ARMS `" "N3 coh unit fb" "N3p coh unit pool" "N2 state unit fb" "NC coh state fb" "NCp coh state pool" "NS state state fb" "'
  local EXTCL unit
}

* ---- truth (zero-noise draw, design is fixed across draws) -----------------------------
file open TR using "`OUTDIR'/truth_`KIND'.txt", write replace
set seed 777
set sortseed 777
mkd `KIND' 1e-6
tempfile td
qui save `td'
foreach a of local ARMS {
  local tag : word 1 of `a'
  local IV  : word 2 of `a'
  local CV  : word 3 of `a'
  local ESM : word 4 of `a'
  use `td', clear
  runarm TR `tag' 0 `IV' `CV' `ESM'
  di "TRUTH_`tag'_RC=`armrc'"
}
* external oracles' truth too
use `td', clear
cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(`EXTCL') long2
if _rc == 0 {
  cap qui estat event, post
  if _rc == 0 wrow TR CS 0
}
use `td', clear
cap qui did_imputation y unit t ttre, horizons(0/3) pretrends(4) minn(0) cluster(`EXTCL')
if _rc == 0 wrow TR DI 0
file close TR

* ---- main loop -------------------------------------------------------------------------
file open MC using "`OUTDIR'/mc_`KIND'_s`SHARD'.txt", write replace
timer clear 1
timer on 1
tempfile dd
forval d = 1/`NDRAWS' {
  cap matrix drop BA
  cap matrix drop BB
  local r = (`SHARD'-1)*100000 + `d'
  local sd = `BASE' + `SHARD'*1000000 + `d'*7 + cond("`KIND'"=="S", 1, 2)
  set seed `sd'
  set sortseed `sd'
  mkd `KIND' 1
  qui save `"`dd'"', replace
  * wooldid arms
  local bkeep
  foreach a of local ARMS {
    local tag : word 1 of `a'
    local IV  : word 2 of `a'
    local CV  : word 3 of `a'
    local ESM : word 4 of `a'
    use `dd', clear
    runarm MC `tag' `r' `IV' `CV' `ESM'
    if `armrc' == 0 {
      if "`tag'" == "S3" | "`tag'" == "N3" matrix BA = e(b)
      if "`tag'" == "S4" | "`tag'" == "NC" matrix BB = e(b)
    }
  }
  * point-estimate invariance to cluster(): same i, different cluster()
  local md = .
  cap local md = mreldif(BA, BB)
  file write MC "inv|`r'|" %21.0g (`md') _n
  * external controls
  use `dd', clear
  cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(`EXTCL') long2
  if _rc == 0 {
    cap qui estat event, post
    if _rc == 0 wrow MC CS `r'
    else file write MC "err|CSEV|`r'|" (_rc) _n
  }
  else file write MC "err|CS|`r'|" (_rc) _n
  use `dd', clear
  cap qui did_imputation y unit t ttre, horizons(0/3) pretrends(4) minn(0) cluster(`EXTCL')
  if _rc == 0 wrow MC DI `r'
  else file write MC "err|DI|`r'|" (_rc) _n

  timer off 1
  qui timer list 1
  file write MC "progress|`r'|" %9.1f (r(t1)) _n
  di "MCPROG kind=`KIND' shard=`SHARD' draw=`d' cum=" %8.1f r(t1)
  timer on 1
}
file write MC "MC_DONE `KIND' `SHARD'" _n
file close MC
di "MC_SHARD_DONE `KIND' `SHARD'"
exit, clear
