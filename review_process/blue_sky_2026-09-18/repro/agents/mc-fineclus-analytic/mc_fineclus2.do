* ---------------------------------------------------------------------------------------
* mc_fineclus2.do — second wave of the finer-than-i analytic-SE coverage MC.
* args: 1 = kind (X | U)   2 = shard   3 = ndraws   4 = seedbase   5 = outdir
*
* DGP-X ("cross-cutting"): 40 states x 5 units x 12 periods, dependence at the STATE level
*   (s_g, v_gt) as in DGP-S, but ONSET IS ASSIGNED AT THE UNIT LEVEL: inside every state
*   unit 1 gets ttre=5, unit 2 gets ttre=8, units 3-5 are never treated. So every state
*   holds all three cohorts: `state' neither nests nor is nested in `coh'. This is the
*   CROSS-CUTTING geometry the task asks about.
*
* DGP-U ("ragged administrative"): onset at the STATE level as in DGP-S, but
*   - state sizes are UNEQUAL and fixed: M_g = 2 + mod(3*g, 9), so 2..10 units per state
*   - the panel is UNBALANCED on a fixed deterministic pattern (mod(unit*t,17)==0 dropped,
*     ~6% of rows), so the design - and therefore the truth - is still fixed across draws
*   This breaks the balanced-panel invariance, so estimand composition actually matters.
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

cap program drop mkd2
program define mkd2
  args kind nsd
  clear
  local G = 40
  local T = 12
  qui set obs `G'
  gen double state = _n
  gen double s_g   = rnormal()*0.8*`nsd'
  if "`kind'" == "X" {
    local M = 5
    qui expand `M'
    bysort state: gen double upos = _n
    gen double unit = (state-1)*`M' + upos
    * onset at the UNIT level: every state holds cohort 5, cohort 8 and never-treated
    gen double ttre = cond(upos==1, 5, cond(upos==2, 8, .))
  }
  else {
    * unequal, fixed state sizes 2..10
    gen double Mg = 2 + mod(3*state, 9)
    qui expand Mg
    bysort state: gen double upos = _n
    gen double unit = state*1000 + upos
    * onset at the STATE level
    gen double ttre = cond(state<=10, 5, cond(state<=20, 8, .))
  }
  gen double a_u = rnormal()*0.5*`nsd'
  qui expand `T'
  bysort unit: gen double t = _n
  tempvar first vv
  bysort state t: gen byte `first' = (_n==1)
  gen double `vv' = rnormal()*0.4*`nsd' if `first'
  bysort state t (`first'): replace `vv' = `vv'[_N]
  gen double v_gt = `vv'
  tempvar eps
  gen double `eps' = rnormal()*0.4*`nsd'
  bysort unit (t): gen double e = `eps'[1]
  by unit (t): replace e = 0.5*e[_n-1] + sqrt(1-0.25)*`eps'[_n] if _n > 1
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_gt + a_u + 0.2*t + 0.5*cos(t) + tau + e
  * ragged panel on a FIXED deterministic pattern (U only): the truth stays draw-invariant
  if "`kind'" == "U" qui drop if mod(upos*t, 17) == 0
  gen double g    = cond(ttre<., ttre, 0)
  gen double coh  = cond(ttre<., ttre, 99)
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

cap program drop runarm
program define runarm
  args H tag draw IV CV ESM
  local esopt esprelength(4) espostlength(4) esfixedbaseperiod
  if "`ESM'" == "pool" local esopt esprelength(0) espostlength(4)
  cap qui wooldid y `IV' t ttre, att cluster(`CV') `esopt'
  local rc = _rc
  if `rc' == 0 wrow `H' `tag' `draw'
  else file write `H' "err|`tag'|`draw'|`rc'" _n
  c_local armrc `rc'
end

if "`KIND'" == "X" {
  local ARMS `" "X3 coh state fb" "X3p coh state pool" "X4 coh unit fb" "X4p coh unit pool" "'
}
else {
  local ARMS `" "U2 state state fb" "U3 coh state fb" "U3p coh state pool" "U4 coh unit fb" "'
}

* ---- truth (near-zero-noise draw; design fixed across draws) ---------------------------
file open TR using "`OUTDIR'/truth_`KIND'.txt", write replace
set seed 777
set sortseed 777
mkd2 `KIND' 1e-4
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
use `td', clear
cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(state) long2
di "TRUTH_CS_RC=" _rc
if _rc == 0 {
  cap qui estat event, post
  if _rc == 0 wrow TR CS 0
}
use `td', clear
cap qui did_imputation y unit t ttre, horizons(0/3) pretrends(4) minn(0) cluster(state)
di "TRUTH_DI_RC=" _rc
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
  local sd = `BASE' + `SHARD'*1000000 + `d'*11 + cond("`KIND'"=="X", 3, 4)
  set seed `sd'
  set sortseed `sd'
  mkd2 `KIND' 1
  qui save `"`dd'"', replace
  foreach a of local ARMS {
    local tag : word 1 of `a'
    local IV  : word 2 of `a'
    local CV  : word 3 of `a'
    local ESM : word 4 of `a'
    use `dd', clear
    runarm MC `tag' `r' `IV' `CV' `ESM'
    if `armrc' == 0 {
      if "`tag'" == "X3" | "`tag'" == "U3" matrix BA = e(b)
      if "`tag'" == "X4" | "`tag'" == "U4" matrix BB = e(b)
    }
  }
  local md = .
  cap local md = mreldif(BA, BB)
  file write MC "inv|`r'|" %21.0g (`md') _n
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
  di "MCPROG kind=`KIND' shard=`SHARD' draw=`d' cum=" %8.1f r(t1)
  timer on 1
}
file write MC "MC_DONE `KIND' `SHARD'" _n
file close MC
di "MC_SHARD_DONE `KIND' `SHARD'"
exit, clear
