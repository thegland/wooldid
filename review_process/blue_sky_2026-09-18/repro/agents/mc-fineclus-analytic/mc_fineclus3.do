* ---------------------------------------------------------------------------------------
* mc_fineclus3.do — THE HEAD-TO-HEAD: the documented defect and the documented remedy on
* the SAME draws, on the exact geometry where the defect was documented.
* args: 1 = shard   2 = ndraws   3 = seedbase   4 = outdir
*
* DGP-P: 40 independent units, T = 12 (the prior 500-draw MC's shape, scaled the same way).
*   units 1-10 ttre=5, 11-20 ttre=8, 21-40 never treated.
*   y = a_u + .2 t + .5 cos t + tau + e ,  a_u ~ N(0,1), e AR(1) rho .5, innovation sd .4
*   tau = 1.0 + .05*rel (ttre=5), 0.4 + .05*rel (ttre=8)
*   Units are independent, so cluster(unit) is EXACTLY the right level: any undercoverage
*   here is the estimator's, not the DGP's.
*
* Arms, all on the same draw, 40 clusters in every one:
*   P1  i=unit  cluster(unit)  esfixedbaseperiod  <- the DOCUMENTED-BAD arm (cluster == i)
*   P1p i=unit  cluster(unit)  esprelength(0)     <- the prior MC's WORST shape (0.661/0.794)
*   P2  i=coh   cluster(unit)  esfixedbaseperiod  <- ** THE REMEDY: i coarse, cluster finer **
*   P2p i=coh   cluster(unit)  esprelength(0)     <- ** THE REMEDY, pooled shape **
*   P3  i=icoh  cluster(icoh)  esfixedbaseperiod  <- the prior MC's "clean" arm, for continuity
*   CS / DI with cluster(unit) as external controls
* ---------------------------------------------------------------------------------------
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local SHARD  = `1'
local NDRAWS = `2'
local BASE   = `3'
local OUTDIR `4'

cap program drop mkp
program define mkp
  args nsd
  clear
  local G = 40
  local T = 12
  qui set obs `G'
  gen double unit = _n
  gen double a_u  = rnormal()*1.0*`nsd'
  gen double ttre = cond(unit<=10, 5, cond(unit<=20, 8, .))
  qui expand `T'
  bysort unit: gen double t = _n
  tempvar eps
  gen double `eps' = rnormal()*0.4*`nsd'
  bysort unit (t): gen double e = `eps'[1]
  by unit (t): replace e = 0.5*e[_n-1] + sqrt(1-0.25)*`eps'[_n] if _n > 1
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
  gen double y = a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double g    = cond(ttre<., ttre, 0)
  gen double coh  = cond(ttre<., ttre, 99)
  gen double icoh = cond(ttre<., ttre, 1000+unit)
  keep y unit t ttre rel g coh icoh
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

local ARMS `" "P1 unit unit fb" "P1p unit unit pool" "P2 coh unit fb" "P2p coh unit pool" "P3 icoh icoh fb" "'

file open TR using "`OUTDIR'/truth_P.txt", write replace
set seed 777
set sortseed 777
mkp 1e-4
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
cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(unit) long2
di "TRUTH_CS_RC=" _rc
if _rc == 0 {
  cap qui estat event, post
  if _rc == 0 wrow TR CS 0
}
use `td', clear
cap qui did_imputation y unit t ttre, horizons(0/3) pretrends(4) minn(0) cluster(unit)
di "TRUTH_DI_RC=" _rc
if _rc == 0 wrow TR DI 0
file close TR

file open MC using "`OUTDIR'/mc_P_s`SHARD'.txt", write replace
timer clear 1
timer on 1
tempfile dd
forval d = 1/`NDRAWS' {
  cap matrix drop BA
  cap matrix drop BB
  local r = (`SHARD'-1)*100000 + `d'
  local sd = `BASE' + `SHARD'*1000000 + `d'*13 + 5
  set seed `sd'
  set sortseed `sd'
  mkp 1
  qui save `"`dd'"', replace
  foreach a of local ARMS {
    local tag : word 1 of `a'
    local IV  : word 2 of `a'
    local CV  : word 3 of `a'
    local ESM : word 4 of `a'
    use `dd', clear
    runarm MC `tag' `r' `IV' `CV' `ESM'
    if `armrc' == 0 {
      if "`tag'" == "P1" matrix BA = e(b)
      if "`tag'" == "P2" matrix BB = e(b)
    }
  }
  * b invariance ACROSS i (unit vs coh) on a balanced panel
  local md = .
  cap local md = mreldif(BA, BB)
  file write MC "inv|`r'|" %21.0g (`md') _n
  use `dd', clear
  cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(unit) long2
  if _rc == 0 {
    cap qui estat event, post
    if _rc == 0 wrow MC CS `r'
    else file write MC "err|CSEV|`r'|" (_rc) _n
  }
  else file write MC "err|CS|`r'|" (_rc) _n
  use `dd', clear
  cap qui did_imputation y unit t ttre, horizons(0/3) pretrends(4) minn(0) cluster(unit)
  if _rc == 0 wrow MC DI `r'
  else file write MC "err|DI|`r'|" (_rc) _n
  timer off 1
  qui timer list 1
  file write MC "progress|`r'|" %9.1f (r(t1)) _n
  di "MCPROG kind=P shard=`SHARD' draw=`d' cum=" %8.1f r(t1)
  timer on 1
}
file write MC "MC_DONE P `SHARD'" _n
file close MC
di "MC_SHARD_DONE P `SHARD'"
exit, clear
