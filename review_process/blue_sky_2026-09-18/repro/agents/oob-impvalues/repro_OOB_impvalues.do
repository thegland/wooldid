*=============================================================================
* SELF-CONTAINED REPRO - oob-impvalues
*
* CLAIM: under `impvalues', when cluster() is FINER than i or CROSSES i, the G
* "per-cluster Ibragimov-Muller estimates" are not G estimates. They are G
* weighted re-reads of ONE pooled fit, so refining the cluster label multiplies
* the reported standard error by a factor that depends only on the LABEL COUNT.
* The 95% interval is then far too narrow and no caution fires.
*
* Three blocks, ~3 minutes:
*   BLOCK 1 (instant) the algebraic identity: b is bit-identical across two
*                     cluster granularities and se scales by sqrt((K-1)/(nK-1)).
*   BLOCK 2 (instant) wooldid itself REFUSES (r459) to estimate an effect from
*                     one unit's data - yet impvalues reports 60 of them.
*   BLOCK 3 (~3 min)  R=40 coverage check at the pairing the help recommends
*                     (i = onset cohort, cluster below i), with the analytic
*                     clustered arm on the SAME draws as the positive control.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_OOB_impvalues.do
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkD
prog def mkD
  args sd
  clear
  set obs 70
  gen long unit = _n
  gen byte  cg  = floor((_n-1)/5)+1          /* 12 treated i-levels + 2 never-treated */
  gen int   reg = mod(_n-1,14)+1             /* regions CROSS the cohorts             */
  gen double ai = 0.8*mod(_n*23,31)/31 - 0.4
  expand 10
  bys unit: gen int t = 2000 + _n
  gen double ttre = .
  replace ttre = 2004 if inrange(cg,1,3)
  replace ttre = 2005 if inrange(cg,4,6)
  replace ttre = 2006 if inrange(cg,7,9)
  replace ttre = 2007 if inrange(cg,10,12)
  gen double truetau = 0
  replace truetau = 0.30 + 0.12*(cg-1) + 0.04*(mod(cg,3)-1)*(t-ttre) if cg<=12 & t>=ttre
  gen double base = ai + 0.35*(t-2000) - 0.012*(t-2000)^2 + truetau
  set seed `sd'
  set sortseed `sd'
  by unit: gen double us = rnormal(0,0.15) if _n==1    /* unit random TREND */
  by unit: replace us = us[1]
  gen double y = base + us*(t-2005.5) + rnormal(0,0.60)
end

*========================= BLOCK 1 : the identity ============================
local K 12
local n 5
di as txt _n "PREDICTED se(cluster=unit)/se(cluster=cg) = sqrt((K-1)/(nK-1)) = " %14.12f sqrt((`K'-1)/(`n'*`K'-1))
foreach s in 4101 4102 4103 {
  mkD `s'
  qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
  mat Q = e(wdidmainresults)
  local bu = Q[1,1]
  local su = Q[1,2]
  local du = e(df_r)
  qui wooldid y cg t ttre, att cluster(cg) impvalues level(95)
  mat Q = e(wdidmainresults)
  local bc = Q[1,1]
  local sc = Q[1,2]
  local dc = e(df_r)
  di as txt "seed `s': b(unit)=" %20.14f `bu' "  b(cg)=" %20.14f `bc' "  |diff|=" %10.3e abs(`bu'-`bc')
  di as txt "          se(unit)=" %12.9f `su' "  se(cg)=" %12.9f `sc' "  ratio=" %14.12f `su'/`sc' "  df posted: " `du' " vs " `dc'
}
di as txt "-> the 60 per-unit 'estimates' are 5 identical copies each of 12 pooled reads."
di as txt "   Refining the label from 12 cohorts to 60 units changes NO data, shrinks the"
di as txt "   reported se by 0.4318, and raises the posted reference df from 11 to 59."

*========================= BLOCK 2 : no such estimate exists ==================
mkD 4101
tempfile full
qui save `full'
di as txt _n "---- can ONE unit support its own DiD estimate? ----"
use `full', clear
keep if unit==3
cap noi wooldid y cg t ttre, att cluster(unit) level(95)
di as txt "rc = " _rc "  (459 = wooldid cannot estimate an effect from one unit's data)"
di as txt "Yet cluster(unit) + impvalues reports 60 such 'per-cluster estimates' and"
di as txt "takes a t test on 59 degrees of freedom across them."

*========================= BLOCK 3 : coverage =================================
mkD 1
replace y = base
qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
mat Q = e(wdidmainresults)
local TIM = Q[1,1]
qui wooldid y cg t ttre, att cluster(unit) level(95)
mat Q = e(wdidmainresults)
local TAN = Q[1,1]
di as txt _n "IM target (noiseless cluster-equally-weighted mean) = " %12.8f `TIM'
di as txt "analytic target (noiseless cell-weighted att)       = " %12.8f `TAN'
local R 40
local cIM 0
local cAN 0
local sIM 0
local sAN 0
local sb 0
local sbb 0
forvalues r = 1/`R' {
  mkD `=880000+`r''
  qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
  mat Q = e(wdidmainresults)
  local sb  = `sb'  + Q[1,1]
  local sbb = `sbb' + Q[1,1]^2
  local sIM = `sIM' + Q[1,2]
  if Q[1,5] <= `TIM' & `TIM' <= Q[1,6] local ++cIM
  qui wooldid y cg t ttre, att cluster(unit) level(95)
  mat Q = e(wdidmainresults)
  local sAN = `sAN' + Q[1,2]
  if Q[1,5] <= `TAN' & `TAN' <= Q[1,6] local ++cAN
}
local sdb = sqrt((`sbb' - `R'*(`sb'/`R')^2)/(`R'-1))
di as txt _n "================= R=`R', i = onset cohort, cluster(unit) ================="
di as txt "true sampling sd of the point estimate        : " %10.6f `sdb'
di as txt "mean reported impvalues (IM) standard error   : " %10.6f (`sIM'/`R')
di as txt "mean reported analytic clustered se           : " %10.6f (`sAN'/`R')
di as txt "impvalues 95% CI coverage                     : " %6.3f (`cIM'/`R')
di as txt "analytic  95% CI coverage (positive control)  : " %6.3f (`cAN'/`R')
di as txt "=========================================================================="
di as txt "REPRO_DONE"
