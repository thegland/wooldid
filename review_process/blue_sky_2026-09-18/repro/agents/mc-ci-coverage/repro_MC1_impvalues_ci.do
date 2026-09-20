*=========================================================================
* SELF-CONTAINED REPRO  MC1
*   The impvalues (Ibragimov-Muller) standard error and confidence interval
*   are wrong whenever SEVERAL CLUSTERS SHARE THE SAME LEVELS OF i.
*   Two blocks, ~2 minutes total:
*     BLOCK 1 (instant)  mechanism: on NOISY data the G per-cluster reads are
*                        numerically IDENTICAL, which is impossible for G
*                        separate fits and necessary for G re-weightings of one.
*     BLOCK 2 (~90 sec)  consequence: 40-replication Monte Carlo at the pairing
*                        the help RECOMMENDS (i = cohort, cluster = unit).
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_MC1_impvalues_ci.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*=========================== BLOCK 1 =====================================
* 100 units, i = cg (3 treated cohorts + never-treated), cluster(reg) = 25
* regions that CUT ACROSS the cohorts and all have the SAME composition
* (exactly one unit from each cohort).  y carries an idiosyncratic sd of 0.60
* AND a unit shock of sd 0.50 -- this data is not remotely deterministic.
clear
set obs 100
gen long unit = _n
gen byte cg  = mod(_n-1,4)+1
gen int  reg = mod(_n-1,25)+1
gen double ai = 0.9*mod(_n*13,17)/17 - 0.45
expand 8
bys unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 if cg==1
replace ttre = 2004 if cg==2
replace ttre = 2005 if cg==3
gen double truetau = 0
replace truetau = 1.0 + 0.30*(t-ttre) if cg==1 & t>=ttre
replace truetau = 0.6 + 0.15*(t-ttre) if cg==2 & t>=ttre
replace truetau = 1.4 - 0.10*(t-ttre) if cg==3 & t>=ttre
set seed 31
set sortseed 31
by unit: gen double ui = rnormal(0,0.50) if _n==1
by unit: replace ui = ui[1]
gen double y = ai + 0.4*(t-1999) + truetau + ui + rnormal(0,0.60)
qui su y
di as txt "sd(y) = " as res r(sd) as txt "   (the data is NOISY)"
di as txt _n "---- impvalues here blanks EVERY se/t/p/CI and blames the data ----"
noi wooldid y cg t ttre, att cluster(reg) impvalues level(95)
di as txt _n "---- the analytic clustered run on the SAME data has no trouble ----"
noi wooldid y cg t ttre, att cluster(reg) level(95)

*=========================== BLOCK 2 =====================================
* i = cohort, cluster(unit): the cluster is FINER than i, which is the pairing
* the help recommends under "What level to cluster at".  Effects 1.00/1.05/0.95.
cap prog drop mkQ
prog def mkQ
  args sd
  clear
  set obs 40
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ai = 0.9*mod(_n*13,17)/17 - 0.45
  expand 8
  bys unit: gen int t = 1999 + _n
  gen double ttre = .
  replace ttre = 2002 if cg==1
  replace ttre = 2004 if cg==2
  replace ttre = 2005 if cg==3
  gen double truetau = 0
  replace truetau = 1.00 if cg==1 & t>=ttre
  replace truetau = 1.05 if cg==2 & t>=ttre
  replace truetau = 0.95 if cg==3 & t>=ttre
  gen double base = ai + 0.4*(t-1999) + truetau
  set seed `sd'
  set sortseed `sd'
  by unit: gen double ui = rnormal(0,0.40) if _n==1
  by unit: replace ui = ui[1]
  gen double y = base + ui + rnormal(0,0.60)
end
mkQ 1
replace y = base
qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
local TI = e(b)[1,1]
di as txt _n "TARGET (noiseless impvalues point estimate) = " as res %12.9f `TI'
local R 40
local nc_im 0
local nc_an 0
local s_im 0
local s_an 0
local sb 0
local sbb 0
forvalues r = 1/`R' {
  mkQ `=1500000+`r''
  qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
  mat Q = e(wdidmainresults)
  local sb  = `sb'  + Q[1,1]
  local sbb = `sbb' + Q[1,1]^2
  local s_im = `s_im' + Q[1,2]
  if Q[1,5] <= `TI' & `TI' <= Q[1,6] local ++nc_im
  qui wooldid y cg t ttre, att cluster(unit) level(95)
  mat Q = e(wdidmainresults)
  local s_an = `s_an' + Q[1,2]
  if Q[1,5] <= `TI' & `TI' <= Q[1,6] local ++nc_an
}
local sdb = sqrt(`sbb'/`R' - (`sb'/`R')^2)
di ""
di as txt "=================== RESULT (R=`R') ==================="
di as txt "true sampling sd of the impvalues point estimate : " as res %10.6f `sdb'
di as txt "mean reported impvalues (IM) standard error      : " as res %10.6f (`s_im'/`R')
di as txt "mean reported analytic clustered standard error  : " as res %10.6f (`s_an'/`R')
di as txt "impvalues 95% CI coverage                        : " as res %6.3f (`nc_im'/`R')
di as txt "analytic  95% CI coverage                        : " as res %6.3f (`nc_an'/`R')
di as txt "======================================================"
