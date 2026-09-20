* r12 (scope control for r11, SELF-CONTAINED): are the e(gteffects) cell standard errors, p-values
* and confidence intervals right at cluster() at-or-above i?
* SAME data as r11 but i = the ONSET COHORT (the help's intended practice: cluster
* BELOW the level of i), cluster(st).  Each treated (i,relyr) cell now holds 7-8 states.
* 60 Monte Carlo draws with a HOMOGENEOUS true effect of 1.  For the cell (i = 1,
* relyr = 0) and for att, compare the MEAN REPORTED standard error against the TRUE
* sampling standard deviation of the estimate across the draws, and count how often the
* cell's own reported 95% interval covers the cell's true effect (= 1).
* Everything is computed inside this file; no external script is needed.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

tempname M
postfile `M' double(b se cest cse clo chi) using r12_res, replace
forvalues d = 1/60 {
  set seed `=310000 + `d'*37'
  set sortseed `=310000 + `d'*37'
  clear
  set obs 30
  gen int st = _n
  gen int grp = cond(st<=8,1,cond(st<=15,2,cond(st<=21,3,4)))
  gen double a_i = rnormal()
  gen int onset = .
  replace onset = 2004 if st <=  8
  replace onset = 2006 if st >  8 & st <= 15
  replace onset = 2008 if st > 15 & st <= 21
  expand 10
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.03*(year-2000) + 1.0*post + rnormal()*0.7
  cap noi qui wooldid y grp year onset, att cluster(st) gteffects
  if _rc == 0 {
    matrix A = e(wdidmainresults)
    matrix G = e(gteffects)
    * e(gteffects) columns: 1 i, 2 onset, 3 caltime, 4 relyr, 5 estimate, 6 se, 7 p, 8 lb, 9 ub
    local ce = .
    local cs = .
    local cl = .
    local ch = .
    forvalues r = 1/`=rowsof(G)' {
      if G[`r',1] == 1 & G[`r',4] == 0 {
        local ce = G[`r',5]
        local cs = G[`r',6]
        local cl = G[`r',8]
        local ch = G[`r',9]
      }
    }
    post `M' (A[1,1]) (A[1,2]) (`ce') (`cs') (`cl') (`ch')
  }
  else dis "FAIL d=`d' rc=" _rc
}
postclose `M'
use r12_res, clear
qui count
dis "DRAWS=" r(N)
qui su b
local sdb = r(sd)
qui su se
local mse = r(mean)
qui su cest
local sdc = r(sd)
qui su cse
local msc = r(mean)
qui gen byte cov = (clo <= 1 & 1 <= chi)
qui su cov
local cvr = r(mean)
qui su cse
local mincse = r(min)
local maxcse = r(max)
dis "att  : TRUE sd=" %9.5f `sdb' "   mean reported se=" %9.5f `mse' "   se/sd=" %7.3f `mse'/`sdb'
dis "cell : TRUE sd=" %9.5f `sdc' "   mean reported se=" %9.5f `msc' "   se/sd=" %7.3f `msc'/`sdc'
dis "cell reported se range across draws: " %9.5f `mincse' " to " %9.5f `maxcse'
dis "cell 95% interval coverage of its own true effect (=1): " %6.3f `cvr' "  (nominal 0.95)"
dis "R12_DONE"
