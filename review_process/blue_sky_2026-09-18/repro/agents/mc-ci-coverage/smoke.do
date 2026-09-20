clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

* ---- build a fixed balanced panel: 120 units in 4 cohorts (3 treated + never) ----
local NU = 120
local NT = 8
set obs `NU'
gen long unit = _n
gen byte cg = mod(_n-1,4)+1          // 1,2,3 treated ; 4 never
gen double ai = rnormal(0,1)
expand `NT'
bys unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 if cg==1
replace ttre = 2004 if cg==2
replace ttre = 2005 if cg==3
gen double gt = 0.4*(t-1999)
* true effect: cohort-and-horizon specific
gen double truetau = 0
replace truetau = 1.0 + 0.30*(t-ttre) if cg==1 & t>=ttre
replace truetau = 0.6 + 0.15*(t-ttre) if cg==2 & t>=ttre
replace truetau = 1.4 - 0.10*(t-ttre) if cg==3 & t>=ttre
gen double base = ai + gt + truetau
save fixture, replace

* ---- design truth at i = cg ----
gen double e = 0
gen double y = base
timer clear
timer on 1
wooldid y cg t ttre, att att_i att_it att_itime cluster(unit)
timer off 1
timer list 1
mat list e(wdidmainresults)
dis "LEVEL=" e(level)
