clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* impvalues at the pairing the help RECOMMENDS: i = cohort, cluster = the unit
* (cluster FINER than i).  Mildly heterogeneous cohort effects.
local R 200
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
  gen double gt = 0.4*(t-1999)
  gen double truetau = 0
  replace truetau = 1.00 if cg==1 & t>=ttre
  replace truetau = 1.05 if cg==2 & t>=ttre
  replace truetau = 0.95 if cg==3 & t>=ttre
  gen double base = ai + gt + truetau
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
qui wooldid y cg t ttre, att cluster(unit) level(95)
local TA = e(b)[1,1]
dis "TRUTH2|im=" %14.9f `TI' "|an=" %14.9f `TA'
tempname fh
file open `fh' using "$WD/mc_imp2.txt", write replace text
file write `fh' "draw|arm|b|se|lb|ub|truth" _n
forvalues r = 1/`R' {
  mkQ `=1500000+`r''
  cap qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
  if _rc==0 {
    mat Q=e(wdidmainresults)
    file write `fh' "`r'|im_clusfiner|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TI') _n
  }
  cap qui wooldid y cg t ttre, att cluster(unit) level(95)
  if _rc==0 {
    mat Q=e(wdidmainresults)
    file write `fh' "`r'|analytic|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TA') _n
  }
}
file close `fh'
dis "MC_IMP2_DONE"
