clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* ================= three cheap analytic-path coverage arms =================
* A  power check : i=cg, cluster(cg)  (cluster AT i; documented to UNDER-cover)
* B  unconditionalse, att only, i=cg cluster(unit)
* C  impvalues at a NESTING cluster : i=unit, cluster(reg), regions nest units
*    -- this bounds the scope of the impvalues finding
local R 250

cap prog drop mkA
prog def mkA
  args sd
  clear
  set obs 60
  gen long unit = _n
  gen byte cg  = mod(_n-1,4)+1
  gen int  reg = ceil(_n/3)
  gen double ai = 0.9*mod(_n*13,17)/17 - 0.45
  expand 8
  bys unit: gen int t = 1999 + _n
  gen double ttre = .
  replace ttre = 2002 if cg==1
  replace ttre = 2004 if cg==2
  replace ttre = 2005 if cg==3
  gen double gt = 0.4*(t-1999)
  gen double truetau = 0
  replace truetau = 1.0 + 0.30*(t-ttre) if cg==1 & t>=ttre
  replace truetau = 0.6 + 0.15*(t-ttre) if cg==2 & t>=ttre
  replace truetau = 1.4 - 0.10*(t-ttre) if cg==3 & t>=ttre
  gen double base = ai + gt + truetau
  set seed `sd'
  set sortseed `sd'
  by unit: gen double ui = rnormal(0,0.40) if _n==1
  by unit: replace ui = ui[1]
  gen double y = base + ui + rnormal(0,0.60)
end

mkA 1
replace y = base
qui wooldid y cg t ttre, att cluster(cg) level(95)
local TA = e(b)[1,1]
qui wooldid y cg t ttre, att cluster(unit) unconditionalse level(95)
local TB = e(b)[1,1]
qui wooldid y unit t ttre, att cluster(reg) impvalues level(95)
local TC = e(b)[1,1]
dis "TRUTH|A=" %12.9f `TA' "|B=" %12.9f `TB' "|C=" %12.9f `TC'

tempname fh
file open `fh' using "$WD/mc_cheap.txt", write replace text
file write `fh' "draw|arm|b|se|lb|ub|truth" _n
forvalues r = 1/`R' {
  mkA `=600000+`r''
  cap qui wooldid y cg t ttre, att cluster(cg) level(95)
  if _rc==0 {
    mat Q=e(wdidmainresults)
    file write `fh' "`r'|A_clusATi|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TA') _n
  }
  cap qui wooldid y cg t ttre, att cluster(unit) unconditionalse level(95)
  if _rc==0 {
    mat Q=e(wdidmainresults)
    file write `fh' "`r'|B_uncond|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TB') _n
  }
  cap qui wooldid y unit t ttre, att cluster(reg) impvalues level(95)
  if _rc==0 {
    mat Q=e(wdidmainresults)
    file write `fh' "`r'|C_impnest|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TC') _n
  }
  if mod(`r',50)==0 dis "CHEAPPROG|`r'"
}
file close `fh'
dis "MC_CHEAP_DONE"
