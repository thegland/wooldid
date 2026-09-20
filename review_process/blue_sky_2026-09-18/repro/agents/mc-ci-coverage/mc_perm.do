clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* usage: -e do mc_perm.do <chunk> <ndraws> <permreps>
local CH `1'
local R  `2'
local PR `3'
if "`CH'"=="" local CH 1
if "`R'"==""  local R 5
if "`PR'"=="" local PR 99

cap prog drop mkP
prog def mkP
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
  gen double truetau = cond(ttre<. & t>=ttre, 1, 0)
  gen double base = ai + gt + truetau
  set seed `sd'
  set sortseed `sd'
  by unit: gen double ui = rnormal(0,0.40) if _n==1
  by unit: replace ui = ui[1]
  gen double y = base + ui + rnormal(0,0.60)
end
mkP 1
replace y = base
qui wooldid y cg t ttre, att cluster(unit) level(95)
local TRUTH = e(b)[1,1]
dis "TRUTH|" %14.9f `TRUTH'
tempname fh
file open `fh' using "$WD/mcperm_`CH'.txt", write replace text
file write `fh' "draw|arm|b|lb|ub|truth" _n
forvalues r = 1/`R' {
  local dr = (`CH'-1)*1000 + `r'
  mkP `=1000000+`dr''
  cap qui wooldid y cg t ttre, att cluster(unit) permute(shuffle) permreps(`PR') ///
       rseed(`=1100000+`dr'') permtstat permtci permonsetcohorts level(95) dots(0)
  if _rc==0 {
    mat Q = e(wdidmainresults)
    file write `fh' "`dr'|permtab|" (Q[1,1]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TRUTH') _n
  }
  cap qui wooldid y cg t ttre, att cluster(unit) level(95)
  if _rc==0 {
    mat Q = e(wdidmainresults)
    file write `fh' "`dr'|analytic|" (Q[1,1]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TRUTH') _n
  }
  dis "MCP|`dr'|done"
}
file close `fh'
dis "MCPERM_DONE_`CH'"
