* p19: bracket the rc-0-with-no-e(b) case for semielasticity+fastnose.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"

cap prog drop b19
prog def b19
  syntax anything(name=tag) [, YV(string) OPTS(string asis)]
  cap noi wooldid `yv' unit t ttre , att cluster(unit) semielasticity fastnose `opts'
  local rc = _rc
  cap confirm matrix e(b)
  local hb = (_rc==0)
  dis "B19|`tag'|RC=" `rc' "|HAS_e(b)=" `hb' "|e(cmd)=`e(cmd)'"
end

* ---- control: the known rc-0 case, base.dta's y ----
use "`W'/base.dta", clear
qui sum y
dis "B19|base_y|min=" r(min) "|max=" r(max) "|mean=" r(mean)
b19 base_y , yv(y)
b19 base_ync , yv(ync)

* ---- inline variants ----
foreach nz in 0.6 1.5 {
 foreach nc in 2 4 {
  clear
  set seed 20260916
  set sortseed 20260916
  set obs 40
  gen long unit = _n
  gen double ttre = .
  if `nc' == 2 {
    replace ttre = 2003 if unit <= 12
    replace ttre = 2005 if unit > 12 & unit <= 24
  }
  else {
    replace ttre = 2002 if unit <= 6
    replace ttre = 2003 if unit > 6  & unit <= 12
    replace ttre = 2004 if unit > 12 & unit <= 18
    replace ttre = 2005 if unit > 18 & unit <= 24
  }
  expand 8
  bysort unit: gen int t = 1999 + _n
  gen double a_i  = rnormal()
  gen double g_t  = 0.4*rnormal()
  bysort t (unit): replace g_t = g_t[1]
  gen double x1 = rnormal()
  gen double x2 = 0.3*a_i + 0.5*rnormal()
  gen double tau = 2.0 + rnormal()*0.5
  bysort unit (t): replace tau = tau[1]
  gen double ync = a_i + g_t + tau*(ttre<. & t>=ttre) + rnormal()*`nz'
  gen double y   = ync + 0.5*x1 - 0.25*x2
  qui sum y
  dis "B19|nz`nz'_nc`nc'|ymin=" r(min) "|ymax=" r(max)
  b19 nz`nz'_nc`nc'_y   , yv(y)
  b19 nz`nz'_nc`nc'_ync , yv(ync)
 }
}
dis "P19_DONE"
