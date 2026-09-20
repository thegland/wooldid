*==============================================================================
* p17: pin the mechanism behind p16.  Is it the SIGN of the outcome?
* y            : crosses zero      (min -3.70, max  6.86)
* yshift       : y + 10            (strictly positive, identical shape)
* ysmall       : y + 3.8           (positive but predictions land near zero)
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
gen double yshift = y + 10
gen double ysmall = y + 3.8
foreach v in y yshift ysmall {
  qui sum `v'
  dis "S17|`v'|min=" r(min) "|max=" r(max)
}
cap prog drop s17
prog def s17
  syntax anything(name=tag) [, OPTS(string asis) YV(string)]
  cap noi wooldid `yv' unit t ttre , att cluster(unit) `opts'
  local rc = _rc
  cap confirm matrix e(b)
  local hb = (_rc==0)
  dis "S17|`tag'|RC=" `rc' "|HAS_e(b)=" `hb' "|bnames=`: colfullnames e(b)'"
end
foreach v in y yshift ysmall {
  s17 `v'_plain    , yv(`v')
  s17 `v'_semi     , yv(`v') opts(semielasticity)
  s17 `v'_semifn   , yv(`v') opts(semielasticity fastnose)
  s17 `v'_semiperm , yv(`v') opts(semielasticity permute(shuffle) permreps(30) rseed(777) dots(0))
}
* is an EMPTY e(b) reachable under fastnose WITHOUT semielasticity?
s17 fn_only , yv(y) opts(fastnose)
dis "P17_DONE"
