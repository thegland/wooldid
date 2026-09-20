*==============================================================================
* p16: semielasticity x permute().  In p02 the OBSERVED run completed (it printed
* its cell-collapse caution) and then the mandatory rep-0 run stopped the whole
* call with r(459) "the observed comparison statistic for att could not be read
* from the rep-0 run".  Is the plain (no-permute) run fine?  And does a strictly
* positive outcome fix it?
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
gen double ypos = exp(y/4) + 1
qui sum y
dis "S16|y_min=" r(min) "|y_max=" r(max)
qui sum ypos
dis "S16|ypos_min=" r(min) "|ypos_max=" r(max)

cap prog drop s16
prog def s16
  syntax anything(name=tag) [, OPTS(string asis) YV(string)]
  if "`yv'" == "" local yv y
  cap noi wooldid `yv' unit t ttre , att cluster(unit) `opts'
  dis "S16|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "S16|`tag'|bnames=`: colfullnames e(b)'"
  cap local bb = _b[main_att]
  dis "S16|`tag'|readable_main_att=" (_rc==0) "|value=" `=cond(_rc==0,"`bb'",".")'
end
dis "### y (crosses zero) ###"
s16 plain
s16 semi            , opts(semielasticity)
s16 semi_fn         , opts(semielasticity fastnose)
s16 semi_fn_nfc     , opts(semielasticity fastnose nofastcollapse)
s16 semi_perm       , opts(semielasticity permute(shuffle) permreps(30) rseed(777) dots(0))
dis "### ypos (strictly positive) ###"
s16 p_plain   , yv(ypos)
s16 p_semi    , yv(ypos) opts(semielasticity)
s16 p_semi_fn , yv(ypos) opts(semielasticity fastnose)
s16 p_semiperm, yv(ypos) opts(semielasticity permute(shuffle) permreps(30) rseed(777) dots(0))
dis "P16_DONE"
