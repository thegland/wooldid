*==============================================================================
* p08: does the permutation NULL DISTRIBUTION depend on safety()/screens?
* A replication whose reassigned onsets make the design inestimable ERRORS and is
* silently dropped from the reference set (a caution fires only above 10% failures).
* On a panel with staggered ENTRY, a unit that enters late and draws an early onset
* has no usable pre-period, so a non-trivial share of reassignments is inestimable.
* If safety(on) and safety(off) drop DIFFERENT replications, the reported p-value
* depends on a screen that is about estimation quality, not about the null.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local gseed   771
local gN      30
local gT      8
local gt0     2000
local gonsets 2003 2005
local gysd    0.30
local gentry  0.45
local gentrymax 3
local gexit   0.25
local gexitmax 2
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
save "p08_data.dta", replace
tab t ttre, missing

cap prog drop schk
prog def schk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(unit) permute(shuffle) permreps(200) rseed(31337) dots(0) `opts'
  dis "SAFE|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "SAFE|`tag'|reps=" e(perm_reps) "|done=" e(perm_completed) "|fail=" e(perm_failed) "|p_att=" %21x e(perm_p_att) "|sd=" %21x e(perm_sd_att) "|b=" %21x _b[main_att]
end
schk safety_on
schk safety_off , opts(safety(off))
schk emc_on     , opts(emptycellsoverride)
schk es_on      , opts(esprelength(2) espostlength(2))
schk es_off     , opts(esprelength(2) espostlength(2) safety(off))
dis "P08_DONE"
