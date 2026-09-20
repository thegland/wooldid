*==============================================================================
* p13: continuous treatment (contreat / the ame family) x permute().
* permctl deliberately ALLOWS contreat() ("the null is timing-as-if-random given
* dose") and refuses only contreatcontrols(untreatedZero) and contreatwithin.
* i = cohort, cluster = cohort, dose varying at the unit level.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local gseed   918
local gN      60
local gT      8
local gt0     2000
local gonsets 2003 2004 2005 2006
local gnevgrp 2
local gysd    0.30
local gdose   unitconst
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
save "p13_data.dta", replace
qui levelsof cohort, local(cc)
dis "P13_NCOHORT=" wordcount("`cc'")
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact/pchk2.do"

cap prog drop cchk
prog def cchk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y cohort t ttre , att cluster(cohort) contreat(dose) permute(shuffle) permreps(30) rseed(777) dots(0) `opts'
  dis "C13|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "C13|`tag'|cfx_blank=`e(perm_cfx_blank)'|gt_blank=`e(perm_gt_blank)'|blank=" e(perm_blank) "|done=" e(perm_completed) "|fail=" e(perm_failed)
  dis "C13|`tag'|bnames=`: colfullnames e(b)'"
  cap confirm matrix e(perm_comparator_b)
  if _rc == 0 mata: wdpcmp2("`tag'")
end
cchk basic        , opts(contreatcontrols(basic))
cchk basic_es     , opts(contreatcontrols(basic) esprelength(2) espostlength(2))
cchk bytime       , opts(contreatcontrols(byTime))
cchk basic_lat    , opts(contreatcontrols(basic) lattice(all 3))
cchk basic_cfx    , opts(contreatcontrols(basic) makecfxplots)
cchk basic_cfxbt  , opts(contreatcontrols(basic) makecfxplots cfxplottypes(post) lattice(byTtre 3))
cchk basic_sg     , opts(contreatcontrols(basic) subgroups(sg2))
cchk basic_semi   , opts(contreatcontrols(basic) contreatelasticitytype(dyex))
cchk basic_poly   , opts(contreatcontrols(basic) contreatpoly(2))
cchk basic_ts     , opts(contreatcontrols(basic) permtstat)
* the two refusals, for the record
cap noi wooldid y cohort t ttre , att cluster(cohort) contreat(dose) contreatcontrols(untreatedZero) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "C13|uz|RC=" _rc
cap noi wooldid y cohort t ttre , att cluster(cohort) contreat(dose) contreatcontrols(basic) contreatwithin permute(shuffle) permreps(30) rseed(777) dots(0)
dis "C13|cw|RC=" _rc
dis "P13_DONE"
