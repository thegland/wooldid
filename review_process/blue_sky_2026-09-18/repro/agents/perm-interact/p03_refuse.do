clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear

cap prog drop rchk
prog def rchk
  syntax anything(name=tag) [, OPTS(string asis) WT(string asis) IFQ(string asis) YV(string) CL(string)]
  if "`yv'" == "" local yv y
  if "`cl'" == "" local cl cluster(unit)
  cap noi wooldid `yv' unit t ttre `wt' `ifq' , att `cl' permute(shuffle) permreps(30) rseed(777) `opts'
  dis "RCHK|`tag'|RC=" _rc
end

rchk impv         , opts(impvalues)
rchk uncond       , opts(unconditionalse robust)
rchk cvce_nocl    , opts(customvce(cluster unit)) cl( )
rchk cvce_cl      , opts(customvce(cluster unit))
rchk ptp          , opts(pretrendpower(0.5) esprelength(2) espostlength(2))
rchk ssg          , opts(suppresssgprimaryeffects subgroups(sg2))
rchk oldsyn       , opts(oldsyntax)
rchk simband      , opts(simultaneousbands esprelength(2) espostlength(2))
rchk ctw          , opts(contreat(dose) contreatwithin)
rchk ctuz         , opts(contreat(dose) contreatcontrols(untreatedZero))
rchk multiclus    , cl(cluster(unit region))
rchk fn_ts        , opts(fastnose permtstat)
rchk fastnose     , opts(fastnose)
rchk permtstat    , opts(permtstat)
rchk permase      , opts(permanalyticse)
rchk permtci      , opts(permtci)
rchk ctreat       , opts(contreat(dose))
rchk ctreat_basic , opts(contreat(dose) contreatcontrols(basic))
rchk gtplus       , opts(gteffects)
rchk makeplots    , opts(makeplots esprelength(2) espostlength(2))
rchk permonset    , opts(permonsetcohorts)
dis "P03_DONE"
