clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
do "`W'/pchk.do"

pchk w_aw       , wt([aw=w])
pchk w_fw       , wt([fw=fwt])
pchk w_pw       , wt([pw=w])
pchk ifq        , ifq(if keepreg==1)
pchk lvl90      , opts(level(90))
pchk esti       , opts(att_i)
pchk estit      , opts(att_it)
pchk estitime   , opts(att_itime)
pchk estall     , opts(att_i att_it att_itime)
pchk ccc        , opts(coarsecohortcontrols(zbig))
pchk icx        , opts(interactivecontrols(x1))
pchk safoff     , opts(safety(off))
pchk regtol6    , opts(regtol(6))
pchk emc        , opts(emptycellsoverride)
pchk pois       , yv(ycount) opts(poisson)
pchk poisexp    , yv(ycount) opts(poisson poisexpresults)
pchk semi       , opts(semielasticity)
pchk jt         , opts(jointtests esprelength(2) espostlength(2))
pchk adv        , opts(adversarial esprelength(2) espostlength(2))
pchk sd         , opts(permstepdown esprelength(2) espostlength(2))
dis "P02_DONE"
