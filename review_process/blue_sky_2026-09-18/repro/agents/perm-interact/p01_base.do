clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
do "`W'/pchk.do"

* ---- 0. assert we are on the pinned build: nofastcollapse must parse ----
cap noi wooldid y unit t ttre , att cluster(unit) nofastcollapse
dis "BUILDASSERT_RC=" _rc

* ---- 1. determinism: same call twice ----
pchk det1
pchk det2

* ---- 2. baseline vs a few point-estimate-shaping options ----
pchk ctrl        , opts(controls(zbig))
pchk fe          , opts(fe(reg5))
pchk tt          , opts(timetrends)
pchk sg          , opts(subgroups(sg2))
pchk es          , opts(esprelength(2) espostlength(2))
pchk esf         , opts(esprelength(2) espostlength(2) esfixedbaseperiod)
pchk esrel       , opts(esprelength(2) espostlength(2) esrelativeto(-2))
pchk antic       , opts(anticipation(1))
pchk oos         , opts(oostreatedcontrols(2))
pchk robustopt   , opts(robust)
pchk nfc         , opts(nofastcollapse)
pchk cwm         , opts(customweightmultiplier(w))
pchk gt          , opts(gteffects)
dis "P01_DONE"
