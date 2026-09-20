clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
gen byte regera = region*2 + (t>=2003)
gen double zt = zbig + 0.7*(t-2000)
do "`W'/pchk2.do"
local ES "esprelength(2) espostlength(2)"

pchk2 es_base   , opts(`ES')
pchk2 es_esf    , opts(`ES' esfixedbaseperiod)
pchk2 es_esf_rel, opts(`ES' esfixedbaseperiod esrelativeto(-2))
pchk2 es_sg     , opts(`ES' subgroups(sg2))
pchk2 es_sg_esf , opts(`ES' subgroups(sg2) esfixedbaseperiod)
pchk2 pre_sg    , opts(esprelength(2) subgroups(sg2))
pchk2 es_antic  , opts(`ES' anticipation(1))
pchk2 es_tt     , opts(`ES' timetrends)
pchk2 es_oos    , opts(`ES' oostreatedcontrols(2))
pchk2 es_ctrl   , opts(`ES' controls(zbig))
pchk2 es_ccc    , opts(`ES' coarsecohortcontrols(c.zt))
pchk2 es_fe     , opts(`ES' fe(regera))
pchk2 es_ic     , opts(`ES' interactivecontrols(x1))
pchk2 es_aw     , opts(`ES') wt([aw=w])
pchk2 es_fw     , opts(`ES') wt([fw=fwt])
pchk2 es_if     , opts(`ES') ifq(if keepreg==1)
pchk2 es_cwm    , opts(`ES' customweightmultiplier(w))
pchk2 es_nfc    , opts(`ES' nofastcollapse)
pchk2 es_all    , opts(`ES' att_i att_it att_itime)
pchk2 es_pois   , yv(ycount) opts(`ES' poisson)
* POSITIVE CONTROL: permanalyticse makes the observed run FULL while rep0 is fastnose,
* so the comparator must differ in the last digits. If this shows 0 too, the detector
* has no power and every clean row above is meaningless.
pchk2 es_ase    , opts(`ES' permanalyticse)
pchk2 es_ts     , opts(`ES' permtstat)
dis "P04_DONE"
