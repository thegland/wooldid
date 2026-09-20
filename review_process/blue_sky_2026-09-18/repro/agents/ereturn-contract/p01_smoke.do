clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vio_smoke.txt"
cap erase "vio_smoke.txt"
include "mkdata.doi"
include "ckc.doi"
mkD1
mkD2
mkD3
mkD4
dis "DATA_OK"

use D1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime
dis "S1_RC=" _rc
ckc S1

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3)
dis "S2_RC=" _rc
ckc S2

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
dis "S3_RC=" _rc
ckc S3

use D3.dta, clear
cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic)
dis "S4_RC=" _rc
ckc S4

dis "SMOKE_DONE"
