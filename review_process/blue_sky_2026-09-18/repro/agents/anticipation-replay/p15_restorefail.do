clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid
local dN 50
local dT 10
local donsets "2004 2006"
include "`W'/mkdata.do"

* The user's session level is 99.9; they ask for a 95% interval on ONE call.
set level 99.9
dis "BEFORE c(level)=" c(level)
cap noi wooldid y i t ttre, att robust espre(2) espost(2) level(95)
dis "RC=" _rc
dis "AFTER c(level)=" c(level)   "   <-- the caller's session level"
dis "e(N)=" e(N) "  e(level)=" e(level) "  e(cmd)=[" e(cmd) "]"
dis "LEVEL0=[$WOOLDID_LEVEL0]  LASTRC=[$WOOLDID_LASTRC]  CLEANUP=[$WOOLDID_CLEANUP]"
cap est dir
dis "---- can the user replay? ----"
cap noi wooldid
dis "REPLAY_rc=" _rc
set level 95
dis "P15_DONE"
