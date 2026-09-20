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

* 1. a good run at the default level
set level 95
wooldid y i t ttre, att robust espre(2) espost(2)
dis "G1_N=" e(N) " level=" e(level)
dis "G1_LASTRC=[$WOOLDID_LASTRC] LEVEL0=[$WOOLDID_LEVEL0]"

* 2. now the user changes the session level and re-runs -- the call FAILS
set level 99.9
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
dis "F_rc=" _rc
dis "F_LASTRC=[$WOOLDID_LASTRC]  <-- documented to be set after ANY failed call"
dis "F_LEVEL0=[$WOOLDID_LEVEL0]  <-- documented to be dropped on every exit path"
dis "F_clevel=" c(level)

* 3. what does a bare `wooldid' now show?  sthlp:1460 says it must REFUSE.
dis as result "==== bare wooldid after the failed call ===="
cap noi wooldid
dis "REPLAY_rc=" _rc
dis "REPLAY_shows_e(N)=" e(N) "  e(level)=" e(level)

* 4. does the leaked global break the NEXT call?
set level 95
cap noi wooldid y i t ttre, att_it robust espre(2) espost(2)
dis "NEXT_rc=" _rc " e(level)=" e(level)
dis "NEXT_LEVEL0=[$WOOLDID_LEVEL0]"
dis "P14_DONE"
