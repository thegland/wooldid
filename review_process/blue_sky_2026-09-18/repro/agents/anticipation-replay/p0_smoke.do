clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid
local dN 60
local dT 10
include "`W'/mkdata.do"
cap prog drop wrap
prog def wrap
  syntax , w(string)
  qui wooldid y i t ttre, att robust espre(2) espost(2)
  dumpe , f("`w'/smoke.txt")
end
include "`W'/dumpe.do"
wrap , w("`W'")
dis "SMOKE_DONE"
