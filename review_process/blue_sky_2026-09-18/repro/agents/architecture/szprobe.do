clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/architecture/sz"
foreach k in 0 1 2 4 8 16 32 {
  cap prog drop _all
  clear mata
  cap noi run "`W'/w_`k'.ado"
  dis "SZPROBE kb=`k' rc=" _rc
}
dis "SZPROBE DONE"
