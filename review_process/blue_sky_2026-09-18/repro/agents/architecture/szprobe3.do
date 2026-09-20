clear all
set more off
set processors 2
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/architecture/sz"
foreach k in 16 18 19 20 {
  cap prog drop _all
  clear mata
  cap noi run "`W'/wide2_`k'.ado"
  dis "SZ3 kb=`k' rc=" _rc
}
dis "SZ3 DONE"
