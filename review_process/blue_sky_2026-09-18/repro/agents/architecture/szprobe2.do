clear all
set more off
set processors 2
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/architecture/sz"
foreach k in 17 18 19 20 22 24 26 28 30 {
  cap prog drop _all
  clear mata
  cap noi run "`W'/w_`k'.ado"
  dis "SZPROBE2 kb=`k' rc=" _rc
}
dis "SZPROBE2 DONE"
