* probe16: pinpoint WHERE a weight expression fails.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
use "p8.dta", clear
set trace on
set tracedepth 4
set traceexpand on
cap noi wooldid y i t ttre [aw=ka+kb], att
set trace off
dis "WHERE_RC=" _rc
dis "PROBE16_DONE"
