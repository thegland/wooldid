clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local P "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-fastnose-doc"

dis "=== T1: renamed vars, only y/unit/year/onset kept ==="
use "`P'/panel1.dta", clear
rename i unitid
rename t year
rename ttre onset
keep y unitid year onset
cap noi wooldid y unitid year onset, att
dis "T1_RC=" _rc

dis "=== T2: original names, only 4 vars kept ==="
use "`P'/panel1.dta", clear
keep y i t ttre
cap noi wooldid y i t ttre, att
dis "T2_RC=" _rc

dis "=== T3: renamed, all vars kept ==="
use "`P'/panel1.dta", clear
rename i unitid
rename t year
rename ttre onset
cap noi wooldid y unitid year onset, att
dis "T3_RC=" _rc
dis "=== NAMES DONE ==="
