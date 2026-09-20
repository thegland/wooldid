clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
use "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-fastnose-doc/panel1.dta", clear
tab ttre, m
tab t if missing(ttre)
dis "--- plain att ---"
cap noi wooldid y i t ttre, att
dis "P1_RC=" _rc
dis "--- att + robust ---"
cap noi wooldid y i t ttre, att robust
dis "P2_RC=" _rc
dis "--- att cluster(i) ---"
cap noi wooldid y i t ttre, att cluster(i)
dis "P3_RC=" _rc
dis "--- att cluster(i) es 3/3 ---"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3)
dis "P4_RC=" _rc
dis "--- att cluster(i) es 3/3 esfixedbaseperiod ---"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod
dis "P5_RC=" _rc
dis "=== DIAG DONE ==="
