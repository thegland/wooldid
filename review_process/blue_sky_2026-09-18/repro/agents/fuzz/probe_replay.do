clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
use "`D'/ds_canon.dta", clear
dis "===== PLAIN RUN ====="
ereturn clear
qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2)
dis "RC=" _rc
dis "--- bare wooldid replay ---"
capture noisily wooldid
dis "BARE_RC=" _rc
dis "--- estimates store/replay ---"
capture noisily estimates store m1
dis "STORE_RC=" _rc
capture noisily estimates replay m1
dis "REPLAY_RC=" _rc
dis "--- second estimation then estimates restore ---"
qui wooldid y unit t ttre, att robust
capture noisily estimates restore m1
dis "RESTORE_RC=" _rc
dis "e(cmdline) after restore: `e(cmdline)'"
capture noisily wooldid
dis "BARE2_RC=" _rc
dis "--- estimates table / esttab-ish consumer ---"
capture noisily estimates table m1
dis "ETAB_RC=" _rc
