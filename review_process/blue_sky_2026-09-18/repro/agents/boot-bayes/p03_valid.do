* p03: bootalpha() validation surface + bayes gates
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
use p01_panel, clear

local n 0
foreach a in equal size EQUAL Size siz eq 1 equal2 "" {
  local ++n
  dis "=== T`n' bootalpha(`a') with bayes ==="
  cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(`a') reps(5) rseed(1) dots(0)
  dis "T`n'_RC=" _rc " alpha=" e(boot_alpha)
}
dis "=== S1 bootalpha with standard ==="
cap noi wooldid y st year onset, att cluster(st) bootstrap(standard) bootalpha(equal) reps(5) rseed(1) dots(0)
dis "S1_RC=" _rc
dis "=== S2 bootalpha with no bootstrap ==="
cap noi wooldid y st year onset, att cluster(st) bootalpha(equal)
dis "S2_RC=" _rc
dis "=== S3 boottstat + bayes ==="
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) boottstat reps(5) rseed(1) dots(0)
dis "S3_RC=" _rc
dis "=== S4 bayes + att_i default safety ==="
cap noi wooldid y st year onset, att_i cluster(st) bootstrap(bayes) reps(5) rseed(1) dots(0)
dis "S4_RC=" _rc
dis "=== S5 bayes + att_i safety(off) ==="
cap noi wooldid y st year onset, att_i cluster(st) bootstrap(bayes) reps(20) rseed(1) dots(0) safety(off)
dis "S5_RC=" _rc
mat li e(wdidmainresults)
dis "=== S6 bayes unclustered (observation-level) ==="
cap noi wooldid y st year onset, att bootstrap(bayes) reps(10) rseed(1) dots(0)
dis "S6_RC=" _rc
mat li e(wdidmainresults)
dis "=== S7 bayes + bootdisplay(bca) ==="
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootdisplay(bca) reps(20) rseed(1) dots(0)
dis "S7_RC=" _rc
dis "=== S8 bayes + poisson ==="
gen double ycount = max(0, round(exp(y/2)*3))
cap noi wooldid ycount st year onset, att cluster(st) poisson bootstrap(bayes) reps(10) rseed(1) dots(0)
dis "S8_RC=" _rc
dis "PROBE_DONE"
