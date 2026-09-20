* probe10 -- is ANYTHING printed on the CCC arm?  Non-quiet, all three runs.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
which wooldid
dis _n "@@@@ 1: no permute, no ccc"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state)
dis "R1=" _rc
dis _n "@@@@ 2: no permute, ccc(c.cohbase)"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state) coarsecohortcontrols(c.cohbase)
dis "R2=" _rc
dis _n "@@@@ 3: permute + permonsetcohorts + ccc(c.cohbase)  -- the FP1 arm, NON-QUIET"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts coarsecohortcontrols(c.cohbase)
dis "R3=" _rc
dis _n "@@@@ 4: permute + permonsetcohorts, no ccc -- NON-QUIET (the comparison)"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts
dis "R4=" _rc
dis "PROBE10 DONE"
