* repro_br2.do -- the reps() resolution Caution is unreachable without cluster().
* sthlp:534: "wooldid prints a Caution naming that grid whenever it is coarser than .05".
* ado:6455-6478 puts the whole block inside `if "`clus'" != "" {`.  Self-contained.
clear all
set more off
set processors 1
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 1
set sortseed 1
clear
set obs 40
gen long unit = _n
gen long ttre = cond(unit <= 12, 2003, cond(unit <= 24, 2005, .))
expand 8
bysort unit: gen int t = 1999 + _n
set seed 7
gen double y = 0.4*unit/40 + 0.1*t + 1.5*(t >= ttre) + rnormal()
dis "=== (1) CLUSTERED reps(12): the Caution IS printed ==="
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(12) rseed(9) dots(0)
dis "RC1=" _rc
dis "=== (2) UNCLUSTERED reps(12): identical resolution, NO Caution ==="
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(12) rseed(9) dots(0)
dis "RC2=" _rc
dis "=== (3) UNCLUSTERED reps(4) boottstat: floor 2/5 = .4, still NO Caution ==="
cap noi wooldid y unit t ttre, att bootstrap(standard) boottstat reps(4) rseed(9) dots(0)
dis "RC3=" _rc
dis "Grep this log for 'Caution: with reps(' -- it appears ONCE, under (1) only."
dis "REPRO_BR2_DONE"
