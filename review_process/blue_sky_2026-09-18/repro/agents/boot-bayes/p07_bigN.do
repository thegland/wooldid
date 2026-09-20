* p07: the Dirichlet draw matrix is a C x R STATA matrix. Unclustered bayes sets
*      C := number of observations. Does a realistically sized micro panel break it?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

dis "matsize=" c(matsize) "  max_matsize=" c(max_matsize)
foreach nn in 5000 11000 11001 15000 {
  cap noi mata: st_matrix("XX", J(`nn', 10, 1))
  dis "STMAT_`nn'_RC=" _rc
  cap mat drop XX
}

*--- realistic micro panel: 1200 individuals x 12 years = 14400 rows ----------
clear
set obs 1200
gen long pid = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2006 if pid <= 300
replace onset = 2009 if pid > 300 & pid <= 600
replace onset = 2012 if pid > 600 & pid <= 800
expand 12
bysort pid: gen int year = 2001 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.03*(year-2001) + 0.8*post + rnormal()
dis "N=" _N
save p07_panel, replace

dis "=== A: unclustered bayes on 14400 rows ==="
cap noi wooldid y pid year onset, att bootstrap(bayes) reps(10) rseed(5) dots(0)
dis "A_RC=" _rc

dis "=== B: same but 9000 rows (750 pid) ==="
use p07_panel, clear
qui keep if pid <= 750
dis "N=" _N
cap noi wooldid y pid year onset, att bootstrap(bayes) reps(10) rseed(5) dots(0)
dis "B_RC=" _rc

dis "=== C: standard bootstrap unclustered on 14400 (control) ==="
use p07_panel, clear
cap noi wooldid y pid year onset, att bootstrap(standard) reps(10) rseed(5) dots(0)
dis "C_RC=" _rc
dis "PROBE_DONE"
