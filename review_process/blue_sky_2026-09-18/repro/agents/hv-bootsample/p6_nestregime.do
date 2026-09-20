* p6_nestregime.do -- the T2/T3 nesting-regime decision (ado:8987-8998, ONE DAY OLD)
* is computed on `smp', which the code's own comment says must be "exactly the row set
* the leg keeps before it draws". It is not: smp holds rows e(sample) drops. A unit that
* straddles two clusters ONLY through a dropped row is therefore misclassified as
* non-nesting, and the comment names the consequence: the keep-i rule "would silently
* take the keep-i rule and collapse att_i's composition variance".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootsample"

clear
set obs 16
gen int state = _n
gen double sfx = rnormal()*0.5
gen int onset = .
replace onset = 2003 + mod(_n-1,4) if _n <= 8
expand 3
bysort state: gen int unit = (state-1)*3 + _n
gen double ufx = rnormal()*0.4
expand 10
bysort unit: gen int year = 1999 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .10*(year-1999) + sfx + ufx + (.3 + .4*sfx)*post + rnormal()*0.6
gen byte stray = 0
tempfile base
save `base', replace

* ONE stray record: unit 19 (state 7) is coded to state 8 for a single year, and that
* year's outcome is missing (a boundary change / crosswalk artefact on an unusable row).
* Inside e(sample) unit 19 is wholly inside state 7, so the cluster NESTS i there.
set obs `=_N+1'
replace unit  = 19      in l
replace state = 8       in l
replace year  = 2009    in l
replace onset = 2005    in l
replace y     = .       in l
replace stray = 1       in l
tempfile full
save `full', replace

qui count if stray
dis "N_STRAY=" r(N)

dis "=== N1: FULL DATA, bootstrap(standard) att_i ==="
cap noi wooldid y unit year onset, cluster(state) att_i bootstrap(standard) reps(300) rseed(21) dots(0)
dis "N1_RC=" _rc
dis "N1_N=" e(N)
mat n1b = e(b)
mat n1v = e(V)
dis "N1_NAMES:"
mat li n1b
dis "N1_SEDIAG:"
mata: st_matrix("sd1", sqrt(diagonal(st_matrix("n1v")))')
mat colnames sd1 = `: colnames n1b'
mat li sd1

dis "=== N2: the ONE stray row dropped (oracle) ==="
use `full', clear
drop if stray
cap noi wooldid y unit year onset, cluster(state) att_i bootstrap(standard) reps(300) rseed(21) dots(0)
dis "N2_RC=" _rc
dis "N2_N=" e(N)
mat n2b = e(b)
mat n2v = e(V)
dis "N2_SEDIAG:"
mata: st_matrix("sd2", sqrt(diagonal(st_matrix("n2v")))')
mat colnames sd2 = `: colnames n2b'
mat li sd2

dis "P6_ALLDONE"
