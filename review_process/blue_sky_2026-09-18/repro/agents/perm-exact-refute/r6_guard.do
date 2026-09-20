* r6_guard.do -- MEASURE (not assume) the guard that closes the one structural route to
* PERMEX-PX-2's exk=1 & memb=0 state: a cluster carrying more than one treatment onset.
* If that were allowed, the rep-0 comparator (which fits on ttrevar as given) and the
* identity replication (which writes ONE cluster-level onset to every row of the cluster)
* would fit different data, so the observed statistic would not be a member of its own
* reference set while the enumeration stayed complete.
* Design: cluster(state), i=unit; state 1 holds units 1 and 2 with onsets 4 and 6.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 8
gen long unit = _n
gen byte state = ceil(unit/2)          /* states 1..4, two units each */
gen double onset = .
replace onset = 4 if unit == 1
replace onset = 6 if unit == 2         /* SAME state as unit 1, DIFFERENT onset */
replace onset = 4 if unit == 3
replace onset = 4 if unit == 4
gen double ue = rnormal()*0.4
expand 9
bysort unit: gen int year = _n
gen double te = rnormal()*0.3
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.12*year + ue + te + 0.9*post
tempfile d
save `d', replace

dis "===== G1: mixed-onset cluster, default safety ====="
capture noisily wooldid y unit year onset, att espostlength(2) cluster(state) ///
  permute(shuffle) permexact dots(0)
dis "G1_RC=" _rc

dis "===== G2: same, safety(off) ====="
use `d', clear
capture noisily wooldid y unit year onset, att espostlength(2) cluster(state) ///
  permute(shuffle) permexact safety(off) dots(0)
dis "G2_RC=" _rc

dis "===== G3: control -- make the state onset-homogeneous, must RUN ====="
use `d', clear
replace onset = 4 if unit == 2
capture noisily wooldid y unit year onset, att espostlength(2) cluster(state) ///
  permute(shuffle) permexact dots(0)
dis "G3_RC=" _rc
dis "G3_perm_exact=" e(perm_exact)
dis "G3_perm_reps=" e(perm_reps)
dis "DONE_R6"
