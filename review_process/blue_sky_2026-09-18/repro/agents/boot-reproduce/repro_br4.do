* repro_br4.do -- SELF-CONTAINED.
* The bootstrap's resampling universe is `smp' (user if & nonmissing cluster), not e(sample).
* Units whose outcome is entirely missing are counted as resampling clusters.  Consequences:
*   (1) the few-cluster CAUTIONS documented at sthlp:532 are computed on the inflated count and
*       go silent on a run that really does have fewer than 10 (and fewer than 6 treated) clusters;
*   (2) e(boot_df) -- the default bootdf(clusters) = G-1 -- is one too large per phantom cluster;
*   (3) the run stops matching Stata's own bootstrap prefix replication for replication,
*       which sthlp:532 says it does.
* ARM A leaves the all-missing-outcome units in the data; ARM B drops them.  Same estimation
* sample (identical e(N), identical e(N_clust)), different cautions and different e(boot_df).
clear all
set more off
set processors 1
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* --- 14 units x 8 periods, 5 treated in 2003, 4 in 2005, 5 never treated ---
clear
set obs 14
gen long unit = _n
gen double ttre = cond(unit <= 5, 2003, cond(unit <= 9, 2005, .))
expand 8
bysort unit: gen int t = 1999 + _n
set seed 31
gen double a_i = rnormal()
gen double y = a_i + 0.08*(t - 2000) + 1.2*(t >= ttre) + 0.4*rnormal()
* five units -- three treated, two never treated -- have NO usable outcome at all.
* This is an ordinary administrative shape: the unit is in the frame, the outcome is not.
gen byte phantom = inlist(unit, 1, 2, 3, 6, 7)
save "br4_full.dta", replace
replace y = . if phantom
save "br4_arma.dta", replace
use "br4_full.dta", clear
drop if phantom
save "br4_armb.dta", replace

dis "=================== ARM A: phantom units present, outcome missing ==================="
use "br4_arma.dta", clear
qui levelsof unit, local(uu)
dis "ARM A: units in the data              = " wordcount("`uu'")
qui levelsof unit if !missing(y), local(uu2)
dis "ARM A: units with any usable outcome  = " wordcount("`uu2'")
qui levelsof unit if !missing(y) & !missing(ttre), local(uu3)
dis "ARM A: TREATED units with an outcome  = " wordcount("`uu3'")
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(4) dots(0)
dis "ARMA_RC=" _rc
dis "ARMA| e(N)=" e(N) " e(N_clust)=" e(N_clust) " e(boot_df)=" e(boot_df) "  [boot_df = G-1, so G = e(boot_df)+1]"

dis "=================== ARM B: the same estimation sample, phantoms dropped ==================="
use "br4_armb.dta", clear
qui levelsof unit, local(uu)
dis "ARM B: units in the data              = " wordcount("`uu'")
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(4) dots(0)
dis "ARMB_RC=" _rc
dis "ARMB| e(N)=" e(N) " e(N_clust)=" e(N_clust) " e(boot_df)=" e(boot_df) "  [boot_df = G-1, so G = e(boot_df)+1]"

dis "=================== WHAT TO LOOK FOR ==================="
dis "Both arms estimate on the SAME 9 units / same e(N) / same e(N_clust)."
dis "Grep this log for: 'clusters are available for the bootstrap'  (the <10 caution)"
dis "                   'treated cluster(s) are available'          (the <6 treated caution)"
dis "They should appear under BOTH arms.  If they appear only under ARM B, the cautions are"
dis "computed on the resampling universe, which counts units the estimation never sees."
dis "REPRO_BR4_DONE"
