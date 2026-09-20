clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid

*==================================================================
* A. bootdf(clusters): is G the estimation sample's cluster count?
*    Realistic: three states whose outcome is suppressed (y all missing).
*==================================================================
local gG 20
local gUnits 3
include "`WD'/gen_small.do"
save base20, replace

gen double ysup = y
replace ysup = . if inlist(state, 4, 9, 14)   /* 3 states suppressed */
qui count if !missing(ysup)
dis "A_NONMISS=" r(N)

dis "=== A1: full y (control) ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(99) bootdf(clusters) rseed(11) att dots(0)
dis "A1_RC=" _rc
dis "A1_NCLUST=" e(N_clust) " BOOTDF=" e(boot_df) " DFR=" e(df_r)

dis "=== A2: three states with y suppressed ==="
cap noi wooldid ysup state t ttre , cluster(state) bootstrap(standard) reps(99) bootdf(clusters) rseed(11) att dots(0)
dis "A2_RC=" _rc
dis "A2_NCLUST=" e(N_clust) " BOOTDF=" e(boot_df) " DFR=" e(df_r)
dis "A2_COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed)
mat A2 = e(wdidmainresults)
mat list A2, format(%21x)

dis "=== A3: same data, dead states physically removed ==="
preserve
drop if inlist(state, 4, 9, 14)
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(99) bootdf(clusters) rseed(11) att dots(0)
dis "A3_RC=" _rc
dis "A3_NCLUST=" e(N_clust) " BOOTDF=" e(boot_df) " DFR=" e(df_r)
mat A3 = e(wdidmainresults)
mat list A3, format(%21x)
restore

*==================================================================
* B. WOOLDID_BOOTNOSUB
*==================================================================
use base20, clear
dis "=== B1: boottstat, substitution ON (normal) ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(99) boottstat rseed(777) att dots(0)
dis "B1_RC=" _rc
mat B1S = e(bootdist_se)
mat B1A = e(boot_ase_obs)
mata:
  S = st_matrix("B1S"); A = st_matrix("B1A")
  for (k=1; k<=cols(S); k++) {
    v = select(S[.,k], S[.,k] :< .)
    printf("B1 col %f ratio_mean=%9.6f ratio_rms=%9.6f  n=%f\n", k, mean(v)/A[1,k], sqrt(mean(v:^2))/A[1,k], rows(v))
  }
end

dis "=== B2: same seed, WOOLDID_BOOTNOSUB set ==="
global WOOLDID_BOOTNOSUB 1
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(99) boottstat rseed(777) att dots(0)
dis "B2_RC=" _rc
mat B2S = e(bootdist_se)
mat B2A = e(boot_ase_obs)
mat B2M = e(wdidmainresults)
mata:
  S = st_matrix("B2S"); A = st_matrix("B2A")
  for (k=1; k<=cols(S); k++) {
    v = select(S[.,k], S[.,k] :< .)
    printf("B2 col %f ratio_mean=%9.6f ratio_rms=%9.6f  n=%f\n", k, mean(v)/A[1,k], sqrt(mean(v:^2))/A[1,k], rows(v))
  }
end
dis "B2_GLOBAL_STILL=" "${WOOLDID_BOOTNOSUB}"

dis "=== B3: does the global survive into the NEXT run untouched? ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(99) boottstat rseed(777) att dots(0)
dis "B3_RC=" _rc
dis "B3_GLOBAL=" "${WOOLDID_BOOTNOSUB}"
mat B3M = e(wdidmainresults)
mat list B2M, format(%21x)
mat list B3M, format(%21x)
global WOOLDID_BOOTNOSUB

*==================================================================
* C. bootdisplay(bc) + boottstat must be refused; bca without cluster refused
*==================================================================
dis "=== C1 ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(20) boottstat bootdisplay(bc) att dots(0)
dis "C1_RC=" _rc
dis "=== C2 ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(20) boottstat bootdisplay(wald) att dots(0)
dis "C2_RC=" _rc
dis "=== C3 (bca no cluster) ==="
cap noi wooldid y state t ttre , bootstrap(standard) reps(20) bootdisplay(bca) att dots(0)
dis "C3_RC=" _rc
dis "=== C4 (bca + boottstat) ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(20) boottstat bootdisplay(bca) att dots(0)
dis "C4_RC=" _rc
dis "P03_DONE"
