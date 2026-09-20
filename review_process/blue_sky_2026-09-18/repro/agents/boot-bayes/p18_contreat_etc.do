* p18: the remaining bayes surfaces - continuous treatment (AME), BCa, if-sample,
* and the 2-cluster boundary.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

prog def ctr
  syntax , tag(string)
  matrix BB = e(b)
  matrix BD = e(bootdist)
  mata {
    B = st_matrix("BB"); D = st_matrix("BD"); nm = st_matrixcolstripe("BB")
    for (k=1; k<=cols(B); k++) {
      d = select(D[.,k], D[.,k] :< .)
      if (rows(d) > 2) {
        if (sqrt(variance(d)) > 0) printf("`tag' %-22s b=%12.6f  sd=%10.6f  (mean-b)/sd=%8.4f  n=%4.0f\n", nm[k,2], B[1,k], sqrt(variance(d)), (mean(d)-B[1,k])/sqrt(variance(d)), rows(d))
        else printf("`tag' %-22s b=%12.6f  ZERO DISPERSION over %4.0f draws\n", nm[k,2], B[1,k], rows(d))
      }
      else printf("`tag' %-22s b=%12.6f  <2 USABLE DRAWS (%4.0f)\n", nm[k,2], B[1,k], rows(d))
    }
  }
end

dis "=== A contreat + bayes (basic) ==="
use p16_panel, clear
cap noi wooldid y st year onset, att cluster(st) contreat(dose) contreatcontrols(basic) bootstrap(bayes) reps(150) rseed(77) dots(0)
dis "A_RC=" _rc
if _rc==0 ctr , tag(A)
dis "=== A2 the same, analytic, for the SE comparison ==="
use p16_panel, clear
cap noi wooldid y st year onset, att cluster(st) contreat(dose) contreatcontrols(basic)
dis "A2_RC=" _rc
if _rc==0 mat li e(wdidmainresults)

dis "=== B bayes + bootdisplay(bca) ==="
use p16_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootdisplay(bca) reps(200) rseed(77) dots(0)
dis "B_RC=" _rc
if _rc==0 {
  mat li e(boot_accel)
  mat li e(boot_z0)
  mat li e(boot_ci_bca)
  mat li e(boot_ci_percentile)
}

dis "=== C bayes under an if that halves the clusters ==="
use p16_panel, clear
cap noi wooldid y st year onset if st <= 15 | st > 22, att cluster(st) bootstrap(bayes) reps(100) rseed(77) dots(0) bootsave(p18_if) replace
dis "C_RC=" _rc
dis "C e(N)=" e(N) " e(N_clust)=" e(N_clust)
use p18_if, clear
qui ds w_*
local wv `r(varlist)'
dis "C number of w_ columns in bootsave = " `: word count `wv''
qui keep if _rep==1
qui egen double rs = rowtotal(`wv')
dis "C rowsum rep1 = " %21x rs[1] " = " rs[1]

dis "=== D the 2-cluster boundary ==="
use p16_panel, clear
qui keep if inlist(st, 1, 20)
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(50) rseed(77) dots(0)
dis "D_RC=" _rc
dis "=== E 1 cluster ==="
use p16_panel, clear
qui keep if st == 1
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(50) rseed(77) dots(0)
dis "E_RC=" _rc
dis "PROBE_DONE"
