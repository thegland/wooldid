* p13: sthlp:557 says the bootalpha(size) _rep=0 row of a bootsave() file holds
* "1/N, N the number of observations in the ESTIMATION SAMPLE". ado:437 divides by
* `wdb_r_nrows', which wooldid_bootrun sets to the MAP-FRAME row count - the user-if
* sample intersected with nonmissing cluster(), which still contains rows the
* estimator drops (missing outcome, singletons, ...). Also: under bootalpha(size)
* the Dirichlet SHAPE alpha_c = n_c counts those same dropped rows.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

use p01_panel, clear
set seed 55
* 15% item nonresponse on the outcome - ordinary administrative-panel missingness
qui gen double u = runiform()
qui replace y = . if u < 0.15
drop u
qui count
dis "ROWS_IN_DATA=" r(N)
qui count if !missing(y)
dis "ROWS_WITH_Y=" r(N)
save p13_panel, replace

cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(20) rseed(1234) bootsave(p13_size) replace dots(0)
dis "RC=" _rc
dis "e(N) estimation sample = " e(N)
local eN = e(N)
use p13_size, clear
qui ds w_*
local wv `r(varlist)'
local w1 : word 1 of `wv'
qui su `w1' if _rep==0
dis "OBS_ROW_VALUE = " %21x r(mean) "  = " r(mean)
dis "1/e(N)        = " %21x 1/`eN' "  = " 1/`eN'
use p13_panel, clear
qui count
dis "1/rows_in_data= " %21x 1/r(N) "  = " 1/r(N)
* and the Dirichlet shapes: does alpha_c count the missing-y rows?
use p13_size, clear
qui keep if _rep >= 1
qui ds w_*
mkmat `r(varlist)', matrix(WF)
use p13_panel, clear
qui bysort st (year): gen long n_all = _N
qui bysort st (year): egen long n_y = total(!missing(y))
qui bysort st (year): keep if _n==1
qui sort st
mkmat n_all, matrix(NALL)
mkmat n_y,   matrix(NY)
mata {
  W = st_matrix("WF")
  na = st_matrix("NALL"); ny = st_matrix("NY")
  printf("sum_c n_all_c * w_c  (rep1) = %20.16f   <- 1 iff alpha/divisor use ALL rows\n", W[1,.]*na)
  printf("sum_c n_y_c   * w_c  (rep1) = %20.16f   <- 1 iff they use the estimated rows\n", W[1,.]*ny)
  // recover the implied Dirichlet shape: E[n_c * w_c] = alpha_c / alpha0
  m = mean(W)'
  printf("implied cluster mass share vs n_all share: max|diff| = %g\n", max(abs((m:*na) :- (na/sum(na)))))
  printf("implied cluster mass share vs n_y   share: max|diff| = %g\n", max(abs((m:*na) :- (ny/sum(ny)))))
}
dis "PROBE_DONE"
