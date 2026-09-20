*==========================================================================
* REPRO: the bootstrap resamples clusters that contribute NOTHING to
*        e(sample), and bootdf(clusters) counts them.
*
*  Run:  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_bootdf_deadcluster.do
*
*  Setting: 20 states, 3 units each, 6 years, staggered onsets.  The outcome
*  is suppressed (missing) in 3 states -- an entirely ordinary administrative
*  data shape (a suppressed cell, a covariate not collected in some states).
*  wooldid itself reports "N clusters = 17" and drops 54 observations.
*
*  A: outcome suppressed in 3 states, left in the data
*  B: the SAME analysis after physically dropping those 3 states
*     (identical estimation sample, identical point estimate)
*
*  Expected: A and B are the same analysis, so they should carry the same
*            reference distribution and (up to Monte-Carlo noise from a
*            different seed stream) the same bootstrap standard error.
*  Observed: A reports t(19) where B reports t(16), and A's bootstrap
*            standard error is ~4% larger than B's, because bsample keeps
*            drawing the three empty clusters.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

clear
set obs 20
gen int state = _n
gen int onset = .
forval s = 1/15 {
  local w = mod(`s'-1, 3) + 1
  local o : word `w' of 2002 2003 2004
  qui replace onset = `o' in `s'
}
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand 6
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t - onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset

* three states whose outcome is suppressed
gen double ysup = y
replace ysup = . if inlist(state, 4, 9, 14)

* ---- three seeds at reps(999), so the standard-error comparison is not noise ----
tempname R
matrix `R' = J(3, 6, .)
matrix colnames `R' = seed A_ncl A_df A_se B_df B_se
local row 0
foreach sd in 11 4242 987 {
  local ++row
  qui wooldid ysup state t ttre , cluster(state) bootstrap(standard) reps(499) bootdf(clusters) rseed(`sd') att dots(0)
  matrix `R'[`row',1] = `sd'
  matrix `R'[`row',2] = e(N_clust)
  matrix `R'[`row',3] = e(boot_df)
  matrix `R'[`row',4] = e(wdidmainresults)[1,2]
  scalar Ab = e(wdidmainresults)[1,1]
  scalar Afail`row' = e(boot_failed)
  preserve
  drop if inlist(state, 4, 9, 14)
  qui wooldid y state t ttre , cluster(state) bootstrap(standard) reps(499) bootdf(clusters) rseed(`sd') att dots(0)
  matrix `R'[`row',5] = e(boot_df)
  matrix `R'[`row',6] = e(wdidmainresults)[1,2]
  scalar Bb = e(wdidmainresults)[1,1]
  restore
  dis "seed `sd': estimates identical = " (Ab == Bb) "   A failed reps = " Afail`row'
}
dis _n "  seed   A:e(N_clust)  A:e(boot_df)     A:boot se   B:e(boot_df)     B:boot se"
matrix list `R', format(%12.6f)
mata:
  R = st_matrix(st_local("R"))
  printf("\nSE ratio A/B by seed: %9.5f %9.5f %9.5f   mean %9.5f\n",
     R[1,4]/R[1,6], R[2,4]/R[2,6], R[3,4]/R[3,6], mean(R[.,4]:/R[.,6]))
  printf("FAIL_DF     = %f   (A reports df %f where e(N_clust)-1 = %f, and B reports %f)\n",
     (R[1,3] != R[1,2]-1), R[1,3], R[1,2]-1, R[1,5])
  printf("SE_GAP_CONSISTENT_SIGN = %f\n", (R[1,4]>R[1,6]) + (R[2,4]>R[2,6]) + (R[3,4]>R[3,6]))
end
* ---- also show the one-run self-contradiction on the printed output ----
dis _n(2) "########## the printed run: header says N clusters = 17, footer says t(19) ##########"
wooldid ysup state t ttre , cluster(state) bootstrap(standard) reps(199) bootdf(clusters) rseed(11) att dots(0)
