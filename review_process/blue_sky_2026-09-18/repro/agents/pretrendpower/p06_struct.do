* ================================================================
* p06_struct.do -- deterministic structural probes of pretrendpower()
* Each block prints a TAG so the log can be grepped.  Everything that
* may legitimately fail is wrapped in capture noisily.
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

capture program drop mkpanel
program define mkpanel
  syntax , [delta(real 0) sd(real 1) anticip(int 0)]
  clear
  set obs 80
  gen long unit = _n
  gen double a_i = 0
  gen int onset = .
  replace onset = 5 if mod(unit,5)==1
  replace onset = 7 if mod(unit,5)==2
  replace onset = 9 if mod(unit,5)==3
  expand 10
  bysort unit: gen int t = _n
  gen double f_t = 0.1*t + 0.3*sin(t)
  gen byte D = (onset<. & t>=onset)
  gen double relyr = t - onset
  gen double trendterm = cond(onset<., relyr + 1, 0)
  gen double y = a_i + f_t + 1.0*D + `delta'*trendterm + `sd'*rnormal()
end

* ---------------------------------------------------------------
* (A0) does a planted linear pre-trend of slope delta really put
*      delta*m into the pre-period event-study coefficients?
*      Near-zero noise makes this an exact algebraic check of the
*      premise the whole power claim rests on.
* ---------------------------------------------------------------
dis "TAG_A0"
mkpanel , delta(0.25) sd(1e-8)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "A0_RC=" _rc
mat list e(b_es_att)
dis "A0 slope (should be .25) = " e(roth_slope_att)
dis "A0 hbar = " e(roth_hbar_att)

* ---------------------------------------------------------------
* (A1) esrelativeto(-2): the pre relyrs straddle the reference, so
*      m must carry BOTH signs.  Planted trend again; b_pre must be
*      delta*(relyr - ref) coefficient by coefficient.
* ---------------------------------------------------------------
dis "TAG_A1"
mkpanel , delta(0.25) sd(1e-8)
cap noi wooldid y unit t onset, esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(2) pretrendpower(80)
dis "A1_RC=" _rc
mat list e(b_es_att)
dis "A1 slope (should be .25) = " e(roth_slope_att)

* ---------------------------------------------------------------
* (A2) anticipation(1): is hbar measured on the SHIFTED onset?
*      hand-compute both candidates and see which e(roth_hbar_att) is.
* ---------------------------------------------------------------
dis "TAG_A2"
mkpanel , delta(0) sd(0.3)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) anticipation(1) pretrendpower(80)
dis "A2_RC=" _rc
dis "A2 e(roth_hbar_att) = " e(roth_hbar_att)
gen byte insamp = e(sample)
gen double rly_shift = t - (onset - 1)
gen double rly_raw   = t - onset
qui sum rly_shift if insamp & onset<. & rly_shift>=0
dis "A2 hbar on SHIFTED onset = " r(mean)
qui sum rly_raw if insamp & onset<. & rly_raw>=0
dis "A2 hbar on RAW onset     = " r(mean)

* ---------------------------------------------------------------
* (A3) the help's own worked example 7 (sthlp:1546), literally:
*      variables named y i t ttre and cluster(i).  Does the run print
*      its own "should not be read as a calibrated test" caution?
* ---------------------------------------------------------------
dis "TAG_A3"
mkpanel , delta(0) sd(0.3)
rename unit i
rename onset ttre
cap noi wooldid y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)
dis "A3_RC=" _rc

* ---------------------------------------------------------------
* (A4) jointtests vs roth_prep: is "the joint pre-trend test above"
*      the same test jointtests reports?
* ---------------------------------------------------------------
dis "TAG_A4"
mkpanel , delta(0.05) sd(0.3)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) jointtests pretrendpower(80)
dis "A4_RC=" _rc
dis "A4 roth_prep_att = " e(roth_prep_att)
dis "A4 roth_W_att    = " e(roth_W_att)
ereturn list

* ---------------------------------------------------------------
* (A5) mdt must not depend on the estimated coefficients (sthlp:791).
*      Same design, same seed, two very different planted trends.
* ---------------------------------------------------------------
dis "TAG_A5"
mkpanel , delta(0) sd(0.3)
cap noi qui wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
mat R0 = e(roth_power_att)
dis "A5 mdt_wald(delta=0)   = " %21x R0[1,2] "   mdt_nis = " %21x R0[1,3]
mkpanel , delta(3.0) sd(0.3)
cap noi qui wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
mat R1 = e(roth_power_att)
dis "A5 mdt_wald(delta=3.0) = " %21x R1[1,2] "   mdt_nis = " %21x R1[1,3]

* ---------------------------------------------------------------
* (A6) a requested pre horizon that no cohort can supply
* ---------------------------------------------------------------
dis "TAG_A6"
mkpanel , delta(0) sd(0.3)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(9) espostlength(2) pretrendpower(80)
dis "A6_RC=" _rc

* ---------------------------------------------------------------
* (A7) bootstrap(): which df does the Roth header carry, and how far
*      does the bootstrap covariance move mdt_wald?
* ---------------------------------------------------------------
dis "TAG_A7"
mkpanel , delta(0) sd(0.3)
cap noi qui wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
mat RA = e(roth_power_att)
dis "A7 analytic mdt_wald = " RA[1,2] "  df_r = " e(df_r) "  slope_se = " e(roth_slope_se_att)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80) bootstrap(standard) reps(99) rseed(4242)
dis "A7_RC=" _rc
mat RB = e(roth_power_att)
dis "A7 bootstrap mdt_wald = " RB[1,2] "  df_r = " e(df_r) "  boot_df = " e(boot_df) "  slope_se = " e(roth_slope_se_att)
dis "A7 boot_method = " e(boot_method)


* ---------------------------------------------------------------
* (A8) BIENNIAL panel: the relative-period axis has gaps.  m must be
*      built in t units.  Planted trend + near-zero noise again.
* ---------------------------------------------------------------
dis "TAG_A8"
clear
set obs 60
gen long unit = _n
gen int onset = .
replace onset = 2008 if mod(unit,3)==1
replace onset = 2012 if mod(unit,3)==2
expand 10
bysort unit: gen int k = _n
gen int t = 2000 + 2*(k-1)
gen double f_t = 0.05*k
gen byte D = (onset<. & t>=onset)
gen double relyr = t - onset
gen double y = f_t + 1.0*D + 0.25*cond(onset<., relyr + 2, 0) + 1e-8*rnormal()
cap noi wooldid y unit t onset, esfixedbaseperiod esrelativeto(-2) esprelength(6) espostlength(4) pretrendpower(80)
dis "A8_RC=" _rc
mat list e(b_es_att)
dis "A8 slope (should be .25) = " e(roth_slope_att)
dis "A8 hbar = " e(roth_hbar_att)
dis "A8 npre = " e(roth_npre_att)

* ---------------------------------------------------------------
* (A9) subgroups() + pretrendpower(): which families get diagnostics?
* ---------------------------------------------------------------
dis "TAG_A9"
mkpanel , delta(0) sd(0.3)
gen byte sg = mod(unit,2)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) subgroups(sg) pretrendpower(80)
dis "A9_RC=" _rc
ereturn list

* ---------------------------------------------------------------
* (A10) poisson: does the Roth block fall back to the chi-squared
*       shortcut its own comment calls "optimistic by ~25 percent at
*       10 clusters"?
* ---------------------------------------------------------------
dis "TAG_A10"
clear
set obs 80
gen long unit = _n
gen int onset = .
replace onset = 5 if mod(unit,5)==1
replace onset = 7 if mod(unit,5)==2
replace onset = 9 if mod(unit,5)==3
expand 10
bysort unit: gen int t = _n
gen byte D = (onset<. & t>=onset)
gen double mu = exp(1.5 + 0.05*t + 0.3*D)
gen double y = rpoisson(mu)
cap noi wooldid y unit t onset, poisson esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "A10_RC=" _rc
dis "A10 df_r = " e(df_r)
dis "A10 roth_keff = " e(roth_keff_att)
dis "A10 roth_prep = " e(roth_prep_att)

* ---------------------------------------------------------------
* (A11) mvnormalcvqp vech ordering: does Mata's COLUMN-major vech()
*       match what mvnormalcvqp expects?  Decisive Monte-Carlo check
*       on a matrix whose row-major and column-major vech differ.
* ---------------------------------------------------------------
dis "TAG_A11"
mata:
  S = (1, .30, .70 \ .30, 4, .20 \ .70, .20, 9)
  sd = sqrt(diagonal(S))'
  c  = 2
  L = -c:*sd
  U =  c:*sd
  pin = mvnormalcvqp(L, U, J(1,3,0), vech(S)', 1000)
  printf("A11 mvnormalcvqp with vech(S)' (column-major)  = %12.8f\n", pin)
  vr = (S[1,1], S[2,1], S[2,2], S[3,1], S[3,2], S[3,3])
  pin2 = mvnormalcvqp(L, U, J(1,3,0), vr, 1000)
  printf("A11 mvnormalcvqp with ROW-major lower triangle = %12.8f\n", pin2)
  rseed(77)
  n = 400000
  X = rnormal(n, 3, 0, 1) * cholesky(S)'
  inb = (abs(X[,1]):<U[1]) :* (abs(X[,2]):<U[2]) :* (abs(X[,3]):<U[3])
  printf("A11 brute-force Monte Carlo P(all inside)      = %12.8f\n", mean(inb))
end



* ---------------------------------------------------------------
* (A12) bootstrap() + cluster(): which degrees of freedom does the
*       Roth block's F test and mdt_wald use?  ado:449 passes
*       e(df_r) (the ANALYTIC residual df) into wooldid_ptp even
*       though the covariance it is about to use is the bootstrap
*       one and every other number in the run is on t(G-1).
* ---------------------------------------------------------------
dis "TAG_A12"
clear
set obs 16
gen long st = _n
gen byte coh = 0
replace coh = 1 if st <= 4
replace coh = 2 if st > 4 & st <= 8
replace coh = 3 if st > 8 & st <= 12
gen int ttre = .
replace ttre = 5 if coh == 1
replace ttre = 7 if coh == 2
replace ttre = 9 if coh == 3
expand 10
bysort st: gen int t = _n
sort st t
gen double f_t = 0.1*t + 0.3*sin(t)
gen byte D = (ttre<. & t>=ttre)
gen double relyr = t - ttre
gen double y = f_t + 1.0*D + rnormal()
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "A12a_RC=" _rc
dis "A12a analytic: df_r = " e(df_r) "  keff = " e(roth_keff_att) "  prep = " e(roth_prep_att)
mat A12A = e(roth_power_att)
dis "A12a analytic mdt_wald(80) = " A12A[1,2] "  nissize = " e(roth_nissize_att)
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80) bootstrap(standard) reps(199) rseed(9191)
dis "A12b_RC=" _rc
dis "A12b bootstrap: df_r = " e(df_r) "  boot_df = " e(boot_df) "  keff = " e(roth_keff_att) "  prep = " e(roth_prep_att)
mat A12B = e(roth_power_att)
dis "A12b bootstrap mdt_wald(80) = " A12B[1,2] "  nissize = " e(roth_nissize_att)
dis "A12b boot_method = " e(boot_method) "  boot_stat = " e(boot_stat)

* ---------------------------------------------------------------
* (A13) unconditionalse and customvce(): what df reaches the Roth block?
* ---------------------------------------------------------------
dis "TAG_A13"
cap noi wooldid y coh t ttre, fe(st) customvce(cluster st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "A13_RC=" _rc
dis "A13 df_r = " e(df_r) "  prep = " e(roth_prep_att)


* ---------------------------------------------------------------
* (A14) att_it request: the shared curve stands in for att_i/att_itime
* ---------------------------------------------------------------
dis "TAG_A14"
mkpanel , delta(0) sd(0.3)
cap noi wooldid y unit t onset, att att_it att_i esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "A14_RC=" _rc
dis "A14 roth_hbar_att = " e(roth_hbar_att) "  roth_hbar_att_it = " e(roth_hbar_att_it)

* ---------------------------------------------------------------
* (A15) k = 1.  With ONE pre-period coefficient the joint Wald test
*       and the any-individually-significant test are THE SAME TEST
*       (|b/se| > tcrit).  So mdt_wald and mdt_nis must coincide and
*       roth_nissize must be 0.05.  Any gap is exactly the cost of
*       evaluating a NORMAL probability at a t critical value.
*       Run it at df_r = 279 (large) and at df_r = 15 (G=16 clusters).
* ---------------------------------------------------------------
dis "TAG_A15"
mkpanel , delta(0) sd(0.3)
cap noi wooldid y unit t onset, esfixedbaseperiod esprelength(1) espostlength(2) pretrendpower(80)
dis "A15a_RC=" _rc
dis "A15a df_r = " e(df_r) "  npre = " e(roth_npre_att) "  nissize = " e(roth_nissize_att)
mat K1A = e(roth_power_att)
dis "A15a mdt_wald = " K1A[1,2] "   mdt_nis = " K1A[1,3] "   ratio = " K1A[1,3]/K1A[1,2]

clear
set obs 16
gen long st = _n
gen byte coh = 0
replace coh = 1 if st <= 4
replace coh = 2 if st > 4 & st <= 8
replace coh = 3 if st > 8 & st <= 12
gen int ttre = .
replace ttre = 5 if coh == 1
replace ttre = 7 if coh == 2
replace ttre = 9 if coh == 3
expand 10
bysort st: gen int t = _n
sort st t
gen double f_t = 0.1*t + 0.3*sin(t)
gen byte D = (ttre<. & t>=ttre)
gen double y = f_t + 1.0*D + rnormal()
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(1) espostlength(2) pretrendpower(80)
dis "A15b_RC=" _rc
dis "A15b df_r = " e(df_r) "  npre = " e(roth_npre_att) "  nissize = " e(roth_nissize_att)
mat K1B = e(roth_power_att)
dis "A15b mdt_wald = " K1B[1,2] "   mdt_nis = " K1B[1,3] "   ratio = " K1B[1,3]/K1B[1,2]
dis "A15b invttail(15,.025) = " invttail(15,.025) "   invnormal(.975) = " invnormal(.975)
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "A15c_RC=" _rc
dis "A15c k=3 small-G: df_r = " e(df_r) "  nissize = " e(roth_nissize_att)
mat K3B = e(roth_power_att)
dis "A15c mdt_wald = " K3B[1,2] "   mdt_nis = " K3B[1,3]

dis "TAG_END"
