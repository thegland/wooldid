clear all
set more off
set linesize 250
set graphics off
qui do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"

*---------------- DGP A (mc_bootse design A) -----------------
cap prog drop dgpA
prog def dgpA
  args seed G U T tau
  clear
  qui set obs `=`G'*`U''
  gen long unit = _n
  gen long cl   = ceil(unit/`U')
  gen byte uw   = unit - (cl-1)*`U'
  gen int ttre = .
  qui replace ttre = 3 if mod(cl,3)==1
  qui replace ttre = 5 if mod(cl,3)==2
  qui expand `T'
  bysort unit: gen int t = _n
  set seed `seed'
  qui by unit: gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.7)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  tempvar clt
  qui egen long `clt' = group(cl t)
  sort `clt'
  gen double s = 0
  qui by `clt': replace s = rnormal(0,0.5) if _n==1
  qui by `clt': replace s = s[1]
  sort unit t
  gen double y = `tau'*(t>=ttre & ttre<.) + 0.05*t + mu + eps + s
end

*---------------- DGP B (unequal clusters) -----------------
* 12 states, state s has s units, T=8, odd states treated at t=5, tau_s = 0.1*n_s
cap prog drop dgpB
prog def dgpB
  args seed T mult
  clear
  qui set obs 12
  gen long state = _n
  qui expand state
  bysort state: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 5 if mod(state,2)==1
  gen double tau = `mult'*0.1*state if mod(state,2)==1
  qui replace tau = 0 if missing(tau)
  qui expand `T'
  bysort unit: gen int t = _n
  set seed `seed'
  qui by unit: gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.7)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  tempvar clt
  qui egen long `clt' = group(state t)
  sort `clt'
  gen double s = 0
  qui by `clt': replace s = rnormal(0,0.5) if _n==1
  qui by `clt': replace s = s[1]
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.05*t + mu + eps + s
end

*===================== TRUTH CHECK: noiseless design B =====================
di "SMK ===== TRUTH B (noiseless) ====="
dgpB 1 8 1
qui replace y = tau*(t>=ttre & ttre<.)
qui wooldid y unit t ttre, att att_it cluster(state) safety(off)
matrix TB = e(wdidmainresults)
di "SMK rownames wdidmainresults: `: rownames TB'"
di "SMK colnames wdidmainresults: `: colnames TB'"
matrix list TB
di "SMK e(b) colnames: `: colnames e(b)'"
matrix list e(b)
di "SMK analytic exact truth = " %18.14f (0.1*286/36)

*===================== structure: design A, bootstrap standard + ES + jt ==========
di "SMK ===== ARM std + ES + jointtests, design A ====="
dgpA 900001 24 2 6 2
timer clear 1
timer on 1
cap noi wooldid y unit t ttre, att cluster(cl) esprelength(2) espostlength(1) jointtests ///
    bootstrap(standard) reps(199) rseed(4242) dots(0) safety(off) level(90)
di "SMK rc=" _rc
timer off 1
qui timer list 1
di "SMK SECS_stdES = " r(t1)
matrix M = e(wdidmainresults)
di "SMK rownames: `: rownames M'"
di "SMK colnames: `: colnames M'"
matrix list M
matrix CP = e(boot_ci_percentile)
di "SMK ci_pct colnames: `: colnames CP'"
matrix list CP
matrix PB = e(boot_p_percentile)
matrix list PB
matrix CB = e(boot_ci_bc)
matrix list CB
matrix PBC = e(boot_p_bc)
matrix list PBC
matrix WM = e(boot_wald_mainresults)
di "SMK wald rownames: `: rownames WM'"
matrix list WM
di "SMK joint_horiz = " e(joint_horiz) " family=" e(boot_jt_family) " boot_N_V=" e(boot_N_V) " level=" e(boot_level)

*===================== plain A without ES: does main_att move? ==========
di "SMK ===== plain A with and without ES ====="
qui wooldid y unit t ttre, att cluster(cl) safety(off) level(90)
matrix P1 = e(wdidmainresults)
di "SMK plain noES b=" %18.12g el(P1,1,1) " se=" %18.12g el(P1,1,2) " rows=" rowsof(P1)
qui wooldid y unit t ttre, att cluster(cl) esprelength(2) espostlength(1) jointtests safety(off) level(90)
matrix P2 = e(wdidmainresults)
di "SMK plain ES rownames: `: rownames P2'"
di "SMK plain ES  b=" %18.12g el(P2,1,1) " se=" %18.12g el(P2,1,2) " rows=" rowsof(P2)

*===================== timings: design A other arms ==========
di "SMK ===== timings design A ====="
dgpA 900001 24 2 6 2
timer clear 2
timer on 2
cap noi qui wooldid y unit t ttre, att cluster(cl) bootstrap(bayes) reps(199) rseed(51) dots(0) safety(off) level(90)
di "SMK bayes rc=" _rc
timer off 2
qui timer list 2
di "SMK SECS_bayesA = " r(t2)

timer clear 3
timer on 3
cap noi wooldid y unit t ttre, att cluster(cl) esprelength(2) espostlength(1) jointtests ///
   bootstrap(standard) boottstat reps(199) rseed(61) dots(0) safety(off) level(90)
di "SMK bts rc=" _rc
timer off 3
qui timer list 3
di "SMK SECS_btsA = " r(t3)
matrix BT = e(wdidmainresults)
di "SMK bts rownames: `: rownames BT'"
matrix list BT
di "SMK bts joint_horiz = " e(joint_horiz) " family=" e(boot_jt_family)

*===================== design B arms ==========
di "SMK ===== timings design B ====="
dgpB 910001 8 1
qui count
di "SMK B nobs=" r(N)
timer clear 4
timer on 4
cap noi wooldid y unit t ttre, att cluster(state) esprelength(2) espostlength(1) jointtests ///
   bootstrap(standard) reps(199) rseed(71) dots(0) safety(off) level(90)
di "SMK B std rc=" _rc
timer off 4
qui timer list 4
di "SMK SECS_stdB = " r(t4)
matrix MB = e(wdidmainresults)
di "SMK B rownames: `: rownames MB'"
matrix list MB
matrix list e(boot_ci_percentile)
di "SMK B joint_horiz=" e(joint_horiz) " fam=" e(boot_jt_family)

timer clear 5
timer on 5
cap noi qui wooldid y unit t ttre, att cluster(state) bootstrap(standard) bootdisplay(bca) reps(199) rseed(81) dots(0) safety(off) level(90)
di "SMK B bca rc=" _rc
timer off 5
qui timer list 5
di "SMK SECS_bcaB = " r(t5)
cap matrix list e(boot_ci_bca)
cap matrix list e(boot_p_bca)
di "SMK B jk_completed=" e(boot_jk_completed)

timer clear 6
timer on 6
cap noi qui wooldid y unit t ttre, att_it cluster(state) bootstrap(standard) reps(199) rseed(91) dots(0) safety(off) level(90)
di "SMK B attit rc=" _rc
timer off 6
qui timer list 6
di "SMK SECS_attitB = " r(t6)
cap matrix list e(wdidmainresults)
cap matrix list e(boot_ci_percentile)
cap di "SMK attit b colnames: `: colnames e(b)'"

timer clear 7
timer on 7
cap noi qui wooldid y unit t ttre, att cluster(state) esprelength(2) espostlength(1) jointtests bootstrap(standard) boottstat reps(199) rseed(101) dots(0) safety(off) level(90)
di "SMK B bts rc=" _rc
timer off 7
qui timer list 7
di "SMK SECS_btsB = " r(t7)
cap di "SMK B bts joint_horiz=" e(joint_horiz) " fam=" e(boot_jt_family)

timer clear 8
timer on 8
cap noi qui wooldid y unit t ttre, att cluster(state) bootstrap(bayes) reps(199) rseed(111) dots(0) safety(off) level(90)
di "SMK B bayes rc=" _rc
timer off 8
qui timer list 8
di "SMK SECS_bayesB = " r(t8)

timer clear 9
timer on 9
cap noi qui wooldid y unit t ttre, att cluster(state) safety(off) level(90)
timer off 9
qui timer list 9
di "SMK SECS_plainB = " r(t9)
di "SMK DONE"
