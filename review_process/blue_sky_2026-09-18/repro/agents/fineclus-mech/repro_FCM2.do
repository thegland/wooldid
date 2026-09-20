*======================================================================
* repro_FCM2.do -- SELF-CONTAINED.
*
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_FCM2.do
*
* CLAIM: wooldid's bootstrap(standard) re-identification leg decides its
* two-regime unit rule on a row set (`smp', ado:8911-8924) that is BROADER
* than the estimation sample, so cluster-variable values on rows the
* regression never sees change the reported bootstrap standard error of
* att_i by about 30% -- silently, with no message and no e() record.
*
* Two datasets, built to be identical wherever the estimator looks:
*   C0 : firm 1's outcome is missing from 2005 on.        -> nesting regime
*   C  : the same, and firm 1's industry code changes on
*        exactly those missing-outcome rows.              -> keep-i regime
* Both give N = 764, the same point estimate to every digit, and the same
* analytic clustered standard error.  Only the bootstrap SE moves.
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

prog def mkpanel
  syntax , [phantom(int 0)]
  clear
  set seed 771
  set sortseed 771
  set obs 96                                   /* 24 industries x 4 firms */
  gen long firm = _n
  gen long ind  = ceil(_n/4)
  gen int onset = .
  replace onset = 2004 if mod(firm,4)==1
  replace onset = 2006 if mod(firm,4)==2
  replace onset = 2008 if mod(firm,4)==3
  gen double indshk = .
  gen double indtau = .
  bysort ind (firm): replace indshk = rnormal()         if _n==1
  bysort ind (firm): replace indtau = 1 + rnormal()*0.8 if _n==1
  bysort ind (firm): replace indshk = indshk[1]
  bysort ind (firm): replace indtau = indtau[1]
  gen double ai = rnormal()
  expand 8
  bysort firm (ind): gen int year = 2001 + _n - 1
  gen double e = rnormal()*0.3
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = ai + indshk + 0.05*(year-2001) + indtau*D + e
  * firm 1 stops reporting its outcome in 2005 (BOTH datasets)
  replace y = . if firm == 1 & year >= 2005
  * ... and in the phantom dataset its industry code also changes, on
  *     exactly the rows whose outcome is missing
  if `phantom' == 1 replace ind = 2 if firm == 1 & year >= 2005
  sort firm year
end

prog def getatti
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  local j 0
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if "`this'" == "main_att_i" local j `q'
  }
  c_local ATTI   = `b'[1,`j']
  c_local ATTISE = sqrt(`V'[`j',`j'])
end

*--- show that the two datasets differ ONLY outside the estimation sample
mkpanel , phantom(1)
qui wooldid y firm year onset, att_i cluster(ind)
gen byte es = e(sample)
count if es == 1 & (missing(y) | missing(ind))
dis "FCM2 rows in e(sample) with missing y or missing ind = " r(N) "  (must be 0)"
count if es == 0
dis "FCM2 rows OUTSIDE e(sample) = " r(N)
list firm year ind y es if es == 0, noobs
* the regime predicate, evaluated the way ado:8989-8991 evaluates it
tempvar smp mn mx
qui gen byte `smp' = 1
qui replace `smp' = 0 if missing(ind)
qui gegen double `mn' = min(ind) if `smp', by(firm)
qui gegen double `mx' = max(ind) if `smp', by(firm)
cap assert `mn' == `mx' if `smp'
dis "FCM2 NEST as wooldid computes it (on smp)   = " (_rc == 0)
drop `mn' `mx'
qui gegen double `mn' = min(ind) if es, by(firm)
qui gegen double `mx' = max(ind) if es, by(firm)
cap assert `mn' == `mx' if es
dis "FCM2 NEST on the ESTIMATION SAMPLE          = " (_rc == 0)

*--- the two runs, same seed, same draws
foreach P in 0 1 {
  mkpanel , phantom(`P')
  qui wooldid y firm year onset, att_i cluster(ind)
  getatti
  dis "FCM2 phantom=`P' ANALYTIC att_i=" %18.14f `ATTI' " se=" %18.14f `ATTISE' " N=" e(N) " Nclust=" e(N_clust)
  cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(200) rseed(4242) dots(0)
  dis "FCM2 phantom=`P' BOOT RC=" _rc
  getatti
  dis "FCM2 phantom=`P' BOOTSTRAP att_i=" %18.14f `ATTI' " se=" %18.14f `ATTISE' " N=" e(N) " Nclust=" e(N_clust)
  * FCM-3: is the regime recorded anywhere in e()?
  qui ereturn list
  local emac : e(macros)
  local hit 0
  foreach m of local emac {
    if strpos(lower("`m'"), "nest") | strpos(lower("`m'"), "reid") local hit 1
  }
  dis "FCM2 phantom=`P' e() macro naming the regime (nest/reid) present = `hit'  (0 = no record at all)"
}
dis "FCM2 DONE  -- compare the two BOOTSTRAP se's; phantom=1 is ~0.71 of phantom=0"
