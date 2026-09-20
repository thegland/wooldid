*======================================================================
* p_lens1.do -- FCM-2, third lens (documented-behaviour). SELF-CONTAINED.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p_lens1.do
*
* THREE ARMS, all with the SAME estimation sample:
*   A : firm 1's outcome missing from 2005 on; its industry never changes.
*   B : the same, and firm 1's industry code changes on exactly the
*       missing-outcome rows (the "phantom").
*   C : arm B with those four rows PHYSICALLY DELETED from the dataset.
* A and C both have a cluster that nests i on every row that exists in
* their data; B differs from C only by four rows the estimator drops.
* If se_boot(A) == se_boot(C) != se_boot(B), the reported bootstrap SE is a
* function of rows outside the estimation sample.
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

prog def mkpanel
  syntax , [phantom(int 0) dropem(int 0)]
  clear
  set seed 771
  set sortseed 771
  set obs 96
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
  replace y = . if firm == 1 & year >= 2005
  if `phantom' == 1 replace ind = 2 if firm == 1 & year >= 2005
  if `dropem'  == 1 drop if firm == 1 & year >= 2005
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

*---------------------------------------------------------------
* PART 1 - regime on smp vs on e(sample), all three arms
*---------------------------------------------------------------
foreach arm in A B C {
  if "`arm'" == "A" mkpanel , phantom(0) dropem(0)
  if "`arm'" == "B" mkpanel , phantom(1) dropem(0)
  if "`arm'" == "C" mkpanel , phantom(1) dropem(1)
  qui count
  local NROW = r(N)
  cap noisily wooldid y firm year onset, att_i cluster(ind)
  dis "L1_`arm'_RC=" _rc
  getatti
  qui gen byte es = e(sample)
  qui count if es
  local NEST_ = r(N)
  * regime as the code computes it: on smp = all rows (no if, ind never missing)
  tempvar cmn cmx
  qui gegen double `cmn' = min(ind), by(firm)
  qui gegen double `cmx' = max(ind), by(firm)
  cap assert `cmn' == `cmx'
  local NSMP = (_rc == 0)
  drop `cmn' `cmx'
  * regime on the ESTIMATION sample
  tempvar cmn2 cmx2
  qui gegen double `cmn2' = min(ind) if es, by(firm)
  qui gegen double `cmx2' = max(ind) if es, by(firm)
  cap assert `cmn2' == `cmx2' if es
  local NEST = (_rc == 0)
  drop `cmn2' `cmx2'
  qui count if es == 0
  local NOUT = r(N)
  dis "L1_`arm'  Nrows=`NROW'  Nest=`NEST_'  rows_outside_esample=`NOUT'"
  dis "L1_`arm'  NEST_on_smp=`NSMP'   NEST_on_esample=`NEST'"
  dis "L1_`arm'  att_i  = " %21x `ATTI'
  dis "L1_`arm'  se_ana = " %21x `ATTISE'
  dis "L1_`arm'  Nclust = " e(N_clust)
}

*---------------------------------------------------------------
* PART 2 - the bootstrap, same reps and same rseed in all three arms
*---------------------------------------------------------------
foreach arm in A B C {
  if "`arm'" == "A" mkpanel , phantom(0) dropem(0)
  if "`arm'" == "B" mkpanel , phantom(1) dropem(0)
  if "`arm'" == "C" mkpanel , phantom(1) dropem(1)
  dis "===== BOOTSTRAP ARM `arm' BEGIN ====="
  cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(100) rseed(90210) dots(0)
  dis "L2_`arm'_RC=" _rc
  dis "===== BOOTSTRAP ARM `arm' END ====="
  getatti
  dis "L2_`arm'  att_i   = " %21x `ATTI'
  dis "L2_`arm'  se_boot = " %21x `ATTISE'
  dis "L2_`arm'  se_boot_dec = " `ATTISE'
  * observability: is the regime recorded anywhere in e()?
  local sc : e(scalars)
  local mc : e(macros)
  local hit 0
  foreach nm of local sc {
    if strpos("`nm'", "nest") | strpos("`nm'", "regime") | strpos("`nm'", "reid") | strpos("`nm'", "rebuil") local hit 1
  }
  foreach nm of local mc {
    if strpos("`nm'", "nest") | strpos("`nm'", "regime") | strpos("`nm'", "reid") | strpos("`nm'", "rebuil") local hit 1
  }
  dis "L2_`arm'  eREGIMEMARKER=`hit'"
  local nsc : word count `sc'
  local nmc : word count `mc'
  dis "L2_`arm'  e_nscalars=`nsc' e_nmacros=`nmc'"
  file open fh using "enames_`arm'.txt", write replace
  file write fh "`sc'" _n "`mc'" _n
  file close fh
}
dis "LENS1_DONE"
