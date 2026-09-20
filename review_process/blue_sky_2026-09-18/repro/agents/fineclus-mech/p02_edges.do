clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P02 -- edge geometries of the two-regime bootstrap unit rule
*   (1)  multi-variable cluster()  (two-way clustering: cluster(ind year))
*   (2)  cluster() EXACTLY EQUAL to i
*   (3)  cluster() STRICTLY FINER than i
*   (4)  one cluster per observation
*   (5)  constant cluster variable (one cluster)
*   (6)  missing cluster ids
*   (7)  i with missing values
*   (8)  a user `if' that removes the straddling rows
*   (9)  phantom straddle through a MISSING CONTROL
*   (10) an entire cluster out of e(sample)
*----------------------------------------------------------------------
prog def mkpanel
  syntax , [misscl(int 0)]
  clear
  set seed 771
  set sortseed 771
  local NI 24
  local NF 4
  set obs `=`NI'*`NF''
  gen long firm = _n
  gen long ind  = ceil(_n/`NF')
  gen int onset = .
  replace onset = 2004 if mod(firm,4)==1
  replace onset = 2006 if mod(firm,4)==2
  replace onset = 2008 if mod(firm,4)==3
  gen double indshk = .
  gen double indtau = .
  bysort ind (firm): replace indshk = rnormal()          if _n==1
  bysort ind (firm): replace indtau = 1 + rnormal()*0.8   if _n==1
  bysort ind (firm): replace indshk = indshk[1]
  bysort ind (firm): replace indtau = indtau[1]
  gen double ai = rnormal()
  expand 8
  bysort firm (ind): gen int year = 2001 + _n - 1
  gen double e = rnormal()*0.3
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = ai + indshk + 0.05*(year-2001) + indtau*D + e
  * a coarser cohort variable: 12 cohort-groups, each holding 2 industries
  gen long cgrp = ceil(ind/2)
  * one cluster per observation
  gen long rowid = _n
  gen byte one = 1
  if `misscl' == 1 replace ind = . if firm == 3 & year >= 2005
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
  if `j' == 0 {
    c_local ATTI = .
    c_local ATTISE = .
  }
  else {
    c_local ATTI = `b'[1,`j']
    c_local ATTISE = sqrt(`V'[`j',`j'])
  }
end

mkpanel
save "p02_base.dta", replace

*--- (1) two-way cluster, plain -----------------------------------------
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att_i cluster(ind year)
dis "T1_PLAIN2WAY_RC=" _rc
getatti
dis "T1 ATTI=" `ATTI' " SE=" `ATTISE'

*--- (1b) two-way cluster, bootstrap COUNT path (att only, no reid) ------
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att cluster(ind year) bootstrap(standard) reps(20) rseed(7) dots(0)
dis "T1B_COUNT2WAY_RC=" _rc

*--- (1c) two-way cluster, bootstrap RE-ID path (att_i) ------------------
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att_i cluster(ind year) bootstrap(standard) reps(20) rseed(7) dots(0)
dis "T1C_REID2WAY_RC=" _rc
getatti
dis "T1C ATTI=" `ATTI' " SE=" `ATTISE'

*--- (2) cluster() == i --------------------------------------------------
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att_i cluster(firm) bootstrap(standard) reps(100) rseed(7) dots(0)
dis "T2_CLUSEQI_RC=" _rc
getatti
dis "T2 ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE' " NCLUST=" e(N_clust)

*--- (3) cluster() strictly FINER than i (i = cohort group, cluster = firm)
use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att_i cluster(firm) bootstrap(standard) reps(100) rseed(7) dots(0)
dis "T3_FINER_RC=" _rc
getatti
dis "T3 ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE' " NCLUST=" e(N_clust)

use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att_i cluster(ind) bootstrap(standard) reps(100) rseed(7) dots(0)
dis "T3B_COARSER_RC=" _rc
getatti
dis "T3B ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE' " NCLUST=" e(N_clust)

*--- (4) one cluster per observation ------------------------------------
use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att_i cluster(rowid) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "T4_ROWCLUS_RC=" _rc
getatti
dis "T4 ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE' " NCLUST=" e(N_clust)

*--- (5) constant cluster variable ---------------------------------------
use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att_i cluster(one) bootstrap(standard) reps(20) rseed(7) dots(0)
dis "T5_ONECLUS_RC=" _rc

*--- (6) missing cluster ids ---------------------------------------------
mkpanel , misscl(1)
save "p02_misscl.dta", replace
cap noisily wooldid y firm year onset, att_i cluster(ind)
dis "T6_PLAIN_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "T6_BOOT_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
getatti
dis "T6 ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE'

*--- (7) i with MISSING values on a few rows -----------------------------
mkpanel
replace firm = . if _n <= 3
save "p02_missi.dta", replace
cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(40) rseed(7) dots(0)
dis "T7_MISSI_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)

*--- (8) a user `if' that removes the straddling rows --------------------
mkpanel
replace ind = 2 if firm == 3 & year >= 2005
save "p02_ifstrad.dta", replace
cap noisily wooldid y firm year onset if year < 2005, att_i cluster(ind) bootstrap(standard) reps(40) rseed(7) dots(0)
dis "T8_IFSTRAD_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)

*--- (9) phantom straddle through a MISSING CONTROL ----------------------
mkpanel
gen double x = rnormal()
replace ind = 2 if firm == 3 & year >= 2005
replace x = .   if firm == 3 & year >= 2005
save "p02_missctl.dta", replace
cap noisily wooldid y firm year onset, att_i controls(x) cluster(ind) bootstrap(standard) reps(100) rseed(7) dots(0)
dis "T9_MISSCTL_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
getatti
dis "T9 ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE'

mkpanel
gen double x = rnormal()
replace x = . if firm == 3 & year >= 2005
save "p02_missctl0.dta", replace
cap noisily wooldid y firm year onset, att_i controls(x) cluster(ind) bootstrap(standard) reps(100) rseed(7) dots(0)
dis "T9B_MISSCTL0_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
getatti
dis "T9B ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE'

*--- (10) an entire cluster out of e(sample) (all y missing) -------------
mkpanel
replace y = . if ind == 5
save "p02_deadclus.dta", replace
cap noisily wooldid y firm year onset, att_i cluster(ind)
dis "T10_PLAIN_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(100) rseed(7) dots(0) bootdf(clusters)
dis "T10_BOOT_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust) " DFR=" e(df_r)
getatti
dis "T10 ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE'

dis "ALLDONE"
