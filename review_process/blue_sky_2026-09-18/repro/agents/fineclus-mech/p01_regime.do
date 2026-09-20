clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P01  -- does the two-regime bootstrap unit rule flip on rows that are
*         NOT in the estimation sample, and how much does the SE move?
*
* Design: i = firm, cluster() = industry (STRICTLY COARSER; firms nest
*         inside industries).  Treatment effects vary BY INDUSTRY, so the
*         composition variance of att_i across industry draws is real.
*----------------------------------------------------------------------
prog def mkpanel
  syntax , [straddle(int 0) missy(int 0)]
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
  * industry-level effect heterogeneity + industry shock
  gen double indshk = .
  gen double indtau = .
  bysort ind (firm): replace indshk = rnormal()      if _n==1
  bysort ind (firm): replace indtau = 1 + rnormal()*0.8 if _n==1
  bysort ind (firm): replace indshk = indshk[1]
  bysort ind (firm): replace indtau = indtau[1]
  gen double ai = rnormal()
  expand 8
  bysort firm (ind): gen int year = 2001 + _n - 1
  gen double e = rnormal()*0.3
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = ai + indshk + 0.05*(year-2001) + indtau*D + e
  * ---- the perturbation: firm 1 is reclassified into industry 2 from 2005 on
  if `straddle' == 1 replace ind = 2 if firm == 1 & year >= 2005
  * ---- and/or: firm 1's outcome is missing from 2005 on (rows leave e(sample))
  if `missy' == 1 replace y = . if firm == 1 & year >= 2005
  sort firm year
end

prog def regimechk
  * faithful replication of ado:8988-8994 on the WHOLE data (no user if):
  *   smp = 1 & !missing(clusv);  nest = 1 iff min(clus)==max(clus) by i, if smp
  syntax varlist(min=2 max=2)
  tokenize `varlist'
  local iv `1'
  local cv `2'
  tempvar smp mn mx
  qui gen byte `smp' = 1
  qui replace `smp' = 0 if missing(`cv')
  qui gegen double `mn' = min(`cv') if `smp', by(`iv')
  qui gegen double `mx' = max(`cv') if `smp', by(`iv')
  cap assert `mn' == `mx' if `smp'
  c_local PREDNEST = (_rc == 0)
end

prog def getatti
  * pull att_i estimate + se out of e(b)/e(V)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  local j 0
  forvalues q = 1/`K' {
    local this : word `q' of `nm'
    if "`this'" == "att_i" local j `q'
  }
  c_local ATTI = `b'[1,`j']
  c_local ATTISE = sqrt(`V'[`j',`j'])
end

local REPS 150

foreach V in A B C C0 {
  if "`V'" == "A"  mkpanel
  if "`V'" == "B"  mkpanel , straddle(1)
  if "`V'" == "C"  mkpanel , straddle(1) missy(1)
  if "`V'" == "C0" mkpanel , missy(1)
  save "p01_`V'.dta", replace
  regimechk firm ind
  dis "VARIANT=`V' PREDNEST=`PREDNEST'"
  timer clear 1
  timer on 1
  cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(`REPS') rseed(4242) dots(0)
  dis "VARIANT=`V' RC=" _rc
  timer off 1
  timer list 1
  getatti
  dis "VARIANT=`V' ATTI=" %14.10f `ATTI' "  BOOTSE=" %14.10f `ATTISE' "  NCLUST=" e(N_clust) "  N=" e(N)
  * analytic reference
  cap noisily wooldid y firm year onset, att_i cluster(ind)
  getatti
  dis "VARIANT=`V' ANALYTIC ATTI=" %14.10f `ATTI' "  ANSE=" %14.10f `ATTISE'
}
dis "ALLDONE"
