clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P08 -- DOSE RESPONSE of the regime rule.
* The regime is a SINGLE GLOBAL FLAG, but "does this unit nest inside one
* cluster?" is a PER-UNIT property.  Reclassify k of the 96 firms into a
* neighbouring industry halfway through the panel and watch what the
* bootstrap standard error of att_i does as k goes 0 -> 1 -> 5 -> 20.
* If SE(k=1) == SE(k=5) == SE(k=20) << SE(k=0), the rule is a cliff: one
* mis-nesting firm out of 96 rewrites the inference for all 96.
*----------------------------------------------------------------------
prog def mkpanel
  syntax , [k(int 0)]
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
  * the FIRST firm of each of the first k industries moves to the NEXT industry
  * in 2005 (a routine SIC / district reclassification)
  if `k' > 0 {
    forvalues q = 1/`k' {
      qui replace ind = cond(`q' == 24, 1, `q' + 1) if firm == `=(`q'-1)*4 + 1' & year >= 2005
    }
  }
  sort firm year
end

prog def regimechk
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
  * how many units actually straddle
  tempvar bad
  qui gen byte `bad' = (`mn' != `mx') if `smp'
  qui gegen double `mn'2 = max(`bad'), by(`iv')
  qui egen byte tagx = tag(`iv')
  qui count if tagx == 1 & `mn'2 == 1
  c_local NSTRAD = r(N)
  drop tagx `mn'2
end

prog def getatti
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  local j 0
  local ja 0
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if "`this'" == "main_att_i" local j `q'
    if "`this'" == "main_att"   local ja `q'
  }
  c_local ATTI   = cond(`j'  > 0, `b'[1,max(`j',1)],  .)
  c_local ATTISE = cond(`j'  > 0, sqrt(`V'[max(`j',1),max(`j',1)]), .)
  c_local ATT    = cond(`ja' > 0, `b'[1,max(`ja',1)], .)
  c_local ATTSE  = cond(`ja' > 0, sqrt(`V'[max(`ja',1),max(`ja',1)]), .)
end

foreach k in 0 1 5 20 {
  mkpanel , k(`k')
  save "p08_k`k'.dta", replace
  regimechk firm ind
  dis "P08 k=`k' PREDNEST=`PREDNEST' NSTRADDLING_UNITS=`NSTRAD' of 96"
  cap noisily wooldid y firm year onset, att att_i cluster(ind) bootstrap(standard) reps(120) rseed(8081) dots(0)
  dis "P08 k=`k' RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
  getatti
  dis "P08 k=`k' ATT  b=" %18.14f `ATT'  " se=" %18.14f `ATTSE'
  dis "P08 k=`k' ATTI b=" %18.14f `ATTI' " se=" %18.14f `ATTISE'
}
dis "ALLDONE"
