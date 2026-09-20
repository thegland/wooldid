clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P04 -- MEASURE the SE consequence of the two-regime flip.
* Variants A/B/C/C0 built by p01_regime.do (saved as p01_*.dta):
*   A  : clean nesting (firm inside industry)            -> nest regime
*   B  : firm 1 really moves industry in 2005            -> keep-i regime
*   C  : firm 1 moves industry in 2005 AND y is missing
*        for firm 1 from 2005 on (those rows are NOT in
*        e(sample))                                      -> keep-i regime
*   C0 : y missing for firm 1 from 2005 on, no move      -> nest regime
* C and C0 have the SAME estimation sample and the SAME point estimate.
*----------------------------------------------------------------------
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

local REPS 300

foreach V in A B C C0 {
  use "p01_`V'.dta", clear
  cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(`REPS') rseed(4242) dots(0)
  dis "P04 VARIANT=`V' RC=" _rc
  getatti
  dis "P04 VARIANT=`V' ATTI=" %18.14f `ATTI' " BOOTSE=" %18.14f `ATTISE' " N=" e(N) " NCLUST=" e(N_clust)
  cap noisily wooldid y firm year onset, att_i cluster(ind)
  getatti
  dis "P04 VARIANT=`V' ANALYTIC ATTI=" %18.14f `ATTI' " ANSE=" %18.14f `ATTISE'
}

*----------------------------------------------------------------------
* INDEPENDENT ORACLE.  Textbook pairs-cluster bootstrap on the C/C0 data
* (identical estimation samples), done two ways by hand:
*   ORA_NEST : draw industries with replacement, mint fresh firm ids
*              newi = group(drawid, firm)              [correct here]
*   ORA_KEEP : draw industries with replacement, KEEP firm ids
*                                                      [what wooldid does
*                                                       in variant C]
*----------------------------------------------------------------------
foreach RULE in NEST KEEP {
  use "p01_C0.dta", clear
  qui wooldid y firm year onset, att_i fastnose
  getatti
  local b0 = `ATTI'
  tempname acc
  mat `acc' = J(`REPS', 1, .)
  set seed 90210
  set sortseed 90210
  forvalues r = 1/`REPS' {
    preserve
    qui bsample, cluster(ind) idcluster(__ns)
    if "`RULE'" == "NEST" qui gegen double __newi = group(__ns firm)
    else qui gen double __newi = firm
    cap qui wooldid y __newi year onset, att_i fastnose
    if _rc == 0 {
      getatti
      mat `acc'[`r',1] = `ATTI'
    }
    restore
  }
  preserve
  clear
  qui svmat double `acc'
  qui sum `acc'1
  dis "P04 ORACLE_`RULE' b0=" %18.14f `b0' " n=" r(N) " mean=" %18.14f r(mean) " SD=" %18.14f r(sd)
  restore
}
dis "ALLDONE"
