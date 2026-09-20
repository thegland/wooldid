clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P06 -- the regime flip should move the CELL-WEIGHTED estimands only.
* Variant C0 (nesting regime) vs variant C (keep-i regime).  The two
* datasets have the SAME estimation sample and the SAME point estimates;
* the same rseed gives them the SAME cluster draws.  If the flip is what
* moves the number, att (observation-weighted) must be unchanged to the
* last digit while att_i / att_it / att_itime move.
*----------------------------------------------------------------------
prog def dumpse
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if inlist("`this'","main_att","main_att_i","main_att_it","main_att_itime") {
      dis "P06 `tag' `this' b=" %18.14f `b'[1,`q'] " se=" %18.14f sqrt(`V'[`q',`q'])
    }
  }
end

foreach V in C0 C {
  use "p01_`V'.dta", clear
  cap noisily wooldid y firm year onset, att att_i att_it att_itime cluster(ind) bootstrap(standard) reps(120) rseed(31337) dots(0)
  dis "P06 VARIANT=`V' RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
  dumpse , tag(`V')
}

*----------------------------------------------------------------------
* e(sample) is a SUBSET of the detection sample smp on variant C, so the
* detection can only err in one direction: report NOT-NESTING where the
* estimation sample nests.  Shown, not asserted in prose.
*----------------------------------------------------------------------
use "p01_C.dta", clear
qui wooldid y firm year onset, att_i cluster(ind)
gen byte es = e(sample)
gen byte smp = 1
replace smp = 0 if missing(ind)
count if es == 1 & smp == 0
dis "P06 ROWS IN e(sample) BUT NOT IN smp = " r(N) "   (must be 0)"
count if es == 0 & smp == 1
dis "P06 ROWS IN smp BUT NOT IN e(sample) = " r(N) "   (the phantom-straddle rows)"
list firm year ind y es smp if es == 0 & smp == 1, noobs sepby(firm)
* and those rows are exactly what makes the detection say NOT-NESTING
tempvar mn mx
qui gegen double `mn' = min(ind) if smp, by(firm)
qui gegen double `mx' = max(ind) if smp, by(firm)
cap assert `mn' == `mx' if smp
dis "P06 NEST ON smp        = " (_rc == 0)
cap assert `mn' == `mx' if es
drop `mn' `mx'
qui gegen double `mn' = min(ind) if es, by(firm)
qui gegen double `mx' = max(ind) if es, by(firm)
cap assert `mn' == `mx' if es
dis "P06 NEST ON e(sample)  = " (_rc == 0)

dis "ALLDONE"
