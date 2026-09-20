clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P10 -- does the REPORTED att standard error depend on which OTHER
*        estimands were requested?
*
* `att' alone with cluster()      -> the frequency-reweighting COUNT leg
* `att att_i' with cluster()      -> the physical RE-ID leg (ado:378)
*
* The help says the two are "the same statistic"; with the same rseed the
* draws are the same (ado:9295: "map-frame draws == full-data draws under
* the cluster()+idcluster() form, measured").  So att's bootstrap SE
* should not move.  Run it under BOTH regimes.
*   p01_A.dta  = clean nesting  -> nest regime on the reid leg
*   p01_B.dta  = one straddler  -> keep-i regime on the reid leg
*----------------------------------------------------------------------
prog def getb
  syntax , nm(string) tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  local j 0
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if "`this'" == "`nm'" local j `q'
  }
  if `j' > 0 dis "P10 `tag' `nm' b=" %18.14f `b'[1,`j'] " se=" %18.14f sqrt(`V'[`j',`j'])
  else dis "P10 `tag' `nm' ABSENT"
end

foreach D in A B {
  use "p01_`D'.dta", clear
  cap noisily wooldid y firm year onset, att cluster(ind) bootstrap(standard) reps(200) rseed(6161) dots(0)
  dis "P10 `D' COUNTLEG RC=" _rc
  getb , nm(main_att) tag(`D'_COUNT)

  use "p01_`D'.dta", clear
  cap noisily wooldid y firm year onset, att att_i cluster(ind) bootstrap(standard) reps(200) rseed(6161) dots(0)
  dis "P10 `D' REIDLEG RC=" _rc
  getb , nm(main_att)   tag(`D'_REID)
  getb , nm(main_att_i) tag(`D'_REID)
}
dis "ALLDONE"
