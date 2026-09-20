clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* ===========================================================================
* MC coverage of the `impvalues' (Ibragimov-Muller) confidence interval when
* the cluster CUTS ACROSS i and the clusters have near-identical composition.
* Reference arm: the ordinary analytic clustered interval on the SAME draws.
* ===========================================================================
local R 250
local ISEED 1
include "$WD/mkimp.doi"
replace y = base
qui wooldid y cg t ttre, att cluster(reg) impvalues level(95)
local TRUTH_IM = e(b)[1,1]
qui wooldid y cg t ttre, att cluster(reg) level(95)
local TRUTH_AN = e(b)[1,1]
dis "TRUTH|im=" %14.9f `TRUTH_IM' "|an=" %14.9f `TRUTH_AN'

tempname fh
file open `fh' using "$WD/mc_imp.txt", write replace text
file write `fh' "draw|arm|b|se|lb|ub|truth" _n
forvalues r = 1/`R' {
  local ISEED = 400000 + `r'
  include "$WD/mkimp.doi"
  cap qui wooldid y cg t ttre, att cluster(reg) impvalues level(95)
  if _rc==0 {
    mat Q = e(wdidmainresults)
    file write `fh' "`r'|im|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TRUTH_IM') _n
  }
  else {
    file write `fh' "`r'|im|.|.|.|.|" (`TRUTH_IM') _n
  }
  cap qui wooldid y cg t ttre, att cluster(reg) level(95)
  if _rc==0 {
    mat Q = e(wdidmainresults)
    file write `fh' "`r'|an|" (Q[1,1]) "|" (Q[1,2]) "|" (Q[1,5]) "|" (Q[1,6]) "|" (`TRUTH_AN') _n
  }
  if mod(`r',50)==0 dis "IMPROG|`r'"
}
file close `fh'
dis "MC_IMP_DONE"
