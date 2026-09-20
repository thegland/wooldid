* p14: the rwgen oracle again, but saved as .dta so the comparison is BITWISE
* (p06 compared through a CSV, which truncates).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
local R 20
local SEED 31415
use p05_panel, clear
qui bysort st (year): gen long nc = _N
qui gen double negnc = -nc
set seed `SEED'
qui rwgen bayes RS, reps(`R') id(st)
set seed `SEED'
qui rwgen bayes RE, reps(`R') id(st) priorpowers(negnc)
qui bysort st (year): keep if _n == 1
qui sort st
qui keep st nc RS* RE*
mkmat RS1-RS`R', matrix(RS)
mkmat RE1-RE`R', matrix(RE)
mkmat nc, matrix(NCv)
foreach am in size equal {
  use p06_`am', clear
  qui keep if _rep >= 1
  qui sort _rep
  qui ds w_*
  mkmat `r(varlist)', matrix(WD`am')
}
mata {
  WS = st_matrix("WDsize"); WE = st_matrix("WDequal")
  S = st_matrix("RS")'; E = st_matrix("RE")'; nn = st_matrix("NCv")
  printf("SIZE  bitwise equal to rwgen bayes id(): %f   max|diff| = %21x\n", (WS==S), max(abs(WS-S)))
  printf("EQUAL bitwise equal to rwgen*n_c      : %f   max|diff| = %21x\n", (WE==(E:*nn')), max(abs(WE-E:*nn')))
  printf("EQUAL max REL diff = %g\n", max(abs((WE-E:*nn'):/(E:*nn'))))
}
dis "PROBE_DONE"
