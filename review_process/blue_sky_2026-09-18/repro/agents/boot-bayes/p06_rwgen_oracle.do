* p06: EXTERNAL ORACLE. wooldid.ado claims bootalpha(size) "matches bayesboot
* BITWISE on unweighted data" (ado:8899-8900) and the help says it reproduces
* "what bayesboot does when passed id(clustervar)" (sthlp:548). Stata 19's
* official rwgen/bayesboot is installed. Test it.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
dis "STATA_VERSION=" c(stata_version)

local R 20
local SEED 31415

*--- wooldid side: bootalpha(size), unweighted ------------------------------
use p05_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(`R') rseed(`SEED') bootsave(p06_size) replace dots(0)
dis "WD_SIZE_RC=" _rc
use p05_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(equal) reps(`R') rseed(`SEED') bootsave(p06_equal) replace dots(0)
dis "WD_EQ_RC=" _rc

*--- rwgen side -------------------------------------------------------------
use p05_panel, clear
qui bysort st (year): gen long nc = _N
qui gen double negnc = -nc
set seed `SEED'
cap noi rwgen bayes RS, reps(`R') id(st)
dis "RWGEN_SIZE_RC=" _rc
set seed `SEED'
cap noi rwgen bayes RE, reps(`R') id(st) priorpowers(negnc)
dis "RWGEN_EQ_RC=" _rc
qui bysort st (year): keep if _n == 1
qui sort st
qui keep st nc RS* RE*
tempfile rw
qui save `rw'

*--- compare ---------------------------------------------------------------
foreach am in size equal {
  use p06_`am', clear
  qui keep if _rep >= 1
  qui sort _rep
  qui ds w_*
  mkmat `r(varlist)', matrix(WD_`am')     // R x C
}
use `rw', clear
mkmat RS1-RS`R', matrix(RS)               // C x R
mkmat RE1-RE`R', matrix(RE)
mkmat nc, matrix(NCv)
mata {
  R = `R'
  WS = st_matrix("WD_size")   // R x C
  WE = st_matrix("WD_equal")
  S  = st_matrix("RS")'       // R x C
  E  = st_matrix("RE")'
  nn = st_matrix("NCv")
  printf("SIZE  : max|wooldid - rwgen(id)|            = %g   bitwise-equal=%f\n", max(abs(WS-S)), (WS==S))
  printf("SIZE  : max rel diff                        = %g\n", max(abs((WS-S):/S)))
  printf("EQUAL : max|wooldid - rwgen(pp=-n_c)*n_c|   = %g   bitwise-equal=%f\n", max(abs(WE - E:*nn')), (WE==(E:*nn')))
  printf("EQUAL : max rel diff                        = %g\n", max(abs((WE - E:*nn'):/(E:*nn'))))
  // distributional fallback: do the two engines agree in DISTRIBUTION?
  printf("SIZE  : colmean wooldid vs rwgen  maxdiff   = %g\n", max(abs(mean(WS)-mean(S))))
  printf("SIZE  : rowsum(n_c*w) wooldid=%g rwgen=%g\n", (WS[1,.]*nn), (S[1,.]*nn))
  printf("EQUAL : rowsum(w) wooldid=%g  rwgen*n=%g\n", sum(WE[1,.]), sum(E[1,.]:*nn'))
}
dis "PROBE_DONE"
