* p06: (a) EXTERNAL ORACLE vs Stata 19's official rwgen/bayesboot
*      (b) the bootalpha(size) arm of the GOF at R=1000
*      (c) massive-R draws from my verified reimplementation, for a high-power GOF
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
dis "STATA_VERSION=" c(stata_version)

mata:
void myb_dw(real colvector n, string scalar amode, real scalar R, string scalar outname)
{
  real matrix E, P
  real scalar C, r, k, a
  C = rows(n)
  E = J(C, R, .)
  for (r = 1; r <= R; r++) {
    for (k = 1; k <= C; k++) {
      a = (amode == "size" ? n[k] : 1)
      E[k, r] = rgamma(1, 1, a, 1)
    }
  }
  P = E :/ colsum(E)
  if (amode == "size") P = P :/ n
  st_matrix(outname, P)
}
end

*==== (b) size arm, R=1000, real engine ====
use p05_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(1000) rseed(20260917) bootsave(p05_size) replace dots(0)
dis "RUN_size_RC=" _rc

*==== (a) rwgen oracle ====
local R 20
local SEED 31415
use p05_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(`R') rseed(`SEED') bootsave(p06_size) replace dots(0)
dis "WD_SIZE_RC=" _rc
use p05_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(equal) reps(`R') rseed(`SEED') bootsave(p06_equal) replace dots(0)
dis "WD_EQ_RC=" _rc

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
qui export delimited using p06_rwgen.csv, replace

*==== (c) massive-R reimplementation draws ====
use p05_panel, clear
qui bysort st (year): gen long nc = _N
qui bysort st (year): keep if _n==1
qui sort st
mkmat nc, matrix(NCv)
foreach am in equal size {
  set seed 90210
  mata: myb_dw(st_matrix("NCv"), "`am'", 50000, "BIG")
  clear
  qui svmat double BIG
  qui export delimited using p06_big_`am'.csv, replace novarnames
  dis "BIG_`am' rows=" _N " cols=" c(k)
}
dis "PROBE_DONE"
