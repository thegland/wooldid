* p02: REIMPLEMENT the Dirichlet draw process and test it bitwise against the
*      multipliers wooldid actually applied (dumped via bootsave()).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

*---------------------------------------------------------------------------
* my OWN implementation of the documented process, written from the help text:
*   per-replication (outer) x per-cluster-ascending (inner) Gamma(alpha,1),
*   normalized within replication; under bootalpha(size) also divided by n_c.
*---------------------------------------------------------------------------
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

use p01_panel, clear
qui bysort st (year): gen long nrow = _N
qui bysort st (year): keep if _n == 1
qui sort st
mkmat nrow, matrix(NC)
local C = _N
dis "C=`C'"
mata: st_matrix("NCv", st_matrix("NC"))

foreach am in equal size {
  use p01_panel, clear
  local rp 60
  cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(`am') reps(`rp') rseed(4242) bootsave(p02_`am') replace dots(0)
  dis "RUN_`am'_RC=" _rc
  dis "Nrows_est=" e(N)

  * my draws, same seed, same order
  set seed 4242
  mata: myb_dw(st_matrix("NCv"), "`am'", `rp', "MYP")

  use p02_`am', clear
  qui ds w_*
  local wv `r(varlist)'
  dis "NWCOLS_`am'=" `: word count `wv''
  * build R x C matrix from file rows _rep>=1
  qui keep if _rep >= 1
  qui sort _rep
  mkmat `wv', matrix(WF)
  mata {
    WF = st_matrix("WF")          // R x C
    MP = st_matrix("MYP")'        // R x C
    d  = WF - MP
    printf("`am': maxabsdiff = %21x\n", max(abs(d)))
    printf("`am': exact bitwise equal = %f\n", (WF == MP))
    printf("`am': rowsum(WF) first 3: %g %g %g\n", sum(WF[1,.]), sum(WF[2,.]), sum(WF[3,.]))
    nn = st_matrix("NCv")
    printf("`am': sum_c n_c*w_c rep1 = %21x\n", WF[1,.]*nn)
  }
  * observed row
  use p02_`am', clear
  qui keep if _rep == 0
  foreach v of local wv {
     local wobs0 = `v'[1]
     continue, break
  }
  dis "`am': OBSROW_w = " %21x `wobs0'
}
dis "DONE"
