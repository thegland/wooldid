* probe06 -- ARE THE DRAWS WHAT THEY CLAIM TO BE?
* The help says each draw "is the estimate wooldid would compute under that
* counterfactual assignment".  Take replications out of the permsave dump, rebuild
* ttre from its onset_* columns, refit wooldid with NO permute(), and compare.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
mata:
void cmprep(real scalar rr, string scalar rbname)
{
  real rowvector rb, dr
  real scalar j, d, worst
  string matrix nm
  rb = st_matrix(rbname)
  dr = st_matrix("DUMPB")[rr+1, .]
  nm = st_matrixcolstripe(rbname)
  worst = 0
  for (j=1;j<=cols(rb);j++) {
    d = reldif(rb[1,j], dr[1,j])
    if (rb[1,j]>=. & dr[1,j]>=.) d = 0
    if (d > worst) worst = d
    printf("  rep %g col %g %-20s refit=%21x dump=%21x reldif=%g\n", rr, j, nm[j,2], rb[1,j], dr[1,j], d)
  }
  printf("REFIT rep=%g WORST reldif = %g\n", rr, worst)
}
end
include "dgp.do"
tempfile base
save `base'
use d_def.dta, clear
qui ds onset_*
local ons `r(varlist)'
dis "dump: " _N " rows"
qui ds onset_* _rep, not
local bc `r(varlist)'
mata: st_matrix("DUMPB", st_data(., "`bc'"))
mata: st_matrix("DUMPO", st_data(., "`ons'"))
local clusnames
foreach v of local ons {
  local clusnames `clusnames' `=subinstr("`v'", "onset_", "", 1)'
}
dis "cluster values: `clusnames'"
foreach rr in 0 1 2 7 23 {
  use `base', clear
  local k 0
  foreach cv of local clusnames {
    local ++k
    local newons = DUMPO[`=`rr'+1', `k']
    qui replace ttre = `newons' if unit == `cv'
  }
  qui count if !missing(ttre)
  local ntr = r(N)
  cap noi wooldid y unit t ttre, cluster(unit) att esprelength(3) espostlength(3) fastnose
  dis "REFIT rep=`rr' rc=" _rc " treatedrows=`ntr'"
  if _rc == 0 {
    matrix RB = e(b)
    mata: cmprep(`rr', "RB")
  }
}
* also: a full (non-fast) refit of the observed assignment, to see whether fast mode
* moves any point estimate
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(3) espostlength(3)
dis "FULLREFIT rc=" _rc
matrix RBF = e(b)
mata: cmprep(0, "RBF")
dis "PROBE06_DONE"
