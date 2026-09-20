* probe12b -- continuous treatment under permute(): the help makes very specific claims
* about the e(cfx*) plot-cell matrices.  Column 1 the observed estimate on every type;
* column 2 the ANALYTIC se under the plain name se; on permuted (assignment-invariant)
* types columns 4-6 the permutation p and shift interval and column 3 the t under
* permtstat; on assignment-DEPENDENT types (baseTreated*, byTreated*, byTtre*) columns
* 4-6 missing always and column 3 missing unless permtstat.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
mata:
void cfxscan(string scalar mn)
{
  real matrix M
  real scalar j
  string matrix cn
  M = st_matrix("e(" + mn + ")")
  if (rows(M) == 0) {
    printf("CFX %s : ABSENT\n", mn)
    return
  }
  cn = st_matrixcolstripe("e(" + mn + ")")
  printf("CFX %s : %g x %g\n", mn, rows(M), cols(M))
  for (j=1;j<=cols(M);j++) printf("   col %g %-20s nonmissing=%g / %g\n", j, cn[j,2], sum(M[.,j] :< .), rows(M))
}
end
* 24 units in 6 groups of 4; i = group (so a treated cell holds four units with
* different doses), cluster = group
clear
set obs 24
gen long unit = _n
gen long grp = ceil(unit/4)
gen double ttre = .
replace ttre = 2004 if inlist(grp,1,2)
replace ttre = 2006 if inlist(grp,3,4)
gen double ai = rnormal()
gen double dose0 = runiform()*3
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double dose = dose0
gen double y = ai + 0.2*(t-2000) + 0.8*dose0*post + 0.5*rnormal()
tempfile base
save `base'
use `base', clear
cap noi wooldid y grp t ttre, cluster(grp) att esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) lattice(all 3) makecfxplots cfxplottypes(post) nofastcollapse permute(shuffle) permreps(25) rseed(43)
dis "RC_CT=" _rc
dis "perm_blank = " e(perm_blank) "  perm_cfx_blank = " e(perm_cfx_blank)
local enames : e(macros)
ereturn list
local ms : e(matrices)
foreach mn of local ms {
  if strpos("`mn'", "cfx") | strpos("`mn'", "permdist_") | strpos("`mn'", "permcmp_") {
    mata: cfxscan("`mn'")
  }
}
local htag CONTREAT
local habsp 0
local hexact 0
include "harvest.do"
matrix list e(perm_n)

* ---- second run: an ASSIGNMENT-DEPENDENT lattice type, which must be blanked ----
use `base', clear
cap noi wooldid y grp t ttre, cluster(grp) att esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) lattice(baseTreated 3) makecfxplots cfxplottypes(post) nofastcollapse permute(shuffle) permreps(25) rseed(43)
dis "RC_CT2=" _rc
dis "perm_blank = " e(perm_blank) "  perm_cfx_blank = " e(perm_cfx_blank)
local ms2 : e(matrices)
foreach mn of local ms2 {
  if strpos("`mn'", "cfx") | strpos("`mn'", "permdist_") | strpos("`mn'", "permcmp_") {
    mata: cfxscan("`mn'")
  }
}
local htag CONTREAT2
local habsp 0
local hexact 0
include "harvest.do"
dis "PROBE12B_DONE"
