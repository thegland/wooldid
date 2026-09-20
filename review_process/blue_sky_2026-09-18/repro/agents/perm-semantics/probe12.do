* probe12 -- the two surfaces the help says carry NO permutation inference:
*  (a) e(gteffects): "keyed by cohort onset - exactly what a replication reassigns - so
*      they have no permutation distribution and their inference columns are blank";
*  (b) the continuous-treatment plot cells: permuted for assignment-invariant lattice
*      types, blanked (cols 4-6 missing always, col 3 missing unless permtstat) for the
*      baseTreated*/byTreated*/byTtre* types, col 2 always the analytic se.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
local dN 24
local dT 8
local dons 2004 2006
local dnever 8
local dtau 1.2
local dsd 0.7
local dseed 606
include "dgp.do"
* a unit-constant dose for the continuous-treatment arm
bysort unit: gen double dose0 = .
by unit: replace dose0 = runiform()*3 if _n == 1
by unit: replace dose0 = dose0[1]
gen double dose = cond(post, dose0, 0)
tempfile base
save `base'

* ---------- (a) gteffects under permute ----------
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) gteffects permute(shuffle) permreps(49) rseed(37)
dis "RC_GT=" _rc
dis "perm_gt_blank = " e(perm_gt_blank)
matrix GT = e(gteffects)
dis "gteffects dims: " rowsof(GT) " x " colsof(GT)
mata:
  G = st_matrix("GT")
  cn = st_matrixcolstripe("GT")
  for (j=1;j<=cols(G);j++) printf("gtcol %g = %s  nonmissing=%g\n", j, cn[j,2], sum(G[.,j] :< .))
end
matrix list GT, format(%10.6f)

* ---------- (b) continuous treatment cells ----------
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) contreat(dose) contreatcontrols(none) lattice(all 3) makecfxplots cfxplottypes(all baseTreated) permute(shuffle) permreps(25) rseed(43)
dis "RC_CT=" _rc
dis "perm_blank = " e(perm_blank) "  perm_cfx_blank = " e(perm_cfx_blank)
ereturn list
dis "PROBE12_DONE"
