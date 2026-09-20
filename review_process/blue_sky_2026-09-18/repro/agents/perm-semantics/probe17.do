* probe17 -- "Every reported permutation quantity is otherwise reconstructable from the
* [permsave] file."  Under permtstat the p ranks b_cmp/se_cmp, where se_cmp is the
* COMPARATOR run's standard error; the file's _rep=0 se_ row is e(perm_ase_obs), the
* observed analytic SE the TABLE shows.  If those two differ, the file cannot reproduce
* the reported p.  Test both, and reconstruct e(perm_p) from the FILE ALONE.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "dgp.do"
tempfile base
save `base'
use `base', clear
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(3) espostlength(3) permute(shuffle) permtstat permreps(99) rseed(7) permsave(r_ts) replace
dis "RC=" _rc
matrix ASE = e(perm_ase_obs)
matrix CSE = e(perm_comparator_se)
matrix CMB = e(perm_comparator_b)
matrix PP  = e(perm_p)
matrix OB  = e(b)
mata:
  a = st_matrix("ASE"); c = st_matrix("CSE")
  printf("e(perm_ase_obs) == e(perm_comparator_se) : %g   maxreldif=%g\n", (a==c), max(reldif(a,c)))
  for (j=1;j<=cols(a);j++) printf("  col %g ase=%21x cmpse=%21x\n", j, a[1,j], c[1,j])
end
* ---- now reconstruct the reported p from the FILE ALONE ----
preserve
use r_ts.dta, clear
qui ds onset_* _rep se_*, not
local bc `r(varlist)'
dis "coef columns from file: `bc'"
local sc
foreach v of local bc {
  local sc `sc' se_`v'
}
mata:
  st_view(Bv=., ., "`bc'")
  st_view(Sv=., ., "`sc'")
  st_view(Rv=., ., "_rep")
  B = Bv[.,.]; S = Sv[.,.]; R = Rv[.,.]
  b0 = select(B, R :== 0); s0 = select(S, R :== 0)
  Br = select(B, R :> 0);  Sr = select(S, R :> 0)
  P = st_matrix("PP")
  nb = 0
  for (k=1;k<=cols(B);k++) {
    t = J(rows(Br),1,.)
    for (r=1;r<=rows(Br);r++) if (Br[r,k]<. & Sr[r,k]<. & Sr[r,k]>0) t[r] = Br[r,k]/Sr[r,k]
    v = select(t, t :< .)
    n = rows(v)
    tobs = b0[1,k]/s0[1,k]
    m0 = sum(v :<= tobs); m1 = sum(v :>= tobs)
    p  = min((1, 2*min(((1+m0)/(n+1), (1+m1)/(n+1)))))
    printf("FILE-RECON col %g n=%g  file p=%21x  reported p=%21x  equal=%g\n", k, n, p, P[1,k], (p==P[1,k]))
    if (p != P[1,k]) nb = nb + 1
  }
  printf("FILE-RECON mismatches: %g\n", nb)
end
restore
dis "PROBE17_DONE"
