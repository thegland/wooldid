*! SKEPTIC repro A -- INDEPENDENT test of PERMSEM-02.
*! Claim under test: under permexact, e(perm_sd) uses the POPULATION SD (divisor n) for
*! elements estimable in all M enumerated assignments and the SAMPLE SD (divisor n-1)
*! for elements that lost assignments, inside one table, with no marker.
*! Own data, own geometry: 9 clusters, 3 treated at one onset -> M = C(9,3) = 84.
*! A subgroup on 3 of the 9 clusters -> SG1 estimable in 84-C(6,3)=64 assignments,
*! SG0 estimable in 84-C(3,3)=83 assignments, main in all 84.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pA_permsd.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4711
set sortseed 4711
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

* ---------------- own DGP: 9 firms, 10 years, 3 treated in 2015 ----------------
clear
set obs 9
gen long firm = _n
gen double onset = cond(inlist(firm,3,6,8), 2015, .)
* subgroup on three firms -- e.g. "large" firms
gen byte big = inlist(firm,2,4,9)
gen double fe_i = rnormal()*2
expand 10
bysort firm: gen int year = 2010 + _n - 1
gen byte post = (year >= onset) & !missing(onset)
gen double y = fe_i + 0.3*(year-2010) + 1.2*post + rnormal()

cap noi wooldid y firm year onset, cluster(firm) att esprelength(2) espostlength(2) ///
    subgroups(big) permute(shuffle) permexact
dis "PA_RC=" _rc
dis "PA_M=" e(perm_npossible) "  completed=" e(perm_completed) "  failed=" e(perm_failed)
dis "PA_incomplete = |" "`e(perm_incomplete)'" "|"
dis "PA_addone     = |" "`e(perm_addone)'" "|"
dis "PA_pconv      = " "`e(perm_pconv)'"

matrix _PD = e(permdist_b)
matrix _SD = e(perm_sd)
matrix _NN = e(perm_n)
matrix _PP = e(perm_p)
matrix _CB = e(perm_comparator_b)
scalar _M = e(perm_npossible)

mata:
void sdcheck(real scalar M)
{
  real matrix D
  real rowvector S, NN, PP, CB
  real colvector v
  real scalar j, n, sp, ss, m0, m1, pex, padd
  string scalar w, tag
  string matrix nm
  D  = st_matrix("_PD")
  S  = st_matrix("_SD")
  NN = st_matrix("_NN")
  PP = st_matrix("_PP")
  CB = st_matrix("_CB")
  nm = st_matrixcolstripe("_SD")
  printf("PA_ROWS_PERMDIST_B = %g\n", rows(D))
  for (j=1;j<=cols(D);j++) {
    v = select(D[.,j], D[.,j] :< .)
    n = rows(v)
    if (n < 2) continue
    sp = sqrt(variance(v)*(n-1)/n)
    ss = sqrt(variance(v))
    w = "NEITHER"
    if (S[1,j] == sp) w = "POP_1_over_n"
    if (S[1,j] == ss) w = "SAMPLE_1_over_nm1"
    if (S[1,j] == sp & S[1,j] == ss) w = "DEGENERATE_EQUAL"
    tag = (n == M ? "COMPLETE" : "INCOMPLETE")
    // p conventions on the same column
    m0 = sum(v :<= CB[1,j]); m1 = sum(v :>= CB[1,j])
    pex  = min((1, 2*min((m0/n, m1/n))))
    padd = min((1, 2*min(((1+m0)/(n+1), (1+m1)/(n+1)))))
    printf("PA_SD %-26s n=%g eN=%g M=%g %-10s posted=%s pop=%s samp=%s %-18s ratio=%12.9f\n",
           nm[j,2], n, NN[1,j], M, tag, strofreal(S[1,j],"%21x"), strofreal(sp,"%21x"), strofreal(ss,"%21x"), w, (ss/sp))
    printf("PA_P  %-26s p=%s exactshare=%s addone=%s %s\n",
           nm[j,2], strofreal(PP[1,j],"%21x"), strofreal(pex,"%21x"), strofreal(padd,"%21x"),
           (PP[1,j]==pex ? "P_IS_EXACT_SHARE" : (PP[1,j]==padd ? "P_IS_ADDONE" : "P_IS_NEITHER")))
  }
}
end
mata: sdcheck(st_numscalar("_M"))

* the printed footer twin for the pooled estimand
dis "PA_FOOTER_SD_att = " %21x e(perm_sd_att)
dis "PA_MATRIX_SD_main_att = " %21x _SD[1, colnumb(_SD, "main_att")]

dis "PA_DONE"
