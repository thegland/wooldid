*! SELF-CONTAINED REPRO -- within ONE permexact run, e(perm_sd) silently switches
*! denominator: 1/n (population) for elements estimable in all M assignments, 1/(n-1)
*! (sample) for elements that lost assignments.  The help documents no such switch,
*! and the p-value convention next to it was deliberately moved OFF completeness (BS-23).
*! Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_permsd_denominator.do
* see which one wooldid posted for a complete column and for an incomplete one.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
clear
set obs 8
gen long unit = _n
gen double ttre = cond(inlist(unit,1,2), 2004, .)
gen byte sg = cond(inlist(unit,1,5), 1, 0)
gen double ai = rnormal()
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ai + 0.25*(t-2000) + 2.0*post + 0.5*rnormal()
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(1) espostlength(1) subgroups(sg) permute(shuffle) permexact
dis "RC=" _rc
matrix _PD = e(permdist_b)
matrix _SD = e(perm_sd)
matrix _N2 = e(perm_n)
mata:
void sdwhich(real scalar M)
{
  real matrix D
  real rowvector S
  real colvector v
  real scalar j, n, sp, ss
  string scalar w
  string matrix nm
  D = st_matrix("_PD")
  S = st_matrix("_SD")
  nm = st_matrixcolstripe("_PD")
  for (j=1;j<=cols(D);j++) {
    v = select(D[.,j], D[.,j] :< .)
    n = rows(v)
    sp = sqrt(variance(v)*(n-1)/n)
    ss = sqrt(variance(v))
    w = "NEITHER"
    if (S[1,j] == sp) w = "POPULATION(1/n)"
    if (S[1,j] == ss) w = "SAMPLE(1/(n-1))"
    printf("SDROW %-24s n=%g M=%g posted=%21x pop=%21x samp=%21x %s\n", nm[j,2], n, M, S[1,j], sp, ss, w)
  }
}
end
mata: sdwhich(28)
dis "PROBE15B_DONE"
