clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* ---- tiny unbalanced staggered panel -------------------------------------
clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 8
replace ttre = 6 if i > 8 & i <= 16
expand 8
bys i: gen t = _n
gen ui = rnormal()
gen ut = rnormal()/3
bys i (t): replace ui = ui[1]
bys t (i): replace ut = ut[1]
gen te = 0
replace te = 0.5 + 0.2*(t - ttre) if t >= ttre & ttre < .
gen y = ui + ut + te + rnormal()
gen w = 1 + mod(i,3)
gen cwm = 1 + mod(i+t,4)
* make it UNBALANCED so the four estimands genuinely differ
drop if i == 3 & t == 7
drop if i == 11 & t == 8
drop if i == 12 & t >= 7
drop if i == 20 & t == 2

* ---- helper: aggregation check -------------------------------------------
cap prog drop aggchk
prog def aggchk
  syntax , tag(string)
  tempname G B
  matrix `G' = e(gteffects)
  matrix `B' = e(b)
  local cn : colnames `B'
  mata: wdaggchk("`G'", "`B'", "`cn'", "`tag'")
end

mata:
void wdaggchk(string scalar gn, string scalar bn, string scalar cn, string scalar tag)
{
  real matrix G, B
  string rowvector nm, wnm
  real scalar j, k, c, num, den, got, want, dd
  G = st_matrix(gn)
  B = st_matrix(bn)
  nm = tokens(cn)
  wnm = ("att","att_it","att_i","att_itime")
  for (k=1;k<=4;k++) {
    c = 11 + k
    num = 0
    den = 0
    for (j=1;j<=rows(G);j++) {
      if (G[j,4] >= 0 & G[j,c] < .) {
        num = num + G[j,c]*G[j,5]
        den = den + G[j,c]
      }
    }
    got = .
    if (den != 0) got = num/den
    want = .
    for (j=1;j<=cols(B);j++) {
      if (nm[j] == "main_" + wnm[k]) want = B[1,j]
    }
    dd = abs(got-want)/max((abs(want),1e-30))
    printf("%s POST %s agg=%21x eb=%21x reldiff=%g\n", tag, wnm[k], got, want, dd)
    num = 0
    den = 0
    for (j=1;j<=rows(G);j++) {
      if (G[j,4] < 0 & G[j,c] < .) {
        num = num + G[j,c]*G[j,5]
        den = den + G[j,c]
      }
    }
    got = .
    if (den != 0) got = num/den
    want = .
    for (j=1;j<=cols(B);j++) {
      if (nm[j] == "pretreatment_" + wnm[k]) want = B[1,j]
    }
    dd = abs(got-want)/max((abs(want),1e-30))
    printf("%s PRE_ %s agg=%21x eb=%21x reldiff=%g\n", tag, wnm[k], got, want, dd)
  }
}
end

dis "SENTINEL_A_START"
cap noi wooldid y i t ttre , att att_it att_i att_itime gteffects esprelength(2) espostlength(6) cluster(i)
dis "rc=" _rc
aggchk , tag(A_unw)
dis "SENTINEL_A_DONE"

dis "SENTINEL_B_START"
cap noi wooldid y i t ttre [aw=w] , att att_it att_i att_itime gteffects esprelength(2) espostlength(6) cluster(i)
dis "rc=" _rc
aggchk , tag(B_aw)
dis "SENTINEL_B_DONE"

dis "SENTINEL_C_START"
cap noi wooldid y i t ttre [aw=w] , att att_it att_i att_itime gteffects esprelength(2) espostlength(6) cluster(i) customweightmultiplier(cwm)
dis "rc=" _rc
aggchk , tag(C_cwm)
dis "SENTINEL_C_DONE"

* ---- ES vs overall, and e(b) vs e(b_es_*) --------------------------------
dis "SENTINEL_D_START"
cap noi wooldid y i t ttre , att att_it gteffects esprelength(2) espostlength(6) cluster(i)
mat list e(b), format(%21x)
cap mat list e(b_es_att), format(%21x)
cap mat list e(b_es_att_it), format(%21x)
tempname V VE
mat `V' = e(V)
cap mat `VE' = e(V_es_att)
cap mat list `VE', format(%21x)
dis "SENTINEL_D_DONE"
