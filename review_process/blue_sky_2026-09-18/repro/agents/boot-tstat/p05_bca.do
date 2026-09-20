clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid

local gG 20
local gUnits 3
include "`WD'/gen_small.do"
save bca20, replace

*==================================================================
* Q1: bootdisplay(bca) -- recompute z0, a and the endpoints from e()
*     plus MY OWN delete-one-cluster jackknife.
*==================================================================
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) rseed(909) bootdisplay(bca) att esprelength(2) espostlength(2) dots(0)
dis "Q1_RC=" _rc
dis "Q1_JK=" e(boot_jk_completed)
mat ACC  = e(boot_accel)
mat CIBCA= e(boot_ci_bca)
mat PBCA = e(boot_p_bca)
mat Z0   = e(boot_z0)
mat CMP  = e(boot_comparator)
mat BD   = e(bootdist)
mat BOB  = e(b)
mat list ACC, format(%21x)
mat list Z0, format(%21x)
mat list CIBCA, format(%21x)
mat list PBCA, format(%21x)
mat MRQ1 = e(wdidmainresults)
mat list MRQ1, format(%21x)
local LVL = e(boot_level)

* --- my own delete-one-cluster jackknife, by aweight 0/1 (the engine's rule) ---
gen double jkw = 1
mat JK = J(20, colsof(BOB), .)
matrix colnames JK = `: colfullnames BOB'
forval g = 1/20 {
  qui replace jkw = (state != `g')
  cap qui wooldid y state t ttre [aweight=jkw], cluster(state) att esprelength(2) espostlength(2) fastnose
  if _rc == 0 {
    mat bb = e(b)
    local cn : colfullnames bb
    forval k = 1/`=colsof(BOB)' {
      local nk : word `k' of `: colfullnames BOB'
      local pos : list posof "`nk'" in cn
      if `pos' > 0 mat JK[`g', `k'] = bb[1, `pos']
    }
  }
  else dis "JKFAIL g=`g' rc=" _rc
}
mata:
  JK = st_matrix("JK"); ACC = st_matrix("ACC")
  for (k=1; k<=cols(JK); k++) {
    th = select(JK[.,k], JK[.,k] :< .)
    if (rows(th) < 3) {
      printf("Q1 col %f: <3 jk fits\n", k)
      continue
    }
    d2 = sum((mean(th) :- th):^2)
    d3 = sum((mean(th) :- th):^3)
    a  = d3/(6*d2^1.5)
    printf("Q1 col %f  a_mine=%18.14f  a_wooldid=%18.14f  reldiff=%g\n", k, a, ACC[1,k], abs(a-ACC[1,k])/max((abs(ACC[1,k]),1e-14)))
  }
end

* --- recompute z0, BCa endpoints and p from e(bootdist) ---
mata:
  A = st_matrix("BD"); ob = st_matrix("BOB"); cmp = st_matrix("CMP")
  ACC = st_matrix("ACC"); Z0 = st_matrix("Z0"); CI = st_matrix("CIBCA"); P = st_matrix("PBCA")
  lvl = strtoreal(st_local("LVL")); al = (100-lvl)/200
  maxd = 0
  for (k=1; k<=cols(A); k++) {
    use = selectindex(A[.,k] :< .)
    nn = length(use)
    v = sort(A[use,k],1)
    z0 = invnormal(sum(v :<= cmp[1,k])/nn)
    a  = ACC[1,k]
    if (a >= . | z0 >= .) {
      printf("Q1 col %f: family absent (a=%g z0=%g) CI=[%g,%g] p=%g\n", k, a, z0, CI[1,k], CI[2,k], P[1,k])
      continue
    }
    w = z0 + invnormal(al)
    plo = normal(z0 + w/(1-a*w))
    w = z0 + invnormal(1-al)
    phi = normal(z0 + w/(1-a*w))
    x = nn*plo
    if (x==floor(x) & x>=1) lo = (v[x]+v[min((x+1,nn))])/2
    else lo = v[max((min((ceil(x),nn)),1))]
    x = nn*phi
    if (x==floor(x) & x>=1) hi = (v[x]+v[min((x+1,nn))])/2
    else hi = v[max((min((ceil(x),nn)),1))]
    r0 = sum(v :<= 0)/nn
    zr = invnormal(r0) - z0
    q  = normal(zr/(1+a*zr) - z0)
    pp = min((1, 2*min((q,1-q))))
    d1 = abs(z0-Z0[1,k]); d2x = abs(lo-CI[1,k]); d3x = abs(hi-CI[2,k]); d4 = abs(pp-P[1,k])
    maxd = max((maxd, d1, d2x, d3x, d4))
    printf("Q1 col %f z0=%12.8f/%12.8f lo=%12.8f/%12.8f hi=%12.8f/%12.8f p=%10.6f/%10.6f  b=%12.8f  contains_b=%f\n",
       k, z0, Z0[1,k], lo, CI[1,k], hi, CI[2,k], pp, P[1,k], ob[1,k], (CI[1,k] <= ob[1,k] & ob[1,k] <= CI[2,k]))
  }
  printf("Q1_BCA_MAXABSDIFF=%g\n", maxd)
end

*==================================================================
* Q2: small G (8 clusters) -- provoke a degenerate BCa family
*==================================================================
clear
local gG 8
local gUnits 4
local gnever 2
local gonsets "2002 2003"
include "`WD'/gen_small.do"
dis "Q2_NOBS=" _N
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(299) rseed(31) bootdisplay(bca) att esprelength(2) espostlength(2) dots(0)
dis "Q2_RC=" _rc
dis "Q2_JK=" e(boot_jk_completed)
mat list e(boot_accel), format(%21x)
mat list e(boot_ci_bca), format(%21x)
mat list e(boot_p_bca), format(%21x)
mat Q2M = e(wdidmainresults)
mat Q2B = e(b)
mata:
  M = st_matrix("Q2M"); B = st_matrix("Q2B")
  nout = 0
  for (r=1; r<=rows(M); r++) {
    if (M[r,5] < . & M[r,6] < . & M[r,1] < .) {
      if (M[r,1] < M[r,5] | M[r,1] > M[r,6]) {
        nout++
        printf("Q2 ROW %f estimate %g OUTSIDE bca CI [%g,%g]\n", r, M[r,1], M[r,5], M[r,6])
      }
    }
  }
  printf("Q2_ESTOUTSIDE=%f\n", nout)
end
dis "P05_DONE"
