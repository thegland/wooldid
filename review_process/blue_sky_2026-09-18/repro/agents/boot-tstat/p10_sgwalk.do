clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid
local gG 24
local gUnits 3
local gT 7
local gnever 6
include "`WD'/gen_small.do"
gen byte reg2 = mod(state,2)
dis "NOBS=" _N

cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(199) boottstat rseed(1717) att esprelength(2) espostlength(2) subgroups(reg2) bootsave(d10) replace dots(0)
dis "M_RC=" _rc
local SGL `e(wdid_sglist)'
dis "SGLIST=`SGL'"
mat MR  = e(wdidmainresults)
mat TS  = e(boot_tstud)
mat ASE = e(boot_ase_obs)
mat BB  = e(b)
local LVL = e(boot_level)
local NMAIN = rowsof(MR)
dis "NMAIN=`NMAIN' K=" colsof(BB)
* stack the subgroup matrices below the main one in sglist order
mat STACK = MR
foreach sg of local SGL {
  mat SG`sg' = e(wdidsgresults_`sg')
  dis "SG`sg' rows=" rowsof(SG`sg')
  mat STACK = STACK \ SG`sg'
}
dis "STACKROWS=" rowsof(STACK)
mat list STACK, format(%21x)

preserve
use d10, clear
mata:
  ST = st_matrix("STACK"); TS = st_matrix("TS"); ASE = st_matrix("ASE"); BB = st_matrix("BB")
  cs = st_matrixcolstripe("BB")
  rs = st_matrixrowstripe("STACK")
  lvl = strtoreal(st_local("LVL")); al = (100-lvl)/200
  K = cols(BB)
  if (rows(ST) != K) printf("ROWCOUNT_MISMATCH stack=%f K=%f\n", rows(ST), K)
  worst_p = 0; worst_lb = 0; worst_ub = 0; worst_t = 0; worst_se = 0
  nblank = 0
  for (k = 1; k <= K; k++) {
    vn = strtoname(cs[k,2])
    sn = substr("se_" + vn, 1, 32)
    if (_st_varindex(vn) >= . | _st_varindex(sn) >= .) {
      printf("MISSINGVAR k=%f %s\n", k, vn)
      continue
    }
    B = st_data(., (vn, sn, "_rep"))
    o = select(B, B[.,3] :== 0)
    bobs = o[1,1]; seobs = o[1,2]
    R = select(B, B[.,3] :> 0)
    ok = (R[.,1] :< .) :& (R[.,2] :< .) :& (R[.,2] :> 0)
    if (sum(ok) < 2) {
      nblank++
      continue
    }
    Rk = select(R, ok)
    tr = sort((Rk[.,1] :- bobs) :/ Rk[.,2], 1)
    n = rows(tr); tobs = bobs/seobs
    x = n*al
    if (x == floor(x) & x >= 1) qlo = (tr[x] + tr[min((x+1,n))])/2
    else qlo = tr[max((min((ceil(x),n)),1))]
    x = n*(1-al)
    if (x == floor(x) & x >= 1) qhi = (tr[x] + tr[min((x+1,n))])/2
    else qhi = tr[max((min((ceil(x),n)),1))]
    m0 = sum(tr :<= tobs); m1 = sum(tr :>= tobs)
    pp = min((1, 2*(min((m0,m1))+1)/(n+1)))
    lb = bobs - qhi*seobs
    ub = bobs - qlo*seobs
    // col 2 of the table under boottstat is the OBSERVED ANALYTIC se
    dse = abs(seobs - ST[k,2])/max((abs(ST[k,2]),1e-12))
    dt  = abs(tobs  - ST[k,3])/max((abs(ST[k,3]),1e-12))
    dp  = abs(pp    - ST[k,4])/max((abs(ST[k,4]),1e-12))
    dlb = abs(lb    - ST[k,5])/max((abs(ST[k,5]),1e-12))
    dub = abs(ub    - ST[k,6])/max((abs(ST[k,6]),1e-12))
    worst_se = max((worst_se,dse)); worst_t = max((worst_t,dt)); worst_p = max((worst_p,dp))
    worst_lb = max((worst_lb,dlb)); worst_ub = max((worst_ub,dub))
    flag = (max((dse,dt,dp,dlb,dub)) > 1e-10 ? "  <<<< MISMATCH" : "")
    printf("k=%2.0f %-14s|%-14s  b=%11.7f se_tab=%11.7f se_dump=%11.7f  p=%9.6f/%9.6f%s\n",
      k, rs[k,1], rs[k,2], bobs, ST[k,2], seobs, pp, ST[k,4], flag)
  }
  printf("SGWALK_WORST se=%g t=%g p=%g lb=%g ub=%g  blank=%f\n", worst_se, worst_t, worst_p, worst_lb, worst_ub, nblank)
end
restore
dis "P10_DONE"
