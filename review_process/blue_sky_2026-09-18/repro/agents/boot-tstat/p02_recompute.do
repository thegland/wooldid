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
dis "NOBS=" _N
qui levelsof state, local(SS)
dis "NSTATES=" `: word count `SS''

timer clear 1
timer on 1
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(199) boottstat rseed(4242) bootsave(d2) replace att esprelength(2) espostlength(2) dots(0)
dis "P02_RC=" _rc
timer off 1
timer list 1

mat MR  = e(wdidmainresults)
mat TS  = e(boot_tstud)
mat ASE = e(boot_ase_obs)
mat BB  = e(b)
local LVL = e(boot_level)
local DONE = e(boot_completed)
local BDF  = e(boot_df)
dis "LEVEL=`LVL' DONE=`DONE' BOOTDF=`BDF'"
dis "BOOT_STAT=" "`e(boot_stat)'"
dis "NCLUS_e=" e(N_clust)
dis "DFR=" e(df_r)
mat list MR, format(%21x)
mat list TS, format(%21x)

* ---------- recompute from the dump ----------
preserve
use d2, clear
dis "DUMPOBS=" _N
ds
restore

mata:
  MR  = st_matrix("MR")
  TS  = st_matrix("TS")
  ASE = st_matrix("ASE")
  BB  = st_matrix("BB")
  cs  = st_matrixcolstripe("BB")
  lvl = strtoreal(st_local("LVL"))
  al  = (100-lvl)/200
  K   = cols(BB)
  nm  = rows(MR)
  printf("K=%f  nmainrows=%f\n", K, nm)
  st_local("dumpfile","d2")
end

preserve
use d2, clear
mata:
  MR  = st_matrix("MR"); TS = st_matrix("TS"); ASE = st_matrix("ASE"); BB = st_matrix("BB")
  cs  = st_matrixcolstripe("BB")
  lvl = strtoreal(st_local("LVL")); al = (100-lvl)/200
  K = cols(BB)
  maxrd_p = 0; maxrd_lb = 0; maxrd_ub = 0; maxrd_t = 0; maxrd_n = 0
  maxrd_b0 = 0; maxrd_se0 = 0
  nbad = 0
  for (k = 1; k <= K; k++) {
    vn = strtoname(cs[k,2])
    sn = substr("se_" + vn, 1, 32)
    if (_st_varindex(vn) >= .) {
      printf("MISSINGVAR %s\n", vn)
      nbad++
      continue
    }
    if (_st_varindex(sn) >= .) {
      printf("MISSINGSEVAR %s\n", sn)
      nbad++
      continue
    }
    B = st_data(., (vn, sn, "_rep"))
    o = select(B, B[.,3] :== 0)
    bobs = o[1,1]; seobs = o[1,2]
    // invariant: dump _rep==0 coefficient == e(b)
    d = abs(bobs - BB[1,k]) / max((abs(BB[1,k]), 1e-12))
    if (d > maxrd_b0) maxrd_b0 = d
    d = abs(seobs - ASE[1,k]) / max((abs(ASE[1,k]), 1e-12))
    if (d > maxrd_se0) maxrd_se0 = d
    R = select(B, B[.,3] :> 0)
    ok = (R[.,1] :< .) :& (R[.,2] :< .) :& (R[.,2] :> 0)
    if (sum(ok) < 2) continue
    Rk = select(R, ok)
    tr = sort((Rk[.,1] :- bobs) :/ Rk[.,2], 1)
    n = rows(tr)
    tobs = bobs/seobs
    // order-statistic rule (wooldid_bootqs)
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
    // compare against TS and MR row k (no subgroups -> row k)
    d = abs(n - TS[1,k]); if (d > maxrd_n) maxrd_n = d
    d = abs(tobs - TS[2,k])/max((abs(TS[2,k]),1e-12)); if (d > maxrd_t) maxrd_t = d
    d = abs(pp - TS[3,k])/max((abs(TS[3,k]),1e-12)); if (d > maxrd_p) maxrd_p = d
    d = abs(lb - MR[k,5])/max((abs(MR[k,5]),1e-12)); if (d > maxrd_lb) maxrd_lb = d
    d = abs(ub - MR[k,6])/max((abs(MR[k,6]),1e-12)); if (d > maxrd_ub) maxrd_ub = d
    printf("col %f %s : n=%f/%f p=%9.6f/%9.6f lb=%12.8f/%12.8f ub=%12.8f/%12.8f\n", k, cs[k,2], n, TS[1,k], pp, TS[3,k], lb, MR[k,5], ub, MR[k,6])
  }
  printf("RECOMP_MAXRD n=%g t=%g p=%g lb=%g ub=%g b0=%g se0=%g nbad=%f\n", maxrd_n, maxrd_t, maxrd_p, maxrd_lb, maxrd_ub, maxrd_b0, maxrd_se0, nbad)
end

* ---------- shape invariants on the dump ----------
qui ds w_*
local WV `r(varlist)'
dis "NWVARS=" `: word count `WV''
qui egen double __wt = rowtotal(`WV') if _rep > 0
qui sum __wt, detail
dis "WROWTOT_min=" r(min) " max=" r(max)
qui ds se_*
dis "NSEVARS=" `: word count `r(varlist)''

* ---------- SE-inflation diagnostic: mean(se_r)/se_obs per column ----------
foreach v of varlist se_* {
    qui sum `v' if _rep > 0
    local mr = r(mean)
    qui sum `v' if _rep == 0
    local mo = r(mean)
    if `mo' > 0 dis "SERATIO `v' = " %8.5f `mr'/`mo'
}
restore
dis "P02_DONE"
