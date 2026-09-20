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
gen double wt = 0.5 + runiform()
save legs20, replace

mata:
void serat(string scalar tg)
{
  real matrix S, A
  real colvector v
  real scalar k
  S = st_matrix("__S"); A = st_matrix("__A")
  for (k=1; k<=cols(S); k++) {
    v = select(S[.,k], S[.,k] :< .)
    if (rows(v) < 2) continue
    printf("%s col %f n=%f ratio_mean=%9.6f ratio_rms=%9.6f\n", tg, k, rows(v), mean(v)/A[1,k], sqrt(mean(v:^2))/A[1,k])
  }
}
end
cap prog drop seratio
prog def seratio
  syntax , tag(string)
  cap confirm matrix e(bootdist_se)
  if _rc != 0 {
    dis "`tag'_NOSEDIST"
    exit
  }
  mat __S = e(bootdist_se)
  mat __A = e(boot_ase_obs)
  mata: serat("`tag'")
end

*=== L1: boottstat, NO cluster, no weight (the fweight fix path) ===
dis "=== L1 unclustered boottstat (fweight leg) ==="
cap noi wooldid y state t ttre , robust bootstrap(standard) reps(199) boottstat rseed(5) att esprelength(2) espostlength(2) dots(0)
dis "L1_RC=" _rc
dis "L1_BOOTDF=" e(boot_df) " DFR=" e(df_r)
seratio , tag(L1)

*=== L2: boottstat, NO cluster, WITH a user aweight (the phys leg) ===
dis "=== L2 unclustered boottstat + aweight (phys leg) ==="
cap noi wooldid y state t ttre [aweight=wt], robust bootstrap(standard) reps(199) boottstat rseed(5) att esprelength(2) espostlength(2) dots(0)
dis "L2_RC=" _rc
seratio , tag(L2)

*=== L3: boottstat + bootstrap(bayes) + cluster  (no reid leg exists for bayes) ===
dis "=== L3 boottstat + bayes + cluster ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(bayes) reps(199) boottstat rseed(5) att esprelength(2) espostlength(2) dots(0)
dis "L3_RC=" _rc
dis "L3_METHOD=" "`e(boot_method)'" " ALPHA=" "`e(boot_alpha)'"
seratio , tag(L3)

*=== L4: boottstat + standard + cluster (reference) ===
dis "=== L4 boottstat + standard + cluster ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(199) boottstat rseed(5) att esprelength(2) espostlength(2) dots(0)
dis "L4_RC=" _rc
seratio , tag(L4)

*=== L5: bootdf() parse edge cases ===
dis "=== L5 bootdf parse ==="
foreach d in 0 -5 2.5 1e9 clusters normal {
  cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(25) bootdf(`d') rseed(5) att dots(0)
  dis "L5_`d'_RC=" _rc " BOOTDF=" e(boot_df) " DFR=" e(df_r)
}
dis "P08_DONE"
