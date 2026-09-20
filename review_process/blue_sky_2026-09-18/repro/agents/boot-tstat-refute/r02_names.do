clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909090
set sortseed 909090
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid
local uG 30
local uT 8
local uRagged 0
local uSGmode bal
include "`WD'/gen_unbal.do"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   bootstrap(standard) reps(99) boottstat simultaneousbands
dis "N_RC=" _rc
dis "--- e(b) colfullnames ---"
mat B = e(b)
local cn : colfullnames B
foreach c of local cn {
  dis "BCOL: `c'"
}
dis "--- e(wdidmainresults) rowfullnames ---"
mat W = e(wdidmainresults)
local rn : rowfullnames W
foreach r of local rn {
  dis "WROW: `r'"
}
dis "--- e(esband_att) ---"
cap mat list e(esband_att), format(%12.8f)
dis "ESBAND_RC=" _rc
dis "esband_crit_att=" %21x e(esband_crit_att)
dis "esband_n_att=" e(esband_n_att)
dis "esband_engine=" "`e(esband_engine)'"
dis "--- e() list ---"
qui ereturn list
dis "--- boot_tstud colfullnames ---"
cap mat T = e(boot_tstud)
if _rc==0 {
  local tn : colfullnames T
  foreach c of local tn {
    dis "TCOL: `c'"
  }
}
dis "--- boot_ase_obs ---"
cap mat A = e(boot_ase_obs)
dis "ASE_RC=" _rc
dis "R02_DONE"
