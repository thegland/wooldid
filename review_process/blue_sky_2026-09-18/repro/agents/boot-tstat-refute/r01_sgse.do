clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid

cap prog drop shose
prog def shose
  syntax , tag(string)
  tempname M A B
  mat `M' = e(wdidmainresults)
  cap mat `A' = e(wdidsgresults_0)
  local ra = _rc
  cap mat `B' = e(wdidsgresults_1)
  local rb = _rc
  dis "`tag'_MAIN b=" %14.10f `M'[1,1] " se=" %21x `M'[1,2]
  if `ra'==0 dis "`tag'_SG0  b=" %14.10f `A'[1,1] " se=" %21x `A'[1,2]
  if `rb'==0 dis "`tag'_SG1  b=" %14.10f `B'[1,1] " se=" %21x `B'[1,2]
  if `ra'==0 & `rb'==0 {
    dis "`tag'_RATIO0=" %18.15f `A'[1,2]/`M'[1,2] "  `tag'_RATIO1=" %18.15f `B'[1,2]/`M'[1,2]
    local nr = rowsof(`M')
    local same 0
    forval r = 1/`nr' {
      if `A'[`r',2] < . & `M'[`r',2] < . {
        if reldif(`A'[`r',2], `M'[`r',2]) < 1e-12 local ++same
      }
    }
    dis "`tag'_ROWS_IDENTICAL=" `same' " of " `nr'
  }
  * pairwise contrast rows
  cap mat `A' = e(pwsgs)
  if _rc == 0 {
    local pr = rowsof(`A')
    forval r = 1/`pr' {
      local rn : word `r' of `: rownames `A''
      dis "`tag'_PW[" `r' "] " "`rn'" "  b=" %14.10f `A'[`r',1] " se=" %14.10g `A'[`r',2]
    }
  }
end

*==================================================================
dis "### CASE A: MY unbalanced DGP, UNBALANCED subgroup, cluster(state) == i"
local uSGmode unbal
include "`WD'/gen_unbal.do"
qui count
dis "A_NOBS=" r(N)
tab sg onset, missing
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "A_RC=" _rc
shose , tag(A)
dis "A_NCLUST=" e(N_clust)

*==================================================================
dis "### CASE B: same data, robust (control)"
cap noi wooldid y state t ttre , robust att esprelength(2) espostlength(2) subgroups(sg)
dis "B_RC=" _rc
shose , tag(B)

*==================================================================
dis "### CASE C: same data, cluster(unit) == BELOW i (control)"
cap noi wooldid y state t ttre , cluster(unit) att esprelength(2) espostlength(2) subgroups(sg)
dis "C_RC=" _rc
shose , tag(C)

*==================================================================
dis "### CASE D: same data+options as A but nofastcollapse"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg) nofastcollapse
dis "D_RC=" _rc
shose , tag(D)

*==================================================================
dis "### CASE E: BALANCED subgroup on the same unbalanced panel, cluster(state)"
local uSGmode bal
include "`WD'/gen_unbal.do"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "E_RC=" _rc
shose , tag(E)

*==================================================================
dis "### CASE F: three-level subgroup, unbalanced, cluster(state)"
local uSGmode unbal
include "`WD'/gen_unbal.do"
gen byte sg3 = mod(state,3)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg3)
dis "F_RC=" _rc
tempname FM F0 F1 F2
mat FM = e(wdidmainresults)
cap mat F0 = e(wdidsgresults_0)
cap mat F1 = e(wdidsgresults_1)
cap mat F2 = e(wdidsgresults_2)
dis "F_MAIN b=" %14.10f FM[1,1] " se=" %21x FM[1,2]
dis "F_SG0  b=" %14.10f F0[1,1] " se=" %21x F0[1,2]
dis "F_SG1  b=" %14.10f F1[1,1] " se=" %21x F1[1,2]
dis "F_SG2  b=" %14.10f F2[1,1] " se=" %21x F2[1,2]

*==================================================================
dis "### CASE G: att_it estimand (different weights), unbalanced sg, cluster(state)"
local uSGmode unbal
include "`WD'/gen_unbal.do"
cap noi wooldid y state t ttre , cluster(state) att att_it esprelength(2) espostlength(2) subgroups(sg)
dis "G_RC=" _rc
tempname GM G0 G1
mat GM = e(wdidmainresults)
mat G0 = e(wdidsgresults_0)
mat G1 = e(wdidsgresults_1)
local nr = rowsof(GM)
local rn : rownames GM
forval r = 1/`nr' {
  local nm : word `r' of `rn'
  dis "G_ROW " `r' " `nm'  main_se=" %21x GM[`r',2] "  sg0_se=" %21x G0[`r',2] "  sg1_se=" %21x G1[`r',2]
}
dis "R01_DONE"
