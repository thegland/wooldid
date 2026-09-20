* ===========================================================================
* p07 -- STRESS: the documented rank-deficiency screen at cluster()==i, weights,
*        missing subgroup values, and a tiny subgroup.
*
*  (A) cluster() AT the level of i.  sthlp "sgdegen" says the contrast covariance
*      goes rank deficient there, that wooldid screens each pairwise row on the
*      implied correlation and, above 1-1e-2, reports SE / t / p / CI and the
*      across-subgroup equality test as ABSENT with a caution.  Does it fire?
*      And is the implied correlation actually near 1?
*  (B) Identity A with analytic weights [aw=w]: only the att weight changes.
*  (C) subgroup with MISSING values on some units: the pooled table must be
*      unchanged (bit-identical %21x) vs a no-subgroups run, and the subgroup
*      aggregation must reproduce att computed on the non-missing subset only.
*  (D) a tiny subgroup (a handful of treated units).
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

local gN 240
local gT 9
local gt0 2000
local gonsets "2003 2004 2006"
local gysd 0.60
local gtsd 0.50
local gaisd 1
local gtaumode hetero
local gtauhet 0.7
local gwtype unit
local gwspread 2
include "`GEN'/gen_core.do"

egen long grp = group(cohort region)
gen byte sgb = (region <= 2)
* sgm: like sgb but a quarter of UNITS carry no subgroup at all
gen byte sgm = sgb
qui replace sgm = . if mod(unit,4)==0
* sgt: level 2 holds only a handful of treated units
gen byte sgt = 1
qui replace sgt = 2 if inlist(grp, 1, 6)

tab sgb, m
tab sgm, m
tab sgt, m

* ===================== (A) cluster() AT the level of i =====================
dis "@@@@@@@@@@ A: cluster(grp) == i @@@@@@@@@@"
cap noi wooldid y grp t ttre , att att_it cluster(grp) subgroups(sgb) ///
    esprelength(2) espostlength(2)
dis "A_RC=" _rc
matrix A_PW = e(pwsgs)
matrix A_S0 = e(wdidsgresults_0)
matrix A_S1 = e(wdidsgresults_1)
scalar A_J  = e(joint_att_allsgssame)
scalar A_DG = e(wdid_sgdegen)
scalar A_DGR = e(wdid_sgdegenr)
scalar A_DGA = e(wdid_sgdegenact)
dis "A joint_att_allsgssame = " A_J
dis "A wdid_sgdegen=" A_DG " sgdegenr=" A_DGR " sgdegenact=" A_DGA
mat list A_PW
dis "A se0(att)=" %20.16f A_S0[1,2] "  se1(att)=" %20.16f A_S1[1,2]
dis "A b0(att) =" %20.16f A_S0[1,1] "  b1(att) =" %20.16f A_S1[1,1]

* ===================== (B) Identity A with [aw=w] =====================
dis "@@@@@@@@@@ B: aweights @@@@@@@@@@"
cap noi wooldid y grp t ttre [aw=w] , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) subgroups(sgb) esprelength(2) espostlength(2)
dis "B_RC=" _rc
qui gen byte insmpB = e(sample)
matrix B_MAIN = e(wdidmainresults)
matrix B_S0 = e(wdidsgresults_0)
matrix B_S1 = e(wdidsgresults_1)

cap drop __*
qui gen byte   __post = insmpB & ttre < . & t >= ttre
qui gen double __w    = w
qui egen double __wit = total(cond(insmpB,__w,0)), by(grp t)
qui egen double __wip = total(cond(__post,__w,0)), by(grp)
qui egen byte   __tag = tag(grp t) if __post
qui egen double __tpi = total(__tag==1), by(grp)
qui gen double __W_att   = __w                    if __post
qui gen double __W_it    = __w / __wit            if __post
qui gen double __W_i     = __w / __wip            if __post
qui gen double __W_itime = __w / (__wit * __tpi)  if __post
foreach e in att it i itime {
  local nm `e'
  if "`e'"=="it"    local nm att_it
  if "`e'"=="i"     local nm att_i
  if "`e'"=="itime" local nm att_itime
  forval v = 0/1 {
    qui summarize __W_`e' if __post==1 & sgb==`v', meanonly
    local WT_`nm'_s`v' = r(sum)
  }
}
dis "--- IDENTITY A under [aw=w] ---"
local r 1
foreach nm in att att_it att_i att_itime {
  local b0 = B_S0[`r',1]
  local b1 = B_S1[`r',1]
  local w0 = `WT_`nm'_s0'
  local w1 = `WT_`nm'_s1'
  local agg = (`w0'*`b0' + `w1'*`b1')/(`w0'+`w1')
  local pooled = B_MAIN[`r',1]
  dis "`nm': W0=" %16.8f `w0' " W1=" %16.8f `w1' "  pooled=" %20.17f `pooled' " rebuilt=" %20.17f `agg' " reldiff=" %12.3e (`agg'-`pooled')/abs(`pooled')
  local r = `r' + 1
}

* ===================== (C) missing subgroup values =====================
dis "@@@@@@@@@@ C: subgroup with missing values @@@@@@@@@@"
cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) subgroups(sgm) esprelength(2) espostlength(2)
dis "C_RC=" _rc
qui gen byte insmpC = e(sample)
matrix C_MAIN = e(wdidmainresults)
matrix C_S0 = e(wdidsgresults_0)
matrix C_S1 = e(wdidsgresults_1)
scalar C_N = e(N)

cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) esprelength(2) espostlength(2)
dis "Cnosg_RC=" _rc
matrix C_NOSG = e(wdidmainresults)
scalar C_N2 = e(N)
dis "C N with sgm=" C_N "  N without subgroups=" C_N2
dis "--- C: main table with vs without subgroups(sgm), %21x ---"
mat list C_MAIN, format(%21x)
mat list C_NOSG, format(%21x)

* aggregation on the NON-MISSING subset only
cap drop __*
qui gen byte   __post = insmpC & ttre < . & t >= ttre
qui gen double __w    = 1
qui egen double __wit = total(cond(insmpC,__w,0)), by(grp t)
qui egen double __wip = total(cond(__post,__w,0)), by(grp)
qui egen byte   __tag = tag(grp t) if __post
qui egen double __tpi = total(__tag==1), by(grp)
qui gen double __W_att   = __w                    if __post
qui gen double __W_it    = __w / __wit            if __post
qui gen double __W_i     = __w / __wip            if __post
qui gen double __W_itime = __w / (__wit * __tpi)  if __post
foreach e in att it i itime {
  local nm `e'
  if "`e'"=="it"    local nm att_it
  if "`e'"=="i"     local nm att_i
  if "`e'"=="itime" local nm att_itime
  forval v = 0/1 {
    qui summarize __W_`e' if __post==1 & sgm==`v', meanonly
    local CWT_`nm'_s`v' = r(sum)
  }
  qui count if __post==1 & sgm >= .
  local CMISS_`nm' = r(N)
}
dis "--- C: does the 2-subgroup aggregation equal the POOLED att? (it must NOT, if some units have missing sg) ---"
local r 1
foreach nm in att att_it att_i att_itime {
  local b0 = C_S0[`r',1]
  local b1 = C_S1[`r',1]
  local w0 = `CWT_`nm'_s0'
  local w1 = `CWT_`nm'_s1'
  local agg = (`w0'*`b0' + `w1'*`b1')/(`w0'+`w1')
  dis "`nm': W0=" %14.4f `w0' " W1=" %14.4f `w1' "  pooled=" %20.17f C_MAIN[`r',1] " sgagg=" %20.17f `agg' " reldiff=" %12.3e (`agg'-C_MAIN[`r',1])/abs(C_MAIN[`r',1])
  local r = `r' + 1
}
dis "post rows with missing sgm: `CMISS_att'"

* ===================== (D) tiny subgroup =====================
dis "@@@@@@@@@@ D: tiny subgroup @@@@@@@@@@"
cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) subgroups(sgt) esprelength(2) espostlength(2)
dis "D_RC=" _rc
cap mat list e(pwsgs)
cap mat list e(wdidsgresults_2)
dis "D joint_att_allsgssame = " e(joint_att_allsgssame)
