* ===========================================================================
* p03 -- INDEPENDENT REIMPLEMENTATION of the subgroup machinery on NOISY data.
*
* Identity A  (aggregation): pooled estimand == subgroup-weight-weighted mean
*                            of the subgroup estimands.  Must hold EXACTLY,
*                            because subgroups() does not touch the regression.
* Identity C  (contrast)   : pwsgs estimate == difference of the two levels.
* Identity D  (2 subgroups): the "all SGs the same" joint p == the single
*                            pairwise p.  One restriction either way.
* Identity S  (covariance) : is se_contrast == sqrt(se0^2+se1^2)?  If yes the
*                            off-diagonal is being ignored.  Also: does the
*                            implied correlation lie in [-1,1]?
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

* ---------------- data: 240 units, i = cohort x region (coarser than unit) ----
local gN 240
local gT 9
local gt0 2000
local gonsets "2003 2004 2006"
local gysd 0.60
local gtsd 0.50
local gaisd 1
local gtaumode hetero
local gtauhet 0.7
include "`GEN'/gen_core.do"

* i = cohort-by-region: 20 levels, several units each.  cluster(unit) is then
* strictly FINER than i, which is the design the help recommends.
egen long grp = group(cohort region)
* subgroup: constant within grp by construction (a function of region)
gen byte sg2b = (region <= 2)
label define SGL 0 "Rest" 1 "Core"
label values sg2b SGL

tab grp sg2b, m

* ---------------- run ----------------
cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) subgroups(sg2b) esprelength(2) espostlength(2)
dis "RUN_RC=" _rc

qui gen byte insmp = e(sample)
local dfr = e(df_r)
dis "DFR=`dfr'"

tempname MAIN SG0 SG1 PW
matrix `MAIN' = e(wdidmainresults)
matrix `SG0'  = e(wdidsgresults_0)
matrix `SG1'  = e(wdidsgresults_1)
matrix `PW'   = e(pwsgs)

foreach s in joint_att_allsgssame joint_att_it_allsgssame joint_att_i_allsgssame joint_att_itime_allsgssame {
  cap local J_`s' = e(`s')
  dis "e(`s') = " e(`s')
}
ereturn list

* ---------------- Identity A: rebuild the aggregation weights ----------------
cap drop __*
qui gen byte   __post = insmp & ttre < . & t >= ttre
qui gen double __w    = 1
qui egen double __wit = total(cond(insmp,__w,0)), by(grp t)
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
    qui summarize __W_`e' if __post==1 & sg2b==`v', meanonly
    local WT_`nm'_s`v' = r(sum)
  }
}

dis ""
dis "============ IDENTITY A: pooled == weight-weighted mean of sg ============"
local r 1
foreach nm in att att_it att_i att_itime {
  local b0 = `SG0'[`r',1]
  local b1 = `SG1'[`r',1]
  local w0 = `WT_`nm'_s0'
  local w1 = `WT_`nm'_s1'
  local agg = (`w0'*`b0' + `w1'*`b1')/(`w0'+`w1')
  local pooled = `MAIN'[`r',1]
  dis "`nm': W0=" %14.6f `w0' " W1=" %14.6f `w1'
  dis "   pooled=" %20.17f `pooled' "  rebuilt=" %20.17f `agg' "  reldiff=" %12.3e (`agg'-`pooled')/abs(`pooled')
  local r = `r' + 1
}

dis ""
dis "============ IDENTITY C + S: contrast and its SE ============"
local pwrn : rownames `PW'
local pwre : roweq `PW'
mat list `PW'
local r 1
foreach nm in att att_it att_i att_itime {
  * find the pwsgs row with eq == `nm'
  local hit 0
  forval pr = 1/`=rowsof(`PW')' {
    local re : word `pr' of `pwre'
    if "`re'" == "`nm'" & `hit' == 0 {
      local hit `pr'
    }
  }
  local b0 = `SG0'[`r',1]
  local b1 = `SG1'[`r',1]
  local s0 = `SG0'[`r',2]
  local s1 = `SG1'[`r',2]
  local bc = `PW'[`hit',1]
  local sc = `PW'[`hit',2]
  local indep = sqrt(`s0'^2+`s1'^2)
  local cov = (`s0'^2 + `s1'^2 - `sc'^2)/2
  local rho = `cov'/(`s0'*`s1')
  dis "`nm' (pwsgs row `hit'):"
  dis "   b1-b0=" %20.17f (`b1'-`b0') "  pw_est=" %20.17f `bc' "  d=" %12.3e (`bc'-(`b1'-`b0'))
  dis "   se0=" %14.10f `s0' " se1=" %14.10f `s1' " se_c=" %14.10f `sc' " indep=" %14.10f `indep'
  dis "   implied rho=" %14.10f `rho'
  local r = `r' + 1
}

dis ""
dis "============ IDENTITY D: 2-subgroup joint p == pairwise p ============"
local r 1
foreach nm in att att_it att_i att_itime {
  local hit 0
  forval pr = 1/`=rowsof(`PW')' {
    local re : word `pr' of `pwre'
    if "`re'" == "`nm'" & `hit' == 0 local hit `pr'
  }
  local pwp = `PW'[`hit',4]
  local jp  = `J_joint_`nm'_allsgssame'
  dis "`nm': pairwise p=" %20.17f `pwp' "  joint p=" %20.17f `jp' "  d=" %12.3e (`jp'-`pwp')
  local r = `r' + 1
}
