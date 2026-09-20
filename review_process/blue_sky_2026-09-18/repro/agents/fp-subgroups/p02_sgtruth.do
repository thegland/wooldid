* ===========================================================================
* p02 -- FIRST-PRINCIPLES ORACLE for subgroups().
* Noiseless DGP (gysd=0, gtrendsd=0) => wooldid MUST recover the design truth
* exactly.  We rebuild every subgroup estimand, the subgroup contrast and the
* homogeneity test from first principles and compare.
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

* ---------------- data: 60 units, 4 onsets, NO noise ----------------
local gN 60
local gT 9
local gt0 2000
local gonsets "2003 2004 2006"
local gysd 0
local gtsd 0.5
local gaisd 1
local gtaumode hetero
local gtauhet 0.7
include "`GEN'/gen_core.do"

* make the two subgroups have genuinely different average effects:
* rebuild the outcome consistently, ync_new = ync_old - tau_old + tau_new
qui gen double __tauold = truetau
qui replace truetau = truetau * (1 + 0.8*sg2)
qui replace ync = ync - __tauold + truetau
drop __tauold

tab sg2 ttre, m

* ---------------- wooldid run ----------------
cap noi wooldid ync unit t ttre , att att_it att_i att_itime cluster(unit) ///
    subgroups(sg2) esprelength(2) espostlength(2) jointtests
dis "RUN_RC=" _rc

qui gen byte insmp = e(sample)

tempname MAIN SG0 SG1 PW
matrix `MAIN' = e(wdidmainresults)
matrix `SG0'  = e(wdidsgresults_0)
matrix `SG1'  = e(wdidsgresults_1)
matrix `PW'   = e(pwsgs)

dis "=== e(pwsgs) rownames ==="
local pwrn : rownames `PW'
local pwre : roweq `PW'
dis "rn: `pwrn'"
dis "eq: `pwre'"
mat list `PW', format(%21x)

dis "=== allsgssame scalars ==="
foreach s in allsgssame_att allsgssame_att_it allsgssame_att_i allsgssame_att_itime {
  cap dis "`s' = " e(`s')
}
ereturn list

* ---------------- rebuild the truth ----------------
* the four wooldid weights, on the REALIZED sample
cap drop __*
qui gen byte   __post = insmp & ttre < . & t >= ttre
qui gen double __w    = 1
qui egen double __wit = total(cond(insmp,__w,0)), by(unit t)
qui egen double __wip = total(cond(__post,__w,0)), by(unit)
qui egen byte   __tag = tag(unit t) if __post
qui egen double __tpi = total(__tag==1), by(unit)
qui gen double __w_att   = __w
qui gen double __w_it    = __w / __wit
qui gen double __w_i     = __w / __wip
qui gen double __w_itime = __w / (__wit * __tpi)

foreach e in att it i itime {
  local nm `e'
  if "`e'"=="it"    local nm att_it
  if "`e'"=="i"     local nm att_i
  if "`e'"=="itime" local nm att_itime
  qui summarize truetau [aw=__w_`e'] if __post==1, meanonly
  local TRUE_`nm' = r(mean)
  forval v = 0/1 {
    qui summarize truetau [aw=__w_`e'] if __post==1 & sg2==`v', meanonly
    local TRUE_`nm'_s`v' = r(mean)
  }
}

dis ""
dis "==================== POINT-ESTIMATE ORACLE ===================="
local r 1
foreach nm in att att_it att_i att_itime {
  local est_main = `MAIN'[`r',1]
  local est_s0   = `SG0'[`r',1]
  local est_s1   = `SG1'[`r',1]
  dis "`nm'   main: wooldid=" %18.15f `est_main' "  truth=" %18.15f `TRUE_`nm'' "  d=" %12.3e (`est_main'-`TRUE_`nm'')
  dis "`nm'   sg0 : wooldid=" %18.15f `est_s0'   "  truth=" %18.15f `TRUE_`nm'_s0' "  d=" %12.3e (`est_s0'-`TRUE_`nm'_s0')
  dis "`nm'   sg1 : wooldid=" %18.15f `est_s1'   "  truth=" %18.15f `TRUE_`nm'_s1' "  d=" %12.3e (`est_s1'-`TRUE_`nm'_s1')
  local r = `r' + 1
}

dis ""
dis "==================== CONTRAST ORACLE ===================="
* pwsgs rows should be one per (estimand x pair).  Check estimate == sg1 - sg0
forval pr = 1/`=rowsof(`PW')' {
  local rn : word `pr' of `pwrn'
  local re : word `pr' of `pwre'
  dis "pwsgs row `pr'  eq=`re' name=`rn'  est=" %18.15f `PW'[`pr',1] "  se=" %18.15f `PW'[`pr',2]
}
local r 1
foreach nm in att att_it att_i att_itime {
  local d10 = `SG1'[`r',1] - `SG0'[`r',1]
  dis "`nm': sg1-sg0 = " %18.15f `d10'
  local r = `r' + 1
}

dis ""
dis "==================== SE INDEPENDENCE PROBE ===================="
* if wooldid ASSUMED independence, se_c == sqrt(se0^2+se1^2) exactly.
local r 1
foreach nm in att att_it att_i att_itime {
  local s0 = `SG0'[`r',2]
  local s1 = `SG1'[`r',2]
  dis "`nm': se0=" %12.9f `s0' " se1=" %12.9f `s1' " indep=" %12.9f sqrt(`s0'^2+`s1'^2)
  local r = `r' + 1
}
