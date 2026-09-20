*==============================================================================
* p1_level.do -- EXACT-IDENTITY probe of level() on every analytic surface.
* Checks:  lb == b - c*se  and  ub == b + c*se  with c = invttail(df,(1-L/100)/2)
*          (or invnormal for z rows), column labels lb_<L>ci, e(level)==L,
*          and c(level) restoration.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---- helper: verify the CI identity inside a results matrix -------------------
cap prog drop ckci
prog def ckci
  syntax , mat(string) lvl(string) [df(string) tag(string)]
  if "`df'" == "" local df .
  local cn : colnames `mat'
  local cf : colfullnames `mat'
  local rn : rownames `mat'
  * locate columns
  local cE = 0
  local cS = 0
  local cL = 0
  local cU = 0
  local ci = 0
  local lbname ""
  foreach c of local cn {
    local ++ci
    if "`c'" == "estimate" local cE = `ci'
    if "`c'" == "se" | "`c'" == "se_analytic" local cS = `ci'
    if substr("`c'",1,3) == "lb_" {
      local cL = `ci'
      local lbname "`c'"
    }
    if substr("`c'",1,3) == "ub_" local cU = `ci'
  }
  dis as text "[`tag'] colnames: `cn'"
  if `cE'==0 | `cS'==0 | `cL'==0 | `cU'==0 {
    dis as error "[`tag'] MISSING one of estimate/se/lb/ub -- skipping identity check"
    exit
  }
  * expected label
  local ltxt = cond(mod(`lvl',1)==0, strofreal(`lvl'), subinstr(strofreal(`lvl',"%9.0g"), ".", "p", .))
  if "`lbname'" != "lb_`ltxt'ci" dis as error "[`tag'] LABEL MISMATCH: got `lbname', expected lb_`ltxt'ci"
  else dis as text "[`tag'] label OK (`lbname')"
  * critical value
  local a = (1-`lvl'/100)/2
  if `df' < . local crit = invttail(`df',`a')
  else        local crit = invnormal(1-`a')
  local nr = rowsof(`mat')
  local worst = 0
  local worstr ""
  forvalues r = 1/`nr' {
    local b  = `mat'[`r',`cE']
    local s  = `mat'[`r',`cS']
    local lo = `mat'[`r',`cL']
    local hi = `mat'[`r',`cU']
    if `b' >= . | `s' >= . | `lo' >= . continue
    local elo = `b' - `crit'*`s'
    local ehi = `b' + `crit'*`s'
    local d1 = abs(`lo'-`elo')/max(abs(`elo'),1e-12)
    local d2 = abs(`hi'-`ehi')/max(abs(`ehi'),1e-12)
    local dd = max(`d1',`d2')
    if `dd' > `worst' {
      local worst = `dd'
      local rw : word `r' of `rn'
      local worstr "`rw'"
    }
  }
  dis as result "[`tag'] crit=" %12.10f `crit' "  worst rel dev = " %12.4e `worst' "  (row `worstr')"
  if `worst' > 1e-8 dis as error "[`tag'] ****** CI IDENTITY FAILS at level `lvl' ******"
end

* ---- data --------------------------------------------------------------------
local dN 150
local dT 12
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay/mkdata.do"
save "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay/d1.dta", replace

dis "SESSION_LEVEL_START=" c(level)

* ============================================================================
* A. level(80), robust, analytic -- main + ES + subgroups + gteffects
* ============================================================================
cap noi wooldid y i t ttre, att att_it att_i att_itime robust espre(3) espost(3) ///
    jointtests gteffects subgroups(sg) level(80)
dis "A_RC=" _rc
dis "A_level=" e(level)
dis "A_df_r=" e(df_r)
dis "A_sessionlevel_after=" c(level)
mat A_main = e(wdidmainresults)
mat list A_main, format(%21x)
ckci , mat(A_main) lvl(80) df(`=e(df_r)') tag(A-main)
cap confirm matrix e(wdidsgresults_1)
if _rc == 0 {
  mat A_sg1 = e(wdidsgresults_1)
  ckci , mat(A_sg1) lvl(80) df(`=e(df_r)') tag(A-sg1)
}
cap confirm matrix e(pwsgs)
if _rc == 0 {
  mat A_pw = e(pwsgs)
  ckci , mat(A_pw) lvl(80) df(`=e(df_r)') tag(A-pwsgs)
}
cap confirm matrix e(gteffects)
if _rc == 0 {
  mat A_gt = e(gteffects)
  local gtc : colnames A_gt
  dis as text "GTEFFECTS colnames: `gtc'"
  mat list A_gt, format(%10.6f) nohalf
}

* ============================================================================
* B. same call at default (95) -- for comparison of which numbers move
* ============================================================================
cap noi wooldid y i t ttre, att att_it att_i att_itime robust espre(3) espost(3) ///
    jointtests gteffects subgroups(sg)
dis "B_RC=" _rc
dis "B_level=" e(level)
mat B_main = e(wdidmainresults)
ckci , mat(B_main) lvl(95) df(`=e(df_r)') tag(B-main)
cap confirm matrix e(gteffects)
if _rc == 0 {
  mat B_gt = e(gteffects)
  mat list B_gt, format(%10.6f) nohalf
}

* ---- do the gteffects cells actually MOVE with level()? ---------------------
cap noi {
  local same = 1
  local nc = colsof(A_gt)
  local nr = rowsof(A_gt)
  forvalues r = 1/`nr' {
    forvalues c = 1/`nc' {
      if A_gt[`r',`c'] != B_gt[`r',`c'] local same = 0
    }
  }
  dis "GTEFFECTS_IDENTICAL_ACROSS_LEVELS=" `same'
}
dis "GTCMP_RC=" _rc

* ============================================================================
* C. non-integer level
* ============================================================================
cap noi wooldid y i t ttre, att robust espre(2) espost(2) level(92.5)
dis "C_RC=" _rc
dis "C_level=" e(level)
dis "C_sessionlevel_after=" c(level)
cap noi mat C_main = e(wdidmainresults)
cap noi ckci , mat(C_main) lvl(92.5) df(`=e(df_r)') tag(C-main)

* ============================================================================
* D. level with more decimals than set level allows
* ============================================================================
cap noi wooldid y i t ttre, att robust espre(2) espost(2) level(92.567)
dis "D_RC=" _rc
dis "D_level=" e(level)
dis "D_sessionlevel_after=" c(level)
cap noi mat D_main = e(wdidmainresults)
cap noi mat list D_main, format(%21x)

* ============================================================================
* E. out-of-range levels
* ============================================================================
foreach L in 5 100 -1 0.5 {
  dis as text "---- trying level(`L') ----"
  cap noi wooldid y i t ttre, att robust level(`L')
  dis "E_`L'_RC=" _rc
  dis "E_`L'_sessionlevel_after=" c(level)
  dis "E_`L'_LEVEL0_global=[$WOOLDID_LEVEL0]"
}

* ============================================================================
* F. restoration after a FAILED run (failure AFTER set level fires)
* ============================================================================
set level 90
dis "F_before=" c(level)
cap noi wooldid y i t ttre if t > 1e9, att robust level(80)
dis "F_failrc=" _rc
dis "F_after=" c(level)
dis "F_LEVEL0_global=[$WOOLDID_LEVEL0]"

* failure caused by a nonexistent control variable (parse stage)
cap noi wooldid y i t ttre, att robust controls(nosuchvar) level(80)
dis "F2_failrc=" _rc
dis "F2_after=" c(level)

* ============================================================================
* G. back-to-back different levels; then a call with no level()
* ============================================================================
set level 90
cap noi wooldid y i t ttre, att robust level(80)
dis "G1_level=" e(level) "  session=" c(level)
cap noi wooldid y i t ttre, att robust level(99)
dis "G2_level=" e(level) "  session=" c(level)
cap noi wooldid y i t ttre, att robust
dis "G3_level=" e(level) "  session=" c(level)
mat G3 = e(wdidmainresults)
local g3c : colnames G3
dis "G3_colnames=`g3c'"
ckci , mat(G3) lvl(90) df(`=e(df_r)') tag(G3-sessionlevel90)
set level 95

dis "PROBE1_DONE"
