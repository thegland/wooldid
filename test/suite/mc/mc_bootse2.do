* mc_bootse2.do - F058 extended bootstrap calibration (author ruling 2026-09-06).
* Feature-scoped; NOT part of run_full.sh. Companion to mc_bootse.do, which measured
* only the Wald (normal-theory) family on one balanced design. This driver measures
* the CI families that mc_bootse.do's own closing note called a zero-cost extension -
* percentile, bias-corrected, BCa, percentile-t - plus the p-value sizes, the two
* joint-test families and the att_it re-identification leg, on mc_bootse.do's design A
* and on a NEW unequal-cluster design B (12 states of 1..12 units).
*
* usage: <stata> -e do mc_bootse2.do <ado> <cfg> <first> <last> <outdir> [reps]
*   cfg in: AN AA AT BA BN BI BJ BT   (see the config table below)
* One CSV row PER ARM RUN, appended as it completes, so a run may be killed at any
* moment and analysed from the rows that landed, or resumed by raising <last>.
clear all
set more off
set linesize 200
set graphics off

local ADO `1'
if `"`ADO'"' == "" local ADO ../../../wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
qui do `"`ADO'"'
local HERE = subinstr(`"`ADO'"', "wooldidinstall/wooldid.ado", "test/suite/mc/", 1)
run `"`HERE'mc_bootse2_dgp.doh"'

local CFG `2'
if "`CFG'" == "" local CFG AN
local FIRST `3'
if "`FIRST'" == "" local FIRST 1
local LAST `4'
if "`LAST'" == "" local LAST 200
local OUTDIR `5'
if `"`OUTDIR'"' == "" local OUTDIR .
local REPS `6'
if "`REPS'" == "" local REPS 199

*----------------------------------------------------------------- config table
* DESIGN A|B ; DGP alt|null ; ARMS ; the ES window flag ; the estimand
local ES esprelength(2) espostlength(1) jointtests
if "`CFG'" == "AN" {
  local DESIGN A
  local DGP    null
  local ARMS   plain std bayes
  local ESOPT  `ES'
}
else if "`CFG'" == "AA" {
  local DESIGN A
  local DGP    alt
  local ARMS   plain std bayes
  local ESOPT  `ES'
}
else if "`CFG'" == "AT" {
  local DESIGN A
  local DGP    null
  local ARMS   plain bts
  local ESOPT  `ES'
}
else if "`CFG'" == "BA" {
  local DESIGN B
  local DGP    alt
  local ARMS   plain stdbca bayes
  local ESOPT
}
else if "`CFG'" == "BN" {
  local DESIGN B
  local DGP    null
  local ARMS   plain stdbca bayes
  local ESOPT
}
else if "`CFG'" == "BI" {
  local DESIGN B
  local DGP    alt
  local ARMS   plainit attit
  local ESOPT
}
else if "`CFG'" == "BJ" {
  local DESIGN B
  local DGP    alt
  local ARMS   plain std
  local ESOPT  `ES'
}
else if "`CFG'" == "BT" {
  local DESIGN B
  local DGP    null
  local ARMS   plain bts
  local ESOPT  `ES'
}
else {
  di as error "MCB2 unknown cfg `CFG'"
  exit 198
}

* the exact truth: design A's effect is the constant tau on every treated post cell,
* so the att is tau under any window; design B's is the unit-weighted mean of tau_s
* over treated units, 0.1*sum_{s odd}s^2 / sum_{s odd}s = 0.1*286/36, identical for
* att and att_it because every treated unit carries the same four post periods
if "`DESIGN'" == "A" {
  local MULT = cond("`DGP'" == "alt", 2, 0)
  local TRUTH = `MULT'
  local CVAR cl
  local SEEDBASE = 7400000
}
else {
  local MULT = cond("`DGP'" == "alt", 1, 0)
  local TRUTH = `MULT'*0.1*286/36
  local CVAR state
  local SEEDBASE = 7600000
}
if "`DGP'" == "null" local SEEDBASE = `SEEDBASE' + 100000

local OUTFILE `"`OUTDIR'/mcb2_`CFG'.csv"'
local VARS cfg,design,dgp,draw,dgpseed,estimand,truth,arm,rc,secs,nobs,b,se,analse,df,nv,ncomp,nfail,main_lb,main_ub,pct_lb,pct_ub,bc_lb,bc_ub,bca_lb,bca_ub,wald_lb,wald_ub,p_pct,p_bc,p_bca,jt_p,jt_fam
di "MCB2_CONFIG cfg=`CFG' design=`DESIGN' dgp=`DGP' arms=`ARMS' es=`ESOPT' reps=`REPS' truth=`TRUTH' draws=`FIRST'..`LAST'"
di "MCB2_OUTFILE `OUTFILE'"

*--------------------------------------------------------- row writer
* reads the surfaces of the run that just finished and appends one CSV line. Not
* eclass and it never ereturns, so e() survives the call untouched.
cap prog drop mcb2row
prog def mcb2row
  syntax , file(string) cfg(string) design(string) dgp(string) draw(int) seed(passthru) ///
           estimand(string) truth(string) arm(string) rc(int) secs(string) nobs(int) [BOOTstrap]
  local dgpseed `=subinstr("`seed'", "seed(", "", 1)'
  local dgpseed = subinstr("`dgpseed'", ")", "", .)
  foreach v in b se analse df nv ncomp nfail main_lb main_ub pct_lb pct_ub bc_lb bc_ub ///
               bca_lb bca_ub wald_lb wald_ub p_pct p_bc p_bca jt_p {
    local `v' = .
  }
  local jt_fam
  if `rc' == 0 {
    tempname MR WM CP CB CA PP PB PA VA
    cap matrix `MR' = e(wdidmainresults)
    if _rc == 0 {
      local b      = el(`MR',1,1)
      local se     = el(`MR',1,2)
      local main_lb = el(`MR',1,5)
      local main_ub = el(`MR',1,6)
    }
    if "`bootstrap'" != "" {
      cap matrix `WM' = e(boot_wald_mainresults)
      if _rc == 0 {
        local b       = el(`WM',1,1)
        local se      = el(`WM',1,2)
        local wald_lb = el(`WM',1,5)
        local wald_ub = el(`WM',1,6)
      }
      cap matrix `VA' = e(V_analytic)
      if _rc == 0 local analse = sqrt(el(`VA',1,1))
      local df    = e(boot_df)
      local nv    = e(boot_N_V)
      local ncomp = e(boot_completed)
      local nfail = e(boot_failed)
      foreach pair in CP:boot_ci_percentile CB:boot_ci_bc CA:boot_ci_bca {
        local nm = substr("`pair'", 1, 2)
        local em = substr("`pair'", 4, .)
        cap matrix ``nm'' = e(`em')
        if _rc continue
        local k = colnumb(``nm'', "`estimand'")
        if `k' >= . continue
        local tgt = cond("`nm'"=="CP", "pct", cond("`nm'"=="CB", "bc", "bca"))
        local `tgt'_lb = el(``nm'',1,`k')
        local `tgt'_ub = el(``nm'',2,`k')
      }
      foreach pair in PP:boot_p_percentile PB:boot_p_bc PA:boot_p_bca {
        local nm = substr("`pair'", 1, 2)
        local em = substr("`pair'", 4, .)
        cap matrix ``nm'' = e(`em')
        if _rc continue
        local k = colnumb(``nm'', "`estimand'")
        if `k' >= . continue
        local tgt = cond("`nm'"=="PP", "p_pct", cond("`nm'"=="PB", "p_bc", "p_bca"))
        local `tgt' = el(``nm'',1,`k')
      }
      if "`e(joint_horiz)'" != "" local jt_p = e(joint_horiz)
      local jt_fam `e(boot_jt_family)'
    }
  }
  local row "`cfg',`design',`dgp',`draw',`dgpseed',`estimand',`truth',`arm',`rc',`secs',`nobs'"
  foreach v in b se analse df nv ncomp nfail main_lb main_ub pct_lb pct_ub bc_lb bc_ub ///
               bca_lb bca_ub wald_lb wald_ub p_pct p_bc p_bca jt_p {
    local x = strtrim(string(``v'', "%18.12g"))
    if "`x'" == "" | "`x'" == "." local x .
    local row "`row',`x'"
  }
  local row "`row',`jt_fam'"
  tempname fh
  file open `fh' using `"`file'"', write append
  file write `fh' `"`row'"' _n
  file close `fh'
end

tempname fh0
cap confirm file `"`OUTFILE'"'
if _rc {
  file open `fh0' using `"`OUTFILE'"', write replace
  file write `fh0' "`VARS'" _n
  file close `fh0'
}

*--------------------------------------------------------------------- the loop
forval m = `FIRST'/`LAST' {
  local dgpseed = `SEEDBASE' + `m'
  if "`DESIGN'" == "A" dgpA `dgpseed' `MULT'
  else                 dgpB `dgpseed' `MULT'
  qui count
  local nobs = r(N)
  local rs = 20000000 + `SEEDBASE' + `m'

  foreach a of local ARMS {
    local est main_att
    if "`a'" == "attit" | "`a'" == "plainit" local est main_att_it
    timer clear 8
    timer on 8
    if "`a'" == "plain" {
      cap qui wooldid y unit t ttre, att cluster(`CVAR') `ESOPT' safety(off)
      local rc = _rc
      local bootflag
    }
    else if "`a'" == "plainit" {
      cap qui wooldid y unit t ttre, att_it cluster(`CVAR') safety(off)
      local rc = _rc
      local bootflag
    }
    else if "`a'" == "std" {
      cap qui wooldid y unit t ttre, att cluster(`CVAR') `ESOPT' ///
          bootstrap(standard) reps(`REPS') rseed(`rs') dots(0) safety(off)
      local rc = _rc
      local bootflag bootstrap
    }
    else if "`a'" == "stdbca" {
      cap qui wooldid y unit t ttre, att cluster(`CVAR') `ESOPT' ///
          bootstrap(standard) bootdisplay(bca) reps(`REPS') rseed(`rs') dots(0) safety(off)
      local rc = _rc
      local bootflag bootstrap
    }
    else if "`a'" == "bayes" {
      cap qui wooldid y unit t ttre, att cluster(`CVAR') `ESOPT' ///
          bootstrap(bayes) reps(`REPS') rseed(`=`rs'+1000000') dots(0) safety(off)
      local rc = _rc
      local bootflag bootstrap
    }
    else if "`a'" == "attit" {
      cap qui wooldid y unit t ttre, att_it cluster(`CVAR') ///
          bootstrap(standard) reps(`REPS') rseed(`=`rs'+2000000') dots(0) safety(off)
      local rc = _rc
      local bootflag bootstrap
    }
    else if "`a'" == "bts" {
      cap qui wooldid y unit t ttre, att cluster(`CVAR') `ESOPT' ///
          bootstrap(standard) boottstat reps(`REPS') rseed(`=`rs'+3000000') dots(0) safety(off)
      local rc = _rc
      local bootflag bootstrap
    }
    timer off 8
    qui timer list 8
    local secs = r(t8)
    mcb2row , file(`"`OUTFILE'"') cfg(`CFG') design(`DESIGN') dgp(`DGP') draw(`m') ///
        seed(`dgpseed') estimand(`est') truth(`TRUTH') arm(`a') rc(`rc') secs(`secs') ///
        nobs(`nobs') `bootflag'
    di "MCB2_ARM `CFG' m=`m' arm=`a' rc=`rc' secs=" %6.1f `secs'
  }
}
di "MCB2_DONE cfg=`CFG' first=`FIRST' last=`LAST'"
