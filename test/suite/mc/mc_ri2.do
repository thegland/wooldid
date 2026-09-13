* mc_ri2.do - F058 extended permutation calibration (author ruling 2026-09-06).
* One shard runs FOUR phases in priority order, each under its own wall-clock budget,
* appending ONE CSV ROW PER (rep, arm) with open/write/close so a shard killed at the
* deadline still leaves every completed rep on disk.
*
*   D1  tau = 0  sharp null, default (equal-tailed) p; jointtests; the step-down record
*                (permstepdown is passed on the main arm: opt-in since PS-45/7a158ea, without
*                it e(perm_stepdown) is not posted and tobs/praw/pss/psd are all '.').
*                Every VERIFY-th rep also runs the package's permabsp / permrawp /
*                permside(left) arms on the SAME rseed (identical draws), so the
*                conventions computed in-driver from e(permdist_b) can be checked
*                against the package bit-for-bit.
*   D2  tau = 1  planted constant effect, default arm: shift-interval coverage of tau.
*   D3  tau = 3  planted constant effect, default arm: the F035 width growth.
*   D4  tau = 1  permtstat arm: does studentizing remove the width growth?
*
* DGP = mc_ri.do's R1: 12 clusters i, onset labels {5 x3, 7 x3, . x6} shuffled uniformly,
* 10 periods, y = a_i + .1 t + e_it (+ tau*post). Reported window esprelength(2)
* espostlength(2), so the table carries pre2 pre1 post0 post1 alongside main/pretreatment.
*
* args: 1=ado 2=shard 3=outdir 4=bD1 5=bD2 6=bD3 7=bD4 8=baseseed 9=maxreps 10=verifyevery
*       11=permreps (production: 199; the smoke pass uses a smaller value to check plumbing)
* Column note (audit 2026-09-13): off_main/off_post0/off_post1 are e(perm_kernel) row 10
*   ("offcentre"), which HEAD sets to 1 only when the level() null band excludes 0 (the shift
*   interval is then blank) and leaves MISSING on a centred cell - so '.' in these columns means
*   "not off-centre", not "unread". analyze_mcri2.py scores blank ll/ul directly, never off_*.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local SHARD  = `2'
local OUTDIR `3'
local BD1 = `4'
local BD2 = `5'
local BD3 = `6'
local BD4 = `7'
local BASE = `8'
local MAXREPS = `9'
local VERIFY = `10'
local PR = `11'

global MCRI2_OUT `"`OUTDIR'/mc_ri2_s`SHARD'.csv"'
global MCRI2_SHARD = `SHARD'

*----------------------------------------------------------------- DGP (R1 + planted tau)
cap program drop mkdgp
program define mkdgp
  args seed tau
  clear
  set seed `seed'
  qui set obs 12
  gen i = _n
  gen u = runiform()
  sort u
  gen ttre = cond(_n <= 3, 5, cond(_n <= 6, 7, .))
  sort i
  gen a = rnormal()
  qui expand 10
  bysort i: gen t = _n
  gen y = a + .1*t + rnormal()
  if `tau' != 0 qui replace y = y + `tau'*(t >= ttre & ttre < .)
  drop a u
  sort i t
  ereturn clear
end

*----------------------------------------------------------------- row writer
* Reads every permutation surface off e() and appends one CSV row. Column positions are
* resolved BY NAME (colnumb) with a positional fallback recorded in namesok.
cap program drop emitrow
program define emitrow
  args phase arm tau rep dgpseed rseed rc secs
  local FLD ncomp nfail npos namesok ///
    p_main p_pre p_pre2 p_pre1 p_post0 p_post1 ///
    plo_main phi_main plo_pre phi_pre plo_pre2 phi_pre2 plo_pre1 phi_pre1 ///
    plo_post0 phi_post0 plo_post1 phi_post1 ///
    ll_main ul_main ll_post0 ul_post0 ll_post1 ul_post1 ///
    off_main off_post0 off_post1 sd_main sd_post0 sd_post1 ///
    b_main b_post0 b_post1 ase_main ase_post0 ase_post1 ///
    jointp tobs0 tobs1 praw0 praw1 pss0 pss1 psd0 psd1 ///
    cpeq cpabs cpraw cpleft ndraw
  foreach nm of local FLD {
    local `nm' .
  }
  if `rc' == 0 {
    tempname P
    cap matrix `P' = e(perm_p)
    if _rc == 0 {
      local i_main  = colnumb(`P', "main_att")
      local i_pre   = colnumb(`P', "pretreatment_att")
      local i_pre2  = colnumb(`P', "ES_att_pre2")
      local i_pre1  = colnumb(`P', "ES_att_pre1")
      local i_post0 = colnumb(`P', "ES_att_post0")
      local i_post1 = colnumb(`P', "ES_att_post1")
      local namesok 1
      local kk 0
      foreach z in main pre pre2 pre1 post0 post1 {
        local ++kk
        if `i_`z'' >= . | `i_`z'' == 0 {
          local i_`z' `kk'
          local namesok 0
        }
      }
      foreach z in main pre pre2 pre1 post0 post1 {
        local p_`z' = `P'[1, `i_`z'']
      }
      * one-sided shares
      tempname PL PH
      cap matrix `PL' = e(perm_p_lo)
      if _rc == 0 {
        foreach z in main pre pre2 pre1 post0 post1 {
          local plo_`z' = `PL'[1, `i_`z'']
        }
      }
      cap matrix `PH' = e(perm_p_hi)
      if _rc == 0 {
        foreach z in main pre pre2 pre1 post0 post1 {
          local phi_`z' = `PH'[1, `i_`z'']
        }
      }
      * shift intervals
      tempname CI
      cap matrix `CI' = e(perm_ci)
      if _rc == 0 {
        local rll = rownumb(`CI', "ll")
        local rul = rownumb(`CI', "ul")
        if `rll' >= . local rll 1
        if `rul' >= . local rul 2
        foreach z in main post0 post1 {
          local ll_`z' = `CI'[`rll', `i_`z'']
          local ul_`z' = `CI'[`rul', `i_`z'']
        }
      }
      * null-distribution SD and the off-centre (blank-interval) flag
      tempname SD KK
      cap matrix `SD' = e(perm_sd)
      if _rc == 0 {
        foreach z in main post0 post1 {
          local sd_`z' = `SD'[1, `i_`z'']
        }
      }
      cap matrix `KK' = e(perm_kernel)
      if _rc == 0 {
        local roff = rownumb(`KK', "offcentre")
        if `roff' >= . local roff 10
        foreach z in main post0 post1 {
          local off_`z' = `KK'[`roff', `i_`z'']
        }
      }
      * point estimates and the observed analytic SE (blank in fast mode)
      tempname BB AS
      cap matrix `BB' = e(b)
      if _rc == 0 {
        foreach z in main post0 post1 {
          local b_`z' = `BB'[1, `i_`z'']
        }
      }
      cap matrix `AS' = e(perm_ase_obs)
      if _rc == 0 {
        foreach z in main post0 post1 {
          local ase_`z' = `AS'[1, `i_`z'']
        }
      }
      * the four p-value conventions recomputed in-driver from the SAME draws
      tempname PD CB
      cap matrix `PD' = e(permdist_b)
      cap matrix `CB' = e(perm_comparator_b)
      cap confirm matrix `PD'
      local hasPD = (_rc == 0)
      cap confirm matrix `CB'
      if `hasPD' & _rc == 0 {
        local kcol = `i_main'
        mata: mcri2_conv("`PD'", "`CB'", `kcol')
      }
    }
    * fixed-metric joint pre-trend test
    cap local jointp = e(joint_horiz)
    * step-down record over the post horizons
    tempname SW
    cap matrix `SW' = e(perm_stepdown)
    if _rc == 0 {
      local r_t = rownumb(`SW', "t_obs")
      local r_r = rownumb(`SW', "p_raw")
      local r_s = rownumb(`SW', "p_singlestep")
      local r_d = rownumb(`SW', "p_stepdown")
      if `r_t' >= . local r_t 1
      if `r_r' >= . local r_r 2
      if `r_s' >= . local r_s 3
      if `r_d' >= . local r_d 4
      local c0 = colnumb(`SW', "ES_att_post0")
      local c1 = colnumb(`SW', "ES_att_post1")
      if `c0' >= . local c0 1
      if `c1' >= . local c1 2
      local tobs0 = `SW'[`r_t', `c0']
      local tobs1 = `SW'[`r_t', `c1']
      local praw0 = `SW'[`r_r', `c0']
      local praw1 = `SW'[`r_r', `c1']
      local pss0  = `SW'[`r_s', `c0']
      local pss1  = `SW'[`r_s', `c1']
      local psd0  = `SW'[`r_d', `c0']
      local psd1  = `SW'[`r_d', `c1']
    }
    cap local ncomp = e(perm_completed)
    cap local nfail = e(perm_failed)
    cap local npos  = e(perm_npossible)
  }
  local row "`phase',`arm',`tau',${MCRI2_SHARD},`rep',`dgpseed',`rseed',`rc',`secs'"
  foreach nm of local FLD {
    local x = strtrim(string(``nm'', "%18.12g"))
    if "`x'" == "" local x .
    local row "`row',`x'"
  }
  tempname fh
  file open `fh' using "${MCRI2_OUT}", write append
  file write `fh' `"`row'"' _n
  file close `fh'
end

*----------------------------------------------------------------- mata: the conventions
* Replicates wooldid_permtail's four p-value conventions on the SAME accumulator column:
*   equal-tailed (default)  min(1, 2 min(plo, phi)),  plo = (1+m0)/(n+1), phi = (1+m1)/(n+1)
*   permabsp                (1 + ma)/(n + 1),         ma = #{|s_r| >= |s_obs|}
*   permrawp                ma/n
*   permside(left)          plo
* Non-exact branch only (permreps(), never permexact), which is what this driver runs.
mata:
void mcri2_conv(string scalar pdname, string scalar cbname, real scalar k)
{
  real matrix A
  real colvector v
  real scalar n, sobs, m0, m1, ma, plo, phi
  A = st_matrix(pdname)
  sobs = st_matrix(cbname)[1, k]
  v = select(A[., k], A[., k] :< .)
  n = rows(v)
  if (n < 2 | sobs >= .) return
  m0 = sum(v :<= sobs)
  m1 = sum(v :>= sobs)
  ma = sum(abs(v) :>= abs(sobs))
  plo = (1 + m0)/(n + 1)
  phi = (1 + m1)/(n + 1)
  st_local("ndraw",  strofreal(n))
  st_local("cpeq",   strofreal(min((1, 2*min((plo, phi)))), "%21x"))
  st_local("cpabs",  strofreal((1 + ma)/(n + 1), "%21x"))
  st_local("cpraw",  strofreal(ma/n, "%21x"))
  st_local("cpleft", strofreal(plo, "%21x"))
}
end

*----------------------------------------------------------------- header (once per shard)
cap confirm file "${MCRI2_OUT}"
if _rc {
  tempname fh0
  file open `fh0' using "${MCRI2_OUT}", write replace
  file write `fh0' "phase,arm,tau,shard,rep,dgpseed,rseed,rc,secs," ///
    "ncomp,nfail,npos,namesok," ///
    "p_main,p_pre,p_pre2,p_pre1,p_post0,p_post1," ///
    "plo_main,phi_main,plo_pre,phi_pre,plo_pre2,phi_pre2,plo_pre1,phi_pre1," ///
    "plo_post0,phi_post0,plo_post1,phi_post1," ///
    "ll_main,ul_main,ll_post0,ul_post0,ll_post1,ul_post1," ///
    "off_main,off_post0,off_post1,sd_main,sd_post0,sd_post1," ///
    "b_main,b_post0,b_post1,ase_main,ase_post0,ase_post1," ///
    "jointp,tobs0,tobs1,praw0,praw1,pss0,pss1,psd0,psd1," ///
    "cpeq,cpabs,cpraw,cpleft,ndraw" _n
  file close `fh0'
}

local ESOPT esprelength(2) espostlength(2)

*================================================================== D1: sharp null
timer clear 91
local el = 0
forval rep = 1/`MAXREPS' {
  if `el' > `BD1' continue, break
  local dgpseed = `BASE' + `SHARD'*100000 + `rep'
  local rseed   = `BASE' + 500000 + `SHARD'*100000 + `rep'
  timer clear 91
  timer on 91
  mkdgp `dgpseed' 0
  cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
      permute(shuffle) permreps(`PR') rseed(`rseed') dots(0) permstepdown
  local rc = _rc
  timer off 91
  qui timer list 91
  local secs = r(t91)
  local el = `el' + `secs'
  emitrow D1 main 0 `rep' `dgpseed' `rseed' `rc' `secs'
  * verification arms: the package's own alternative conventions on the SAME draws
  if mod(`rep', `VERIFY') == 0 {
    foreach av in absp rawp left {
      local opt = cond("`av'" == "absp", "permabsp", cond("`av'" == "rawp", "permrawp", "permside(left)"))
      timer clear 91
      timer on 91
      mkdgp `dgpseed' 0
      cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
          permute(shuffle) permreps(`PR') rseed(`rseed') dots(0) `opt'
      local rc2 = _rc
      timer off 91
      qui timer list 91
      local secs2 = r(t91)
      local el = `el' + `secs2'
      emitrow D1 `av' 0 `rep' `dgpseed' `rseed' `rc2' `secs2'
    }
  }
  di "MCRI2_PROGRESS D1 shard=`SHARD' rep=`rep' rc=`rc' el=" %7.1f `el' " / `BD1'"
}
di "MCRI2_PHASE_DONE D1 shard=`SHARD' elapsed=" %8.1f `el'

*================================================================== D2: planted tau = 1
local el = 0
forval rep = 1/`MAXREPS' {
  if `el' > `BD2' continue, break
  local dgpseed = `BASE' + `SHARD'*100000 + `rep'
  local rseed   = `BASE' + 500000 + `SHARD'*100000 + `rep'
  timer clear 91
  timer on 91
  mkdgp `dgpseed' 1
  cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
      permute(shuffle) permreps(`PR') rseed(`rseed') dots(0)
  local rc = _rc
  timer off 91
  qui timer list 91
  local secs = r(t91)
  local el = `el' + `secs'
  emitrow D2 main 1 `rep' `dgpseed' `rseed' `rc' `secs'
  di "MCRI2_PROGRESS D2 shard=`SHARD' rep=`rep' rc=`rc' el=" %7.1f `el' " / `BD2'"
}
di "MCRI2_PHASE_DONE D2 shard=`SHARD' elapsed=" %8.1f `el'

*================================================================== D3: planted tau = 3
local el = 0
forval rep = 1/`MAXREPS' {
  if `el' > `BD3' continue, break
  local dgpseed = `BASE' + `SHARD'*100000 + `rep'
  local rseed   = `BASE' + 500000 + `SHARD'*100000 + `rep'
  timer clear 91
  timer on 91
  mkdgp `dgpseed' 3
  cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
      permute(shuffle) permreps(`PR') rseed(`rseed') dots(0)
  local rc = _rc
  timer off 91
  qui timer list 91
  local secs = r(t91)
  local el = `el' + `secs'
  emitrow D3 main 3 `rep' `dgpseed' `rseed' `rc' `secs'
  di "MCRI2_PROGRESS D3 shard=`SHARD' rep=`rep' rc=`rc' el=" %7.1f `el' " / `BD3'"
}
di "MCRI2_PHASE_DONE D3 shard=`SHARD' elapsed=" %8.1f `el'

*================================================================== D4: tau = 1, permtstat
local el = 0
forval rep = 1/`MAXREPS' {
  if `el' > `BD4' continue, break
  local dgpseed = `BASE' + `SHARD'*100000 + `rep'
  local rseed   = `BASE' + 500000 + `SHARD'*100000 + `rep'
  timer clear 91
  timer on 91
  mkdgp `dgpseed' 1
  cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
      permute(shuffle) permreps(`PR') permtstat rseed(`rseed') dots(0)
  local rc = _rc
  timer off 91
  qui timer list 91
  local secs = r(t91)
  local el = `el' + `secs'
  emitrow D4 tstat 1 `rep' `dgpseed' `rseed' `rc' `secs'
  di "MCRI2_PROGRESS D4 shard=`SHARD' rep=`rep' rc=`rc' el=" %7.1f `el' " / `BD4'"
}
di "MCRI2_PHASE_DONE D4 shard=`SHARD' elapsed=" %8.1f `el'

di "MCRI2_SHARD_DONE `SHARD'"
exit, clear
