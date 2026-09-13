* mc_ri.do — RI size-calibration MC driver (PREREG_NEWFEAT_2026-08-24 §M1).
* One shard: loops reps of one config, one wooldid permute() run per rep, appends one
* CSV row per rep. The TRUE assignment is drawn from the SAME mechanism the config
* hands to permute(), and the outcome never depends on the assignment (sharp null).
* args: 1=path to wooldid.ado (file or dir)  2=config (R1..R7)  3=shard (int)
*       4=nreps  5=baseseed  6=outdir
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local CFG `2'
local SHARD = `3'
local NREPS = `4'
local BASE = `5'
local OUTDIR `6'

* ---- per-config data-generating programs; the seed is set by the rep loop ----
cap program drop mkdgp
program define mkdgp
  args cfg
  clear
  if inlist("`cfg'", "R1", "R2", "R3", "R7") {
    qui set obs 12
    gen i = _n
    * truth = uniform shuffle of the label multiset {5 x3, 7 x3, . x6}
    gen u = runiform()
    sort u
    gen ttre = cond(_n <= 3, 5, cond(_n <= 6, 7, .))
    sort i
    gen a = rnormal()
    if "`cfg'" == "R3" gen s = 1 + 3*(i-1)/11
    if "`cfg'" == "R7" gen w = 1 + mod(i, 4)
    qui expand 10
    bysort i: gen t = _n
    if "`cfg'" == "R2" {
      bysort i (t): gen e = rnormal()
      bysort i (t): replace e = .7*e[_n-1] + sqrt(1 - .49)*rnormal() if _n > 1
      gen y = a + .1*t + e
      drop e
    }
    else if "`cfg'" == "R3" gen y = a + .1*t + s*rnormal()
    else gen y = a + .1*t + rnormal()
    drop a u
    cap drop s
  }
  else if "`cfg'" == "R4" {
    qui set obs 16
    gen i = _n
    gen str = 1 + (_n > 8)
    gen u = runiform()
    sort str u
    by str: gen rk = _n
    gen ttre = cond(str == 1 & rk <= 3, 5, cond(str == 2 & rk <= 3, 7, .))
    sort i
    gen a = rnormal()
    qui expand 10
    bysort i: gen t = _n
    gen y = a + .1*t + rnormal()
    drop a u rk
  }
  else if "`cfg'" == "R5" {
    qui set obs 10
    gen i = _n
    gen u = runiform()
    gen ttre = cond(u < .3, 5, cond(u < .6, 7, .))
    gen a = rnormal()
    qui expand 10
    bysort i: gen t = _n
    gen y = a + .1*t + rnormal()
    drop a u
  }
  else if "`cfg'" == "R6" {
    qui set obs 6
    gen i = _n
    gen u = runiform()
    sort u
    gen ttre = cond(_n <= 3, 4, .)
    sort i
    gen a = rnormal()
    qui expand 8
    bysort i: gen t = _n
    gen y = a + .1*t + rnormal()
    drop a u
  }
  sort i t
  ereturn clear
end

tempname fh
file open `fh' using `"`OUTDIR'/`CFG'_s`SHARD'.csv"', write replace
file write `fh' "config,shard,rep,seed,attempts,rc,p,nfail,ncomp,npos,lb,ub,comp,tlb,tub" _n

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  local attempts 0
  local rc -1
  * R5's mechanism can draw an inestimable observed assignment (e.g. no controls);
  * conditioning on estimability (redraw, seed advanced by 1000000) matches the
  * engine's own exclusion of inestimable permutation draws — see the prereg.
  while `rc' != 0 & `attempts' < 5 {
    local ++attempts
    set seed `=`seed' + (`attempts'-1)*1000000'
    mkdgp `CFG'
    if "`CFG'" == "R3" {
      cap qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(99) permtstat permtci
    }
    else if "`CFG'" == "R4" {
      cap qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(99) permstrata(str)
    }
    else if "`CFG'" == "R5" {
      cap qui wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(99) permprobs((ttre 5 7 .) (1 .3 .3 .4) (2 .3 .3 .4) (3 .3 .3 .4) (4 .3 .3 .4) (5 .3 .3 .4) (6 .3 .3 .4) (7 .3 .3 .4) (8 .3 .3 .4) (9 .3 .3 .4) (10 .3 .3 .4))
    }
    else if "`CFG'" == "R6" {
      cap qui wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
    }
    else if "`CFG'" == "R7" {
      cap qui wooldid y i t ttre [aw=w], att cluster(i) permute(shuffle) permreps(99)
    }
    else {
      cap qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(99)
    }
    local rc = _rc
    * only R5 redraws; every other config records the first attempt whatever its rc
    if "`CFG'" != "R5" continue, break
  }
  foreach nm in p nfail ncomp npos lb ub comp tlb tub {
    local `nm' .
  }
  if `rc' == 0 {
    local p = e(perm_p_att)
    local nfail = e(perm_failed)
    local ncomp = e(perm_completed)
    local npos = e(perm_npossible)
    * CS-57: e(perm_bandlo/bandhi_att) are the canonical names of the level() null band;
    * e(perm_lb/ub_att) are their deprecated twins (identical values, kernel rows 11/12)
    local lb = e(perm_bandlo_att)
    local ub = e(perm_bandhi_att)
    cap local tlb = e(perm_tlb_att)
    cap local tub = e(perm_tub_att)
    tempname CMP
    cap matrix `CMP' = e(perm_comparator)
    if _rc == 0 local comp = `CMP'[1, 1]
    cap matrix drop `CMP'
  }
  local row "`CFG',`SHARD',`rep',`seed',`attempts',`rc'"
  foreach nm in p nfail ncomp npos lb ub comp tlb tub {
    local v : di %18.0g ``nm''
    local row "`row',`=trim("`v'")'"
  }
  file write `fh' "`row'" _n
}
file close `fh'
dis "MC_RI_SHARD_DONE `CFG' `SHARD'"
exit, clear
