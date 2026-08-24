* mc_adv.do — adversarial line-test size-calibration MC driver
* (PREREG_NEWFEAT_2026-08-24 §M2). One shard of one config; one CSV row per rep.
* Null DGPs: the true ES profile lies on a line (flat zero, or a planted line through
* the reference period), so the composite null of the adversarial test is TRUE.
* args: 1=adopath  2=config (A1 A2 A3a A3b A5 A6)  3=shard  4=nreps  5=baseseed  6=outdir
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

* config geometry: planted-line configs use t=1..9 (relyr -5..3) so the 4/4 ES window
* spans the support exactly (no endcap pooling — verified by the null-fidelity pretest)
local G 40
local T 12
local SLOPE 0
local ALT
if "`CFG'" == "A2" {
  local T 9
  local SLOPE .15
}
if "`CFG'" == "A3a" local ALT altmethod_adversarial
if "`CFG'" == "A3b" {
  local T 9
  local SLOPE .15
  local ALT altmethod_adversarial
}
if "`CFG'" == "A5" local G 10
local AR1 = ("`CFG'" == "A6")

tempname fh
file open `fh' using `"`OUTDIR'/`CFG'_s`SHARD'.csv"', write replace
file write `fh' "config,shard,rep,seed,rc,padv,pjoint" _n

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  set seed `seed'
  clear
  qui set obs `G'
  gen i = _n
  gen ttre = cond(mod(_n, 2) == 0, 6, .)
  gen a = rnormal()
  qui expand `T'
  bysort i: gen t = _n
  if `AR1' {
    bysort i (t): gen e = rnormal()
    bysort i (t): replace e = .7*e[_n-1] + sqrt(1 - .49)*rnormal() if _n > 1
  }
  else gen e = rnormal()
  * the planted line runs through the reference period (relyr -1): value slope*(relyr+1)
  gen y = a + .05*t + `SLOPE'*(t - 6 + 1)*(ttre < .) + e
  drop a e
  sort i t
  ereturn clear
  cap qui wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod adversarial `ALT' jointtests
  local rc = _rc
  local padv .
  local pjoint .
  if `rc' == 0 {
    local padv = e(adversarial_pval_att)
    local pjoint = e(joint_horiz)
  }
  local v1 : di %18.0g `padv'
  local v2 : di %18.0g `pjoint'
  file write `fh' "`CFG',`SHARD',`rep',`seed',`rc',`=trim("`v1'")',`=trim("`v2'")'" _n
}
file close `fh'
dis "MC_ADV_SHARD_DONE `CFG' `SHARD'"
exit, clear
