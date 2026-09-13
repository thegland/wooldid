* mc_adv2.do — A4: continuous-treatment adversarial size under a SOUND VCE
* (PREREG_NEWFEAT_2026-08-24 §M2 AMENDMENT A4). Cohort-level cluster(i) is excluded by
* the INV-1 finding (rank-1 V_es_ame: CRVE scores of cell dose slopes identically
* zero); robust on iid errors is a probe-verified full-rank configuration.
* WARNING (Stage A audit 2026-09-13, critic C1): the ame arm at 6 obs per cohort-by-period cell
* (`qui expand 6`) sits BELOW the ADVBOUND cell-size gate for the continuous adversarial test
* (REPORT_ADVBOUND_2026-08-26.md: 31.55% at 6 obs/cell -> 9.92% at 20 -> 5.95% PASS at 50, AB22;
* ruling GATE ON CELL SIZE), so the A4 ame band FAILs on HEAD by ruling (08-24: 159/504 =
* 31.55%). The att side (robust, 6 units/cell) is AB2's PASS regime and remains a valid size
* test. NOT rewritten pending an author ruling; see PREREG_NEWFEAT_2026-08-24.md NOTE (M2, M3b).
* args: 1=adopath  2=config (A4)  3=shard  4=nreps  5=baseseed  6=outdir
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

tempname fh
file open `fh' using `"`OUTDIR'/`CFG'_s`SHARD'.csv"', write replace
file write `fh' "config,shard,rep,seed,rc,padv,pjoint,pame,pjointC" _n

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  set seed `seed'
  clear
  qui set obs 30
  gen i = _n
  gen ttre = cond(mod(_n, 2) == 0, 6, .)
  gen a = rnormal()
  qui expand 6
  bysort i: gen unit = _n
  gen dose = exp(0.8*runiform() + 0.2*mod(i, 5))
  qui expand 12
  bysort i unit: gen t = _n
  gen y = a + .05*t + rnormal()
  drop a
  sort i unit t
  ereturn clear
  cap qui wooldid y i t ttre, att robust contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial jointtests
  local rc = _rc
  foreach nm in padv pjoint pame pjointC {
    local `nm' .
  }
  if `rc' == 0 {
    local padv = e(adversarial_pval_att)
    local pjoint = e(joint_horiz)
    local pame = e(adversarial_pval_ame)
    local pjointC = e(joint_Choriz)
  }
  local row "`CFG',`SHARD',`rep',`seed',`rc'"
  foreach nm in padv pjoint pame pjointC {
    local v : di %18.0g ``nm''
    local row "`row',`=trim("`v'")'"
  }
  file write `fh' "`row'" _n
}
file close `fh'
dis "MC_ADV2_SHARD_DONE `CFG' `SHARD'"
exit, clear
