* mc_advb.do — adversarial line-test SIZE BOUNDARY driver (AB campaign;
* PREREG_ADVBOUND_2026-08-25). One shard of one config; one CSV row per rep.
*
* The question: the adversarial test is known to be severely oversized under
* cohort-level clustering with singleton cells (A1-A6: 11-48% at nominal 5%), and
* correctly sized in exactly one measured configuration (A4 att: 5.75% under robust
* with 6-unit cells). This driver maps the boundary between those two regimes by
* taking the panel geometry, the VCE and the error structure as ARGUMENTS.
*
* The null is TRUE by construction in every arm: y carries no treatment effect at
* all, so the true event-study profile is flat zero, which lies on a line.
*
* args: 1=adopath 2=config 3=shard 4=nreps 5=baseseed 6=outdir
*       7=G (cohorts) 8=U (units per cohort) 9=T (periods)
*       10=vce token, one of: robust | cluster(i) | cluster(unitid)
*       11=errtype: iid | ar1 | cohortshock | unitre
*       12=mode: binary | contreat
*       13,14=extra wooldid options (e.g. altmethod_adversarial, safety(off))
*
* NOTE the unit identifier: `unitid' is GLOBALLY unique ((i-1)*U + unit). The
* within-cohort counter 1..U repeats across cohorts, so clustering on it would
* silently cluster on U meaningless cross-cohort groups - and nothing would error.
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
local G = `7'
local U = `8'
local T = `9'
local VCE `10'
local ERR `11'
local MODE `12'
local EXTRA `13' `14'

local ONSET = ceil(`T'/2)
local CONTOPT
if "`MODE'" == "contreat" local CONTOPT contreat(dose) contreatcontrols(basic)

tempname fh
file open `fh' using `"`OUTDIR'/`CFG'_s`SHARD'.csv"', write replace
file write `fh' "config,shard,rep,seed,rc,padv,pjoint,pame,pjointC" _n

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  set seed `seed'
  clear
  qui set obs `G'
  gen i = _n
  * half the cohorts treated at a common onset, half never-treated (comparison group)
  gen ttre = cond(mod(_n, 2) == 0, `ONSET', .)
  * cohort-level intercept: absorbed by the cohort FE, kept for comparability with
  * the A-series (inert by construction)
  gen a = rnormal()
  qui expand `U'
  bysort i: gen unit = _n
  gen unitid = (i - 1)*`U' + unit
  * unit-level random intercept is NOT absorbed (no unit FE), so it persists in the
  * residuals and invalidates robust: it belongs only to the unitre violation arm
  gen ure = cond("`ERR'" == "unitre", rnormal(), 0)
  * dose varies WITHIN each cohort-period cell (across units), so the ame is
  * identified and the SI-A4 void screen cannot fire
  gen dose = exp(0.8*runiform() + 0.2*mod(i, 5))
  qui expand `T'
  bysort i unitid: gen t = _n

  if "`ERR'" == "ar1" {
    bysort unitid (t): gen e = rnormal()
    bysort unitid (t): replace e = .7*e[_n-1] + sqrt(1 - .49)*rnormal() if _n > 1
  }
  else gen e = rnormal()
  * cohort-by-period common shock: correlated ACROSS units within a cohort-period,
  * which is exactly what unit-level clustering and robust both fail to capture
  if "`ERR'" == "cohortshock" {
    tempvar cs
    bysort i t: gen `cs' = rnormal() if _n == 1
    bysort i t (`cs'): replace `cs' = `cs'[1]
    qui replace e = e + `cs'
  }

  * NO treatment effect and no dose effect: the null of the adversarial test (all
  * relative-time coefficients on a single line) is exactly true
  gen y = a + ure + .05*t + e
  drop a e ure
  sort i unitid t
  ereturn clear

  cap qui wooldid y i t ttre, att `VCE' `CONTOPT' esprelength(4) espostlength(4) ///
      esfixedbaseperiod adversarial jointtests `EXTRA'
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
dis "MC_ADVB_SHARD_DONE `CFG' `SHARD'"
exit, clear
