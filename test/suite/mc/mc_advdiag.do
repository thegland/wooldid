* mc_advdiag.do — adversarial oversize diagnostics (post-A1 finding, 2026-08-24).
* ARCHIVED (inventory note 2026-09-13, critic C2): a diagnostic of the PRE-BS-12 (free-intercept,
* q=6) adversarial statistic on cluster(i) containment designs. No runner, no analyzer, not in
* any README. It runs on HEAD (the SI-A1 error is now the PS-64 caution, so the
* wooldid_solved_issues.md note that it "needs safety(off)" is stale) but its diagnosis is of the
* old statistic. Kept as a record, not a component.
* A1 measured ~20.5% rejection at nominal .05 (G=40, singleton cohorts, 4/4 window);
* joint_horiz ~11.3%. These arms separate the candidate mechanisms:
*   D1: G=200, singleton cohorts, 4/4  — finite-G? (big G, V structure unchanged)
*   D2: G=40, 8 units/cohort, 4/4     — singleton-cell V amplification? (INV-1/O16)
*   D3: G=40, singleton, 2/2 window   — test-dimension scaling?
*   D4: G=100, 8 units/cohort, 4/4    — near-asymptotic clean reference
* All flat-null, iid errors, cluster(i), esfixedbaseperiod, adversarial + jointtests.
* args: 1=adopath  2=config (D1..D4)  3=shard  4=nreps  5=baseseed  6=outdir
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

local G 40
local U 1
local PRE 4
local POST 4
if "`CFG'" == "D1" local G 200
if "`CFG'" == "D2" local U 8
if "`CFG'" == "D3" {
  local PRE 2
  local POST 2
}
if "`CFG'" == "D4" {
  local G 100
  local U 8
}

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
  if `U' > 1 qui expand `U'
  qui expand 12
  bysort i: gen t = mod(_n - 1, 12) + 1
  gen y = a + .05*t + rnormal()
  drop a
  sort i t
  ereturn clear
  cap qui wooldid y i t ttre, att cluster(i) esprelength(`PRE') espostlength(`POST') esfixedbaseperiod adversarial jointtests
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
dis "MC_ADVDIAG_SHARD_DONE `CFG' `SHARD'"
exit, clear
