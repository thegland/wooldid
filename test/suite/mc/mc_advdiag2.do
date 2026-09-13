* mc_advdiag2.do — direct V-understatement verification for the adversarial oversize
* finding.
* ARCHIVED (inventory note 2026-09-13, critic C2): 2026-08-24 diagnostic of the PRE-BS-12
* statistic on containment designs; no runner, no analyzer, not in any README; runs on HEAD
* (SI-A1 error is now the PS-64 caution) but diagnoses the old statistic. Record, not component. Records per-rep e(b_es_att) entries, diag(e(V_es_att)), and the eigenvalue
* range of V_es_att, so the analyzer can compare MC-truth Var(b_h) and contrast
* variances against the reported V, per design.
*   V1: A1 shape (G=40, singleton cohorts, 4/4)
*   V2: D2 shape (G=40, 8 units/cohort, 4/4)
*   D5: G=40, 8 units/cohort, 2/2 window — rates only (padv/pjoint), tests the
*       prediction that multi-unit contrast-direction damage persists at low q
* args: 1=adopath  2=config (V1 V2 D5)  3=shard  4=nreps  5=baseseed  6=outdir
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
if inlist("`CFG'", "V2", "D5") local U 8
if "`CFG'" == "D5" {
  local PRE 2
  local POST 2
}
local K = `PRE' + `POST'

tempname fh
file open `fh' using `"`OUTDIR'/`CFG'_s`SHARD'.csv"', write replace
local hdr "config,shard,rep,seed,rc,padv,pjoint,evmin,evmax"
forval j = 1/`K' {
  local hdr "`hdr',b`j',v`j'"
}
file write `fh' "`hdr'" _n

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
  local evmin .
  local evmax .
  forval j = 1/`K' {
    local b`j' .
    local v`j' .
  }
  if `rc' == 0 {
    local padv = e(adversarial_pval_att)
    local pjoint = e(joint_horiz)
    tempname B V
    cap matrix `B' = e(b_es_att)
    cap matrix `V' = e(V_es_att)
    if _rc == 0 & colsof(`B') == `K' {
      forval j = 1/`K' {
        local b`j' = `B'[1, `j']
        local v`j' = `V'[`j', `j']
      }
      mata: st_local("evmin", strofreal(min(symeigenvalues(st_matrix(st_local("V")))), "%18.0g"))
      mata: st_local("evmax", strofreal(max(symeigenvalues(st_matrix(st_local("V")))), "%18.0g"))
    }
    cap matrix drop `B' `V'
  }
  local row "`CFG',`SHARD',`rep',`seed',`rc'"
  foreach nm in padv pjoint evmin evmax {
    local v : di %18.0g ``nm''
    local row "`row',`=trim("`v'")'"
  }
  forval j = 1/`K' {
    local vb : di %18.0g `b`j''
    local vv : di %18.0g `v`j''
    local row "`row',`=trim("`vb'")',`=trim("`vv'")'"
  }
  file write `fh' "`row'" _n
}
file close `fh'
dis "MC_ADVDIAG2_SHARD_DONE `CFG' `SHARD'"
exit, clear
