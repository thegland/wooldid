*=======================================================================
* repro_joint_size.do -- SELF-CONTAINED.  Size of the JOINT PRE-TREND TEST
* (and of the adversarial line test) under a TRUE SHARP NULL at G=10, under
* the three ways wooldid can form them:
*    analytic                     (documented oversized: NF-1b 65.7% at G=10)
*    permute(shuffle)             (rewired to a permutation rank test 2026-09)
*    bootstrap(standard)          (rewired to a fixed-metric rank test 2026-09)
* Under this DGP the PERMUTATION version is EXACT by theorem (the onset
* assignment is drawn uniformly and y does not depend on it), so a rejection
* rate materially above 5% there is a defect, not an approximation error.
*   stata-mp -e do repro_joint_size.do          (~35 min at NDRAW 40)
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAW 2

cap prog drop mkdgp
prog def mkdgp
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `ngrp'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  local nc = round(`ngrp'*0.2)
  gen double ttre = .
  replace ttre = 2004 if rk <= `nc'
  replace ttre = 2005 if rk >    `nc' & rk <= 2*`nc'
  replace ttre = 2006 if rk >  2*`nc' & rk <= 3*`nc'
  gen double ai = rnormal()
  drop u rk
  expand `nper'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75) if _n>1
  gen double y = ai + ft + e         /* no treatment term anywhere */
  keep st t ttre y
end

tempname P
postfile `P' int draw str10 engine double pjoint double padv int rc using rt3_out.dta, replace
forvalues d = 1/`NDRAW' {
  local sx = 2000000 + 10*100000 + 10*1000 + `d'
  local rs = 700000 + `d'*7 + 10
  mkdgp , ngrp(10) nper(8) seedx(`sx')
  tempfile dat
  qui save "`dat'", replace
  foreach e in analytic permute bootstrap {
    qui use "`dat'", clear
    ereturn clear
    local o
    if "`e'" == "permute"   local o permute(shuffle) permreps(79) rseed(`rs')
    if "`e'" == "bootstrap" local o bootstrap(standard) reps(79) rseed(`rs')
    cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) ///
        esfixedbaseperiod jointtests adversarial `o'
    local rc = _rc
    local pj = .
    local pa = .
    if `rc' == 0 {
      local pj = e(joint_horiz)
      local pa = e(adversarial_pval_att)
    }
    post `P' (`d') ("`e'") (`pj') (`pa') (`rc')
  }
  if mod(`d',10)==0 dis "... draw `d' done"
}
postclose `P'
use rt3_out.dta, clear
gen byte rj = (pjoint <= .05 + 1e-9) if pjoint < .
gen byte ra = (padv   <= .05 + 1e-9) if padv   < .
dis "=== rejection at nominal 5% of a TRUE sharp null ==="
table engine, stat(mean rj) stat(count pjoint) stat(mean ra) stat(count padv)
