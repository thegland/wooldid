clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioE.txt"
cap erase "vioE.txt"
include "ckc.doi"
include "ckc2.doi"

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"

* ============ 1. predict / postestimation after a SUCCESSFUL wooldid ============
dis "########## 1. POSTESTIMATION AFTER A SUCCESSFUL RUN ##########"
use D1.dta, clear
qui wooldid y id t ttre, att att_i att_it att_itime `ES'
dis "OK_CMD=`e(cmd)' K=" colsof(e(b))
cap noi predict double yh1
dis "PREDICT_RC=" _rc
cap drop yh1
cap noi predict double yh2, xb
dis "PREDICTXB_RC=" _rc
cap drop yh2
cap noi estat vce
dis "ESTATVCE_RC=" _rc
cap noi estat summarize
dis "ESTATSUM_RC=" _rc
cap noi margins
dis "MARGINS_RC=" _rc
cap noi lincom main_att - main_att_i
dis "LINCOM_RC=" _rc " est=" r(estimate) " se=" r(se)
cap noi test main_att = main_att_i
dis "TESTCROSS_RC=" _rc " F=" r(F) " p=" r(p)
cap noi estimates store SS1
dis "ESTSTORE_RC=" _rc
cap noi estimates table SS1
dis "ESTTABLE_RC=" _rc

* ============ 2. level() interval column naming ============
dis "########## 2. LEVEL NAMING ##########"
foreach L in 90 95 97.5 99.9 99.99 {
  cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' level(`L')
  dis "LVL`L'_RC=" _rc
  local cn : colnames e(wdidmainresults)
  dis "LVLCOLS|`L'|`cn'|e(level)=" e(level)
  ckc "L`=subinstr("`L'",".","p",.)'"
}
* set level globally then omit level()
set level 90
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES'
dis "SETLVL90_RC=" _rc
local cn : colnames e(wdidmainresults)
dis "SETLVLCOLS|`cn'|e(level)=" e(level)
ckc SETLVL90
set level 95

* ============ 3. the 32-character ledger: 7- vs 8-digit subgroup values ============
dis "########## 3. NAME LENGTH LEDGER ##########"
use D1.dta, clear
gen long sg7 = 1000000 + mod(id,2)
gen long sg8 = 10000000 + mod(id,2)
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sg7)
dis "SG7_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b))
ckc SG7
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sg8)
dis "SG8_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b))
* 8-digit with att only (recon says rc 0)
cap noi wooldid y id t ttre, att subgroups(sg8)
dis "SG8ATT_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b))
ckc SG8ATT
* adversarial-by-sg names overflow at len(v)>=8 (ado 4349-4350) - unguarded per recon
cap noi wooldid y id t ttre, att `ES' subgroups(sg8) jointtests adversarial cluster(clu)
dis "SG8ADV_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b))
cap noi wooldid y id t ttre, att `ES' subgroups(sg7) jointtests adversarial cluster(clu)
dis "SG7ADV_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b))
ckc SG7ADV

* ============ 4. subgroup-value screens (are the regex assumptions guarded?) ============
dis "########## 4. SUBGROUP VALUE SCREENS ##########"
use D4.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sgpfx) cluster(clu)
dis "SGPFX_RC=" _rc
ckc SGPFX
ckc2 SGPFX
cap noi wooldid y id t ttre, att subgroups(sgneg)
dis "SGNEG_RC=" _rc
cap noi wooldid y id t ttre, att subgroups(sgfrac)
dis "SGFRAC_RC=" _rc
* negative / fractional ONLY on never-treated units (outside the screen's domain)
use D4.dta, clear
gen double sgnt = mod(id,2)
replace sgnt = -1   if missing(ttre) & mod(id,4)==0
replace sgnt = 0.5  if missing(ttre) & mod(id,4)==2
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sgnt)
dis "SGNT_RC=" _rc " sglist=`e(wdid_sglist)'"
ckc SGNT
ckc2 SGNT

* ============ 5. the sentinel's reach: matrix consumers ============
dis "########## 5. SENTINEL CONSUMERS ##########"
use D3.dta, clear
cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(untreatedZero)
dis "CTUZ_RC=" _rc
ckc CTUZ
ckc2 CTUZ
tempname VV
cap matrix `VV' = e(V)
local nsent = 0
local nzero = 0
forval a = 1/`=colsof(`VV')' {
  if `VV'[`a',`a'] <= -1e99 local ++nsent
  if `VV'[`a',`a'] == 0 local ++nzero
}
dis "CTUZ_SENT=`nsent' ZERO=`nzero'"
cap mata: st_local("nm", strofreal(sum(st_matrix("e(V)") :>= .)))
dis "CTUZ_VMISSCELLS=`nm'"
cap mata: st_local("ch", strofreal(sum(cholesky(makesymmetric(st_matrix("e(V)"))) :>= .)))
dis "CTUZ_CHOLMISS=`ch' rc=" _rc
* the twin copies: e(V_es_ame) vs e(V)
cap confirm matrix e(V_es_ame)
dis "HAS_V_ES_AME=" (_rc==0)
cap confirm matrix e(b_es_ame)
dis "HAS_B_ES_AME=" (_rc==0)

* ============ 6. bootstrap twin matrices must never diverge ============
dis "########## 6. BOOTSTRAP TWINS ##########"
use D1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(40) rseed(777)
dis "BT_RC=" _rc
cap confirm matrix e(boot_pctl_mainresults)
local h1 = (_rc==0)
cap confirm matrix e(boot_percentile_mainresults)
local h2 = (_rc==0)
dis "TWINS pctl=`h1' percentile=`h2'"
if `h1' & `h2' {
  tempname T1 T2
  matrix `T1' = e(boot_pctl_mainresults)
  matrix `T2' = e(boot_percentile_mainresults)
  cap mata: st_local("md", strofreal(max(abs(st_matrix("e(boot_pctl_mainresults)") - st_matrix("e(boot_percentile_mainresults)")))))
  dis "TWIN_MAXDIFF=`md'"
}
* rank of the replication VCE at reps 40 vs K
cap mata: st_local("rk", strofreal(rank(st_matrix("e(V)"))))
dis "BOOT40_RANK=`rk' K=" colsof(e(b)) " completed=" e(boot_completed)

dis "EDGE_DONE"
