* ===== the STRUCTURAL rank of the bootstrap e(V), and e(df_r) vs the printed table =====
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioK.txt"
cap erase "vioK.txt"
include "ckc.doi"

cap prog drop RANKLINE
prog def RANKLINE
  args tag
  local rk = .
  local rka = .
  cap mata: st_local("rk", strofreal(rank(st_matrix("e(V)"))))
  cap mata: st_local("rka", strofreal(rank(st_matrix("e(V_analytic)"))))
  local bk = .
  cap confirm matrix e(b)
  if _rc == 0 local bk = colsof(e(b))
  dis "RANK|`tag'|K=`bk'|rank(e(V))=`rk'|rank(e(V_analytic))=`rka'|reps=" e(boot_reps) "|completed=" e(boot_completed) "|boot_N_V=" e(boot_N_V) "|v_posted=" e(boot_v_posted) "|df_r=" e(df_r) "|boot_df=" e(boot_df)
end

use D1.dta, clear
local ES "esfixedbaseperiod esprelength(3) espostlength(3)"

* ---- the rank of the ANALYTIC e(V) and of the bootstrap e(V), by estimand set ----
* the point: the coefficient vector is structurally lower-dimensional than K, so the
* replication covariance is rank-deficient AT EVERY reps(), while the analytic
* diag(se^2) always claims full rank K.
cap noi wooldid y id t ttre, att `ES' cluster(clu)
dis "S1_RC=" _rc
RANKLINE ANALYTIC_att_only
cap noi wooldid y id t ttre, att att_i `ES' cluster(clu)
dis "S2_RC=" _rc
RANKLINE ANALYTIC_att_atti
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu)
dis "S3_RC=" _rc
RANKLINE ANALYTIC_all4

cap noi wooldid y id t ttre, att `ES' cluster(clu) bootstrap(standard) reps(60) rseed(11)
dis "B1_RC=" _rc
RANKLINE BOOT60_att_only
cap noi ckc BOOT60_att_only
cap noi wooldid y id t ttre, att att_i `ES' cluster(clu) bootstrap(standard) reps(60) rseed(11)
dis "B2_RC=" _rc
RANKLINE BOOT60_att_atti
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(60) rseed(11)
dis "B3_RC=" _rc
RANKLINE BOOT60_all4
cap noi ckc BOOT60_all4
* what a user's joint test does on it
cap noi test ES_att_post0 ES_att_post1 ES_att_post2
dis "B3_TEST3 p=" r(p) " df=" r(df) " F=" r(F)
cap noi test main_att main_att_it main_att_i main_att_itime
dis "B3_TESTEST p=" r(p) " df=" r(df) " F=" r(F)
cap noi lincom main_att - main_att_it
dis "B3_LINCOM est=" r(estimate) " se=" r(se) " rc=" _rc
cap noi estat vce
dis "B3_ESTATVCE_RC=" _rc
cap mata: st_local("iv", strofreal(sum(invsym(st_matrix("e(V)")) :>= .)))
dis "B3_INVSYM_MISSINGCELLS=`iv'"
* subgroups: K=60
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) subgroups(sg) bootstrap(standard) reps(40) rseed(11)
dis "B4_RC=" _rc
RANKLINE BOOTSG40
cap noi ckc BOOTSG40
cap noi test main_SG0_att main_SG1_att main_SG0_att_i main_SG1_att_i
dis "B4_TEST4 p=" r(p) " df=" r(df) " F=" r(F)

* ---- e(df_r) vs the printed table: lincom must reproduce the table row ----
cap prog drop CMP
prog def CMP
  args tag row
  tempname MM
  cap confirm matrix e(wdidmainresults)
  if _rc exit
  matrix `MM' = e(wdidmainresults)
  local rn : rowfullnames `MM'
  local i = 0
  local hit = 0
  foreach r of local rn {
    local ++i
    if `"`r'"' == "main:`row'" local hit = `i'
  }
  if `hit' == 0 exit
  local nm = subinstr("main:`row'",":","_",.)
  cap noi lincom `nm'
  dis "CMP|`tag'|`row'|table est=" `MM'[`hit',1] " se=" `MM'[`hit',2] " lb=" `MM'[`hit',5] " ub=" `MM'[`hit',6] "|lincom est=" r(estimate) " se=" r(se) " lb=" r(lb) " ub=" r(ub) "|df_r=" e(df_r)
end

use D1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu)
dis "C1_RC=" _rc " df_r=" e(df_r)
CMP clustered att
CMP clustered att_itime
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) impvalues
dis "C2_RC=" _rc " df_r=" e(df_r) " imdf_min=" e(imdf_min) " imdf_max=" e(imdf_max) " df_r_analytic=" e(df_r_analytic)
CMP impvalues att
CMP impvalues att_itime
cap noi wooldid ypos id t ttre, att att_i att_it att_itime `ES' poisson
dis "C3_RC=" _rc " df_r=" e(df_r)
CMP poisson att
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(40) rseed(11)
dis "C4_RC=" _rc " df_r=" e(df_r) " boot_df=" e(boot_df)
CMP boot att

dis "RANK_DONE"
