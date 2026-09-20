* ===== level() handling, the interval-column naming contract, and weights =====
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioL.txt"
cap erase "vioL.txt"
include "ckc.doi"
include "ckc2.doi"

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"
use D1.dta, clear

* ---- 1. the level() grid: wooldid vs an ordinary Stata estimation command ----
dis "########## LEVEL GRID ##########"
foreach L in 90 95 99 97.5 99.5 68.5 99.9 90.1 99.99 99.95 68.3 92.5 {
  cap noi regress y t
  cap noi regress y t, level(`L')
  local rreg = _rc
  cap noi set level `L'
  local rset = _rc
  set level 95
  cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' level(`L')
  local rwd = _rc
  local cn ""
  if `rwd' == 0 local cn : colnames e(wdidmainresults)
  dis "LEVELGRID|`L'|regress=`rreg'|setlevel=`rset'|wooldid=`rwd'|elevel=" e(level) "|cols=`cn'"
  if `rwd' == 0 ckc LVL`=subinstr("`L'",".","p",.)'
}
set level 95

* the same via `set level` then no level() option (the other route into wdlvltxt)
foreach L in 99.9 90.1 {
  cap noi set level `L'
  local rset = _rc
  if `rset' == 0 {
    cap noi wooldid y id t ttre, att att_i att_it att_itime `ES'
    local cn : colnames e(wdidmainresults)
    dis "SETLEVEL|`L'|wooldid=" _rc "|elevel=" e(level) "|cols=`cn'"
  }
  else dis "SETLEVEL|`L'|set level itself rc=`rset'"
  set level 95
}

* ---- 2. weights (correct Stata placement) ----
dis "########## WEIGHTS ##########"
cap noi wooldid y id t ttre [aw=awt], att att_i att_it att_itime `ES'
dis "AW_RC=" _rc " wtype=`e(wtype)' wexp=`e(wexp)' vce=`e(vce)'"
if _rc == 0 ckc AWT
cap noi wooldid y id t ttre [pw=pwt], att att_i att_it att_itime `ES'
dis "PW_RC=" _rc " wtype=`e(wtype)' wexp=`e(wexp)' vce=`e(vce)' vcetype=`e(vcetype)'"
if _rc == 0 ckc PWT
gen int fwx = 1 + mod(id,3)
cap noi wooldid y id t ttre [fw=fwx], att `ES'
dis "FW_RC=" _rc
cap noi wooldid y id t ttre [iw=awt], att `ES'
dis "IW_RC=" _rc

* ---- 3. the subgroup coefficient-name contract the help states ----
dis "########## SUBGROUP COEFFICIENT NAMES ##########"
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sg) cluster(clu)
dis "SGRUN_RC=" _rc
dis "BCOLNAMES|`: colfullnames e(b)'"
* what sthlp:1306 tells the user to type:
cap noi lincom main_att_s1
dis "LINCOM_HELPFORM_RC=" _rc
cap noi lincom main_att_s0
dis "LINCOM_HELPFORM0_RC=" _rc
* what the code actually produces:
cap noi lincom main_SG1_att
dis "LINCOM_CODEFORM_RC=" _rc " est=" r(estimate)
cap noi test main_SG1_att = main_SG0_att
dis "TEST_CODEFORM_RC=" _rc
* the OTHER subgroup-keyed names DO use the documented suffix:
foreach n in b_es_att_s1 numpre_att_s1 wdidsgresults_1 escomp_s1 {
  cap confirm matrix e(`n')
  local m = (_rc==0)
  cap confirm scalar e(`n')
  local s = (_rc==0)
  dis "SGNAME|`n'|matrix=`m'|scalar=`s'"
}
ckc SGNAMES
ckc2 SGNAMES

dis "LEVEL_DONE"
