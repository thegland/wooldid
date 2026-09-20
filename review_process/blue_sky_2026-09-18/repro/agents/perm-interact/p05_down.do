clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear

* ============ A. gteffects x permute ============
dis "==== A1: plain gteffects, no permute (reference) ===="
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) gteffects
dis "A1_RC=" _rc
cap confirm matrix e(gteffects)
dis "A1_HASGT=" (_rc==0)
if _rc == 0 {
  matrix GT0 = e(gteffects)
  dis "A1_GTDIM=" rowsof(GT0) "x" colsof(GT0)
  matrix list GT0, format(%9.4f) nohalf
}

dis "==== A2: gteffects UNDER permute ===="
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) gteffects permute(shuffle) permreps(30) rseed(777)
dis "A2_RC=" _rc
dis "A2_GTBLANK=" "`e(perm_gt_blank)'"
cap confirm matrix e(gteffects)
dis "A2_HASGT=" (_rc==0)
if _rc == 0 {
  matrix GT1 = e(gteffects)
  dis "A2_GTDIM=" rowsof(GT1) "x" colsof(GT1)
  matrix list GT1, format(%9.4f) nohalf
}

* ============ B. e(V) surface under the three permute modes ============
cap prog drop vchk
prog def vchk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) `opts'
  dis "VCHK|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "VCHK|`tag'|vcvstate=" "`e(wdid_vcvstate)'" "|vce=" "`e(vce)'" "|vcetype=" "`e(vcetype)'"
  cap confirm matrix e(V)
  dis "VCHK|`tag'|HASV=" (_rc==0)
  cap confirm matrix e(V_es_att)
  dis "VCHK|`tag'|HASVESATT=" (_rc==0)
  cap confirm matrix e(V_es_analytic_att)
  dis "VCHK|`tag'|HASVESANALYTIC=" (_rc==0)
  cap noi test main_att = 0
  dis "VCHK|`tag'|TEST_RC=" _rc "|p=" r(p)
  cap noi lincom main_att
  dis "VCHK|`tag'|LINCOM_RC=" _rc "|se=" r(se)
  cap noi estat vce
  dis "VCHK|`tag'|ESTATVCE_RC=" _rc
end
vchk fastmode
vchk ase        , opts(permanalyticse)
vchk ts         , opts(permtstat)

* ============ C. honestdid downstream after permanalyticse ============
dis "==== C: honestdid on a permutation run (permanalyticse) ===="
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) permanalyticse
dis "C_RC=" _rc
cap confirm matrix e(V_es_att)
if _rc == 0 {
  matrix VESP = e(V_es_att)
  dis "C_VES_DIM=" rowsof(VESP) "x" colsof(VESP)
  matrix list VESP, format(%12.8f)
}
* the analytic reference run
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2)
dis "C2_RC=" _rc
cap confirm matrix e(V_es_att)
if _rc == 0 {
  matrix VESA = e(V_es_att)
  dis "C2_VES_DIM=" rowsof(VESA) "x" colsof(VESA)
  matrix list VESA, format(%12.8f)
}
dis "P05_DONE"
