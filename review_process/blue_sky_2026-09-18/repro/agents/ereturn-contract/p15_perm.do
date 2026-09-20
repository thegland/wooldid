* ===== permute(): does e(perm_p_<est>) agree with the table and with e(perm_p)? =====
* recon sec 4b item 3 asked whether a SUBSET estimand set can mis-map the column.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioP.txt"
cap erase "vioP.txt"
include "ckc.doi"

clear
set obs 40
gen long id = _n
gen byte grp = mod(id,4)
gen double ttre = .
replace ttre = 3 if grp==1
replace ttre = 5 if grp==2
replace ttre = 7 if grp==3
bysort grp (id): gen long clu = grp*1000 + ceil(_n/3)
sort id
gen double ufe = rnormal()
expand 9
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()

cap prog drop PERMX
prog def PERMX
  args tag
  if "`e(perm_method)'" == "" {
    dis "PERMX|`tag'|NO PERMUTE"
    exit
  }
  tempname M P
  cap confirm matrix e(wdidmainresults)
  if _rc exit
  matrix `M' = e(wdidmainresults)
  local rn : rowfullnames `M'
  local cn : colnames `M'
  local pc : list posof "p" in cn
  local hasP = 0
  cap confirm matrix e(perm_p)
  if _rc == 0 {
    local hasP = 1
    matrix `P' = e(perm_p)
  }
  local pn ""
  if `hasP' local pn : colnames `P'
  dis "PERMSET|`tag'|estimands=`e(wdid_estimands)'|perm_p cols=`pn'|completed=" e(perm_completed)
  foreach est in `e(wdid_estimands)' {
    local sc = .
    cap local sc = e(perm_p_`est')
    local trow = 0
    local i = 0
    foreach r of local rn {
      local ++i
      if `"`r'"' == "main:`est'" local trow = `i'
    }
    local tval = .
    if `trow' > 0 & `pc' > 0 local tval = `M'[`trow',`pc']
    local vval = .
    if `hasP' {
      local q : list posof "main_`est'" in pn
      if `q' > 0 local vval = `P'[1,`q']
    }
    local flag "OK"
    if !missing(`sc') & !missing(`tval') {
      if abs(`sc'-`tval') > 1e-12 local flag "SCALAR_MISMATCH"
    }
    if !missing(`vval') & !missing(`tval') {
      if abs(`vval'-`tval') > 1e-12 local flag "`flag'/VECTOR_MISMATCH"
    }
    dis "PERMX|`tag'|`est'|e(perm_p_`est')=`sc'|table p=`tval'|e(perm_p)[main_`est']=`vval'|`flag'"
  }
  * the two near-identical distribution matrices
  local d1 = .
  local d2 = .
  cap local d1 = colsof(e(permdist))
  cap local d2 = colsof(e(permdist_b))
  dis "PERMDIST|`tag'|permdist cols=`d1'|permdist_b cols=`d2'|K=" colsof(e(b))
end

foreach S in "att" "att att_i" "att_itime" "att_i att_itime" "att att_i att_it att_itime" {
  cap noi wooldid y id t ttre, `S' cluster(clu) permute(shuffle) permreps(40) rseed(7)
  local rc = _rc
  dis "PERM_RC[`S']=`rc'"
  if `rc' == 0 {
    PERMX "`=subinstr("`S'"," ","_",.)'"
    ckc "PERM_`=subinstr("`S'"," ","_",.)'"
  }
}
* with an event study and subgroups
gen byte sg = mod(id,2)
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu) permute(shuffle) permreps(40) rseed(7)
dis "PERM_ES_RC=" _rc
if _rc == 0 {
  PERMX es
  ckc PERM_ES
}
dis "PERM_DONE"
