* honestdid_crosscheck.do — F2 -> HonestDiD integration check.
* On dta_binary: wooldid ES window + esfixedbaseperiod; feed e(b_es_att)/e(V_es_att) to
* honestdid exactly as the wooldid help file documents; assert rc 0 and that honestdid's
* "Original" row CI reproduces the wooldid post0 CI (est +/- invnormal(.975)*se — probe
* showed honestdid uses invnormal(.975), not literal 1.96) at 1e-6 relative. Also assert
* diag(e(V_es_att)) == reported ES SEs squared at 1e-12 (the F2 identity).
* Args: 1 = path to wooldid.ado (file or dir); 2 = suite data dir (contains dta_binary.dta).
* Run:  "<stata>" -e do honestdid_crosscheck.do <adopath> <datadir>
*    -> expect HONESTDID_CROSSCHECK_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 8642
set sortseed 8642

local ADO `1'
if `"`ADO'"' == "" local ADO /tmp/wdsuite_build/newfeat/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local DATA `2'
if `"`DATA'"' == "" local DATA /Users/hegland/Documents/Research/wooldid_dev/test/data

global pf 0

use `"`DATA'/dta_binary.dta"', clear

* ---- H0: the ES run (help-file example 10 configuration + cluster) ----
cap noi wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "H0 FAIL: wooldid ES run rc=" _rc
  global pf = $pf + 1
  dis as error "HONESTDID_CROSSCHECK_FAILS=$pf"
  exit, clear
}
tempname BES VES WMR
matrix `BES' = e(b_es_att)
matrix `VES' = e(V_es_att)
matrix `WMR' = e(wdidmainresults)
local npre = e(numpre_att)
local eslbl : colnames `BES'
local wmrlbl : rownames `WMR'
local kes = colsof(`BES')

* ---- H1: the F2 identity — diag(e(V_es_att)) == reported ES SEs squared (1e-12) ----
* the ES rows of e(wdidmainresults) are rows 3.., in the same order as e(b_es_att);
* guard the positional match by name before comparing
forval j = 1/`kes' {
  local esl : word `j' of `eslbl'
  local wml : word `=`j'+2' of `wmrlbl'
  if "`esl'" != "`wml'" {
    dis as error "H1 FAIL: row alignment broken at `j': b_es label [`esl'] vs wdidmainresults row [`wml']"
    global pf = $pf + 1
    continue
  }
  local sej = `WMR'[`=`j'+2', 2]
  if reldif(`VES'[`j', `j'], `sej'^2) > 1e-12 {
    dis as error "H1 FAIL at `esl': diag(V_es)=" %21.0g `VES'[`j',`j'] " vs se^2=" %21.0g `sej'^2 " (reldif " %12.4e reldif(`VES'[`j',`j'], `sej'^2) ")"
    global pf = $pf + 1
  }
  * bonus: the coefficient itself matches the reported estimate
  if reldif(`BES'[1, `j'], `WMR'[`=`j'+2', 1]) > 1e-12 {
    dis as error "H1b FAIL at `esl': b_es=" %21.0g `BES'[1,`j'] " vs reported " %21.0g `WMR'[`=`j'+2', 1]
    global pf = $pf + 1
  }
}

* ---- H2: honestdid consumes the matrices the way the sthlp documents (example 10) ----
matrix b = `BES'
matrix V = `VES'
cap noi honestdid, numpre(`npre') b(b) vcov(V) delta(rm) mvec(0.5(0.5)2) parallel(0)
if _rc != 0 {
  dis as error "H2 FAIL: honestdid rc=" _rc
  global pf = $pf + 1
  dis as error "HONESTDID_CROSSCHECK_FAILS=$pf"
  exit, clear
}
local hobj `s(HonestEventStudy)'
tempname HCI
cap noi mata: st_matrix(st_local("HCI"), `hobj'.CI)
local mrc = _rc
cap confirm matrix `HCI'
if `mrc' != 0 | _rc != 0 {
  dis as error "H2b FAIL: could not extract the CI matrix (mata rc=`mrc', confirm rc=" _rc ")"
  global pf = $pf + 1
  dis as error "HONESTDID_CROSSCHECK_FAILS=$pf"
  exit, clear
}
if rowsof(`HCI') != 5 | colsof(`HCI') != 3 {
  dis as error "H2b FAIL: CI matrix is " rowsof(`HCI') "x" colsof(`HCI') " (want 5x3)"
  global pf = $pf + 1
  dis as error "HONESTDID_CROSSCHECK_FAILS=$pf"
  exit, clear
}

* ---- H3: the "Original" row (row 1) reproduces the wooldid post0 CI at 1e-6 relative ----
* the first post-period entry of b_es_att is post0 (guard by label)
local p0 = `npre' + 1
local p0l : word `p0' of `eslbl'
if "`p0l'" != "post0" {
  dis as error "H3 FAIL: entry `p0' of b_es_att is [`p0l'], expected post0"
  global pf = $pf + 1
}
tempname EST SE LB UB
scalar `EST' = `BES'[1, `p0']
scalar `SE'  = sqrt(`VES'[`p0', `p0'])
scalar `LB'  = `EST' - invnormal(.975)*`SE'
scalar `UB'  = `EST' + invnormal(.975)*`SE'
if reldif(`HCI'[1, 2], `LB') > 1e-6 | reldif(`HCI'[1, 3], `UB') > 1e-6 {
  dis as error "H3 FAIL: Original CI (" `HCI'[1,2] ", " `HCI'[1,3] ") vs wooldid post0 CI (" `LB' ", " `UB' ")"
  dis as error "         reldif lb = " %12.4e reldif(`HCI'[1,2], `LB') ", ub = " %12.4e reldif(`HCI'[1,3], `UB')
  global pf = $pf + 1
}
else dis "H3 info: Original CI reldif lb = " %12.4e reldif(`HCI'[1,2], `LB') ", ub = " %12.4e reldif(`HCI'[1,3], `UB')

* ---- H4: the robust rows exist for each requested M (values only; no reproducibility bet) ----
forval r = 2/5 {
  if reldif(`HCI'[`r', 1], (`r'-1)*0.5) > 1e-12 {
    dis as error "H4 FAIL: row `r' M = " `HCI'[`r',1] " (expected " (`r'-1)*0.5 ")"
    global pf = $pf + 1
  }
  if missing(`HCI'[`r', 2]) | missing(`HCI'[`r', 3]) | `HCI'[`r', 2] >= `HCI'[`r', 3] {
    dis as error "H4 FAIL: row `r' CI (" `HCI'[`r',2] ", " `HCI'[`r',3] ") is empty or inverted"
    global pf = $pf + 1
  }
}

if $pf == 0 dis "HONESTDID_CROSSCHECK_ALL_PASS"
else dis as error "HONESTDID_CROSSCHECK_FAILS=$pf"
exit, clear
