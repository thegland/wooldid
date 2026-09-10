* test_impvalues.do - the Ibragimov-Muller reference distribution in e() (bug sweep BS-10,
* author ruling 2026-09-07). Run from test/: stata-mp -e do test_impvalues.do
* -> expect IMPVALUES_ALL_PASS
*
* Under impvalues every printed t, p and interval comes from a per-row ttest on t(q-1),
* where q is the number of clusters contributing to that row. e(df_r) used to carry the
* stage-1 regression's residual df instead, which is always larger, so a hand-run lincom,
* test or esttab reported p-values too small and intervals too narrow against the table
* wooldid had just printed. e(df_r) now carries the MINIMUM IM df across the reported rows.
*
* I1  the printed p-values reproduce from e(df_r) when the df are common across rows
* I1b lincom agrees with the printed p and interval (the ruling's own sentence)
* I2  a wide window makes the df not common: e(imdf_min)/e(imdf_max) bracket them and
*     e(df_r) is the minimum, so the posted reference is conservative, not exact, for the
*     rows that used a larger q
* I2b the footer names that, and replay reproduces it
* I3  the analytic df survives as hidden e(df_r_analytic)
* I5  e(df_r) is never missing, including on a degeneracy fixture
clear all
set more off
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop scanIM
program define scanIM, rclass
  args lf needle
  * Stata wraps long display lines and continues them after "> "; rejoin before searching
  tempname fh
  local found 0
  local acc
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if substr(`"`macval(line)'"', 1, 2) == "> " {
      local acc `"`macval(acc)'`=substr(`"`macval(line)'"', 3, .)'"'
    }
    else {
      if strpos(`"`macval(acc)'"', `"`needle'"') != 0 local found 1
      local acc `"`macval(line)'"'
    }
    file read `fh' line
  }
  if strpos(`"`macval(acc)'"', `"`needle'"') != 0 local found 1
  file close `fh'
  return scalar found = `found'
end

* ---- I1 / I1b / I3: common degrees of freedom
use data/dta_binary.dta, clear
qui wooldid y i t ttre, att cluster(i) impvalues esprelength(2) espostlength(2)
local d1 = e(df_r)
local a1 = e(df_r_analytic)
local mn1 = e(imdf_min)
local mx1 = e(imdf_max)
matrix M1 = e(wdidmainresults)
local bad1 0
local n1 = rowsof(M1)
forval r = 1/`n1' {
  if M1[`r',2] < . & M1[`r',2] > 0 {
    local pr = 2*ttail(`d1', abs(M1[`r',1]/M1[`r',2]))
    if reldif(`pr', M1[`r',4]) > 1e-10 local ++bad1
  }
}
qui lincom main_att
local lp = r(p)
local llb = r(lb)
local lub = r(ub)
di "I1 BS-10: df_r=`d1' analytic=`a1' imdf=[`mn1',`mx1'] rows checked=`n1' mismatches=`bad1' | lincom p=" %12.4e `lp' " vs printed " %12.4e M1[1,4] " lb=" %10.7f `llb' " vs " %10.7f M1[1,5]
if `d1' != 49 | `bad1' != 0 | `mn1' != 49 | `mx1' != 49 {
  di as error "I1 FAIL: e(df_r) does not reproduce the printed impvalues p-values (df_r=`d1', `bad1' mismatched row(s))"
  local ++fails
}
if reldif(`lp', M1[1,4]) > 1e-10 | abs(`llb' - M1[1,5]) > 1e-8 | abs(`lub' - M1[1,6]) > 1e-8 {
  di as error "I1b FAIL: lincom does not agree with the printed table"
  local ++fails
}
if `a1' != 79 {
  di as error "I3 FAIL: the stage-1 residual df is not preserved as e(df_r_analytic) (got `a1', want 79)"
  local ++fails
}

* ---- I2 / I2b: the degrees of freedom are not common across rows
use data/dta_binary.dta, clear
qui log using impv_i2.txt, text replace name(iv)
cap noi wooldid y i t ttre, att cluster(i) impvalues esprelength(6) espostlength(6)
local rc2 = _rc
qui log close iv
local d2 = e(df_r)
local mn2 = e(imdf_min)
local mx2 = e(imdf_max)
matrix M2 = e(wdidmainresults)
scanIM impv_i2.txt "Ibragimov-Muller degrees of freedom are not common across the rows above"
local c2 = r(found)
cap erase impv_i2.txt
* rows that used the minimum q reproduce exactly under e(df_r); rows that used a larger q
* must NOT - that is what makes the posted value conservative rather than universal
local exact2 0
local loose2 0
forval r = 1/`=rowsof(M2)' {
  if M2[`r',2] < . & M2[`r',2] > 0 {
    local pr = 2*ttail(`d2', abs(M2[`r',1]/M2[`r',2]))
    if reldif(`pr', M2[`r',4]) < 1e-10 local ++exact2
    else if reldif(`pr', M2[`r',4]) > 1e-6 local ++loose2
  }
}
di "I2 BS-10: rc=`rc2' df_r=`d2' imdf=[`mn2',`mx2'] rows on the minimum=`exact2' rows on a larger q=`loose2' footer=`c2'"
if `rc2' != 0 | `mn2' != 24 | `mx2' != 49 | `d2' != 24 | `exact2' == 0 | `loose2' == 0 {
  di as error "I2 FAIL: the not-common-df case does not post the bracketing minimum (rc=`rc2' df_r=`d2' imdf=[`mn2',`mx2'] exact=`exact2' loose=`loose2')"
  local ++fails
}
if `c2' != 1 {
  di as error "I2b FAIL: the not-common-df footer did not print"
  local ++fails
}

* ---- I5: never missing, including where the IM rows degenerate
use suite/data/dta_s_imdegen.dta, clear
cap qui wooldid y i t ttre, att cluster(clus) impvalues esprelength(1) espostlength(1)
local rc5 = _rc
local d5 = e(df_r)
di "I5 BS-10: degeneracy fixture rc=`rc5' df_r=`d5'"
if `rc5' == 0 & (`d5' >= . | `d5' <= 0) {
  di as error "I5 FAIL: e(df_r) is missing or non-positive on the degeneracy fixture (got `d5')"
  local ++fails
}

if `fails' == 0 di as result "IMPVALUES_ALL_PASS"
else di as error "IMPVALUES_FAILS=`fails'"
