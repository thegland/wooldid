* amese_structure.do — deterministic STRUCTURE PASS over the 28-case AME-SE grid.
* One run per case on a fixed seed: validates that each case's intended absorb
* structure materializes (verbose echoes the internal regression), then records rc,
* V_es_ame rank/eigenvalues, the displayed ame SE pattern (how many distinct values),
* and e(es_ndegen_ame). Output: grep-friendly STRUCT| lines.
* args: 1 = path to wooldid.ado (file or dir)
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local ADO `1'
if `"`ADO'"' == "" local ADO wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
do "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh"

local CASES ///
  d1_basic_ci d1_basic_cu d1_basic_rb d1_bycoh_ci d1_bycoh_cu d1_bycoh_rb ///
  d1_bytime_ci d1_bytime_cu d1_bytime_rb ///
  d2_basic_ci d2_basic_cu d2_basic_rb d2_bycoh_ci d2_bycoh_cu d2_bycoh_rb ///
  d2_bytime_ci d2_bytime_cu d2_bytime_rb ///
  d2_within_ci d2_within_cu d2_within_rb d1_within_ci d1_within_cu d1_within_rb ///
  d3_basic_ci d3_basic_cu d3_basic_rb d1c_basic_ci

foreach case of local CASES {
  tokenize "`case'", parse("_")
  local dc `1'
  local cc `3'
  local vc `5'
  set seed 60000001
  amese_mkdata `dc'
  amese_opts `cc' `vc'
  local opts `r(opts)'
  dis _n as result "STRUCT|`case'|CALL|wooldid y i t ttre, `opts'"
  cap noi wooldid y i t ttre, `opts' verbose
  local rc = _rc
  if `rc' != 0 {
    dis "STRUCT|`case'|RC|`rc'"
    continue
  }
  dis "STRUCT|`case'|RC|0"
  local nd = e(es_ndegen_ame)
  tempname V W
  cap matrix `V' = e(V_es_ame)
  if _rc != 0 dis "STRUCT|`case'|VES|absent|ndegen|`nd'"
  else {
    mata: V = st_matrix(st_local("V")); ev = symeigenvalues(V); st_local("evmax", strofreal(max(ev), "%12.0g")); st_local("evmin", strofreal(min(ev), "%12.0g")); st_local("rk", strofreal(sum(ev :> 1e-10*max(ev)), "%3.0f"))
    dis "STRUCT|`case'|VES|`=colsof(`V')'x`=colsof(`V')'|rank|`=trim("`rk'")'|evmax|`evmax'|evmin|`evmin'|ndegen|`nd'"
  }
  cap matrix `W' = e(wdidmainresults)
  if _rc == 0 {
    local ses
    foreach rn in "main:ame" "pretreatment:ame" "ES_ame:pre5" "ES_ame:pre4" "ES_ame:pre3" "ES_ame:pre2" "ES_ame:post0" "ES_ame:post1" "ES_ame:post2" "ES_ame:post3" {
      local r = rownumb(`W', "`rn'")
      if `r' == . local ses "`ses' ."
      else {
        local v : di %11.0g `W'[`r', 2]
        local ses "`ses' `=trim("`v'")'"
      }
    }
    dis "STRUCT|`case'|AMESE|`ses'"
  }
}
dis _n "STRUCTURE_PASS_DONE"
exit, clear
