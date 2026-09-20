* Is the reported interval actually at level()?  For every analytic row with a finite SE,
* (ub-lb)/(2*se) must equal the critical value implied by e(df_r) and level().
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_civ
program define wdfz_civ, rclass
  syntax , TAG(string) LEV(real)
  local alpha = (100-`lev')/200
  local df = e(df_r)
  local tcrit = invttail(`df', `alpha')
  local zcrit = invnormal(1-`alpha')
  matrix M = e(wdidmainresults)
  local R = rowsof(M)
  local rn : rownames M
  local re : roweq M
  local badt 0
  local badz 0
  local n 0
  local detail ""
  forvalues r = 1/`R' {
    local se = M[`r',2]
    local lb = M[`r',5]
    local ub = M[`r',6]
    if missing(`se') | missing(`lb') | missing(`ub') continue
    if `se' <= 0 continue
    local n = `n' + 1
    local hw = (`ub'-`lb')/(2*`se')
    if abs(`hw' - `tcrit') > 1e-7*max(1,`tcrit') {
      local badt = `badt' + 1
      if abs(`hw' - `zcrit') < 1e-7*max(1,`zcrit') local badz = `badz' + 1
      local nn : word `r' of `rn'
      local ee : word `r' of `re'
      if `badt' <= 4 local detail "`detail' (`ee':`nn' hw=`hw')"
    }
  }
  dis "CIV|`tag'|lev=`lev'|df=`df'|tcrit=`tcrit'|zcrit=`zcrit'|rows=`n'|mismatch=`badt'|matches_z=`badz'|`detail'"
end

foreach ds in canon unbal rcs fewG10 clustonset longT {
  foreach v in "" "robust" "cluster(region)" "cluster(unit)" "cluster(cohort)" {
    foreach L in 80 90 95 99 {
      use "`D'/ds_`ds'.dta", clear
      ereturn clear
      capture qui wooldid y unit t ttre, att att_i att_it att_itime `v' esprelength(2) espostlength(2) level(`L')
      if _rc == 0 wdfz_civ , tag(`ds'-`=cond("`v'"=="","iid","`v'")') lev(`L')
      else dis "CIV|`ds'-`v'|lev=`L'|RC=`=_rc'"
    }
  }
}
dis "CIV_DONE"
