clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/assert.do"

local specs `" "att robust" "att att_i att_it att_itime cluster(region) esprelength(3) espostlength(3)" "att cluster(unit) subgroups(sg2) jointtests" "att robust contreat(dose) contreatcontrols(basic)" "att cluster(region) permute(shuffle) permreps(10)" "att cluster(region) bootstrap(standard) reps(10) rseed(7)" "att robust esprelength(2) espostlength(2) makeplots" "'"
local k 0
foreach sp of local specs {
  local k = `k' + 1
  use "`D'/ds_canon.dta", clear
  wdfz_snap
  local n0 `r(nobs)'
  local v0 `r(nvar)'
  local s0 `r(sort)'
  local g0 `r(sig)'
  timer clear 1
  timer on 1
  ereturn clear
  capture noisily quietly wooldid y unit t ttre, `sp'
  local rc = _rc
  timer off 1
  qui timer list 1
  dis "TIME|`k'|`=r(t1)'|rc=`rc'|`sp'"
  if `rc' == 0 {
    wdfz_assert, id(`k') nobs0(`n0') nvar0(`v0') sort0(`s0') sig0(`g0')
    dis "VIOL|`k'|`r(viol)'"
  }
}
