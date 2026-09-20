* CONVEX-HULL INVARIANT.  Every att-family estimand is a NONNEGATIVELY weighted average of
* the post-treatment cohort-by-relative-period cell effects in e(gteffects).  So it must lie
* between the smallest and the largest of them.  A point estimate outside that hull means
* some aggregation weight went negative (or the wrong cells were pooled).
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_hull
program define wdfz_hull, rclass
  syntax , TAG(string)
  capture matrix G = e(gteffects)
  if _rc {
    dis "HULL|`tag'|NO_GTEFFECTS"
    exit
  }
  local gc : colnames G
  local ce = 0
  local cr = 0
  local k = 0
  foreach c of local gc {
    local k = `k' + 1
    if "`c'" == "estimate" local ce = `k'
    if "`c'" == "relyr"    local cr = `k'
  }
  if `ce' == 0 {
    dis "HULL|`tag'|NOCOL estimate; colnames=`gc'"
    exit
  }
  local RG = rowsof(G)
  local lo = .
  local hi = .
  local npost 0
  forvalues r = 1/`RG' {
    local rel = .
    if `cr' > 0 local rel = G[`r',`cr']
    if `cr' > 0 & (missing(`rel') | `rel' < 0) continue
    local v = G[`r',`ce']
    if missing(`v') continue
    local npost = `npost' + 1
    if missing(`lo') | `v' < `lo' local lo = `v'
    if missing(`hi') | `v' > `hi' local hi = `v'
  }
  matrix M = e(wdidmainresults)
  local R = rowsof(M)
  local rn : rownames M
  local re : roweq M
  local bad ""
  local nchk 0
  forvalues r = 1/`R' {
    local ee : word `r' of `re'
    local nn : word `r' of `rn'
    if "`ee'" != "main" continue
    if !inlist("`nn'","att","att_it","att_i","att_itime") continue
    local v = M[`r',1]
    if missing(`v') continue
    local nchk = `nchk' + 1
    local tol = 1e-8 * max(1, abs(`hi'-`lo'))
    if `v' < `lo' - `tol' local bad "`bad' BELOW(`nn'=`v'<`lo')"
    if `v' > `hi' + `tol' local bad "`bad' ABOVE(`nn'=`v'>`hi')"
  }
  dis "HULL|`tag'|postcells=`npost'|lo=`lo'|hi=`hi'|checked=`nchk'|`bad'"
end

* one diagnostic print of the gteffects layout
use "`D'/ds_canon.dta", clear
qui wooldid y unit t ttre, att gteffects cluster(region)
dis "GTCOLS: `: colnames e(gteffects)'"
dis "GTDIM: " rowsof(e(gteffects)) " x " colsof(e(gteffects))
matrix list e(gteffects)

local SPECS `"att att_i att_it att_itime cluster(region) gteffects"'
foreach ds in canon unbal rotate rcs biennial fewG10 pw_survey miss_y nevergrp longT hetero_unit onsetgap clustonset sg_gaps {
  foreach ex in "PLAIN" "robust" "cluster(unit)" "AWEIGHT" "esprelength(2) espostlength(2)" "controls(x1 x2)" "fe(region)" "IFCOND" "anticipation(1)" "subgroups(sg2)" {
    local wc ""
    local ifc ""
    local oo "`ex'"
    if "`ex'" == "PLAIN" local oo ""
    if "`ex'" == "AWEIGHT" {
      local wc "[aw=w]"
      local oo ""
    }
    if "`ex'" == "IFCOND" {
      local ifc "if region<=3"
      local oo ""
    }
    use "`D'/ds_`ds'.dta", clear
    ereturn clear
    capture qui wooldid y unit t ttre `ifc' `wc', att att_i att_it att_itime cluster(region) gteffects `oo'
    if _rc == 0 wdfz_hull , tag(`ds'/`ex')
    else dis "HULL|`ds'/`ex'|RC=`=_rc'"
  }
}
dis "HULL_DONE"
