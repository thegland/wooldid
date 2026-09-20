* WEIGHTED-RECONSTRUCTION INVARIANT (new).  e(gteffects) ships, per cohort-by-relative-period
* cell, both the cell coefficient (`estimate') and the four aggregation weights
* (w_att, w_att_it, w_att_i, w_att_itime).  Each reported estimand must therefore equal
*     sum_cells w_e * estimate / sum_cells w_e
* over the post cells.  This is an EXACT identity on REALISTIC (noisy) data: it needs no
* oracle DGP.  It also implies the convex-hull property (min <= estimand <= max), which is
* checked alongside so that a negative aggregation weight is visible as such.
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_gtwt
program define wdfz_gtwt, rclass
  syntax , TAG(string)
  capture matrix G = e(gteffects)
  if _rc {
    dis "GTWT|`tag'|NO_GTEFFECTS"
    exit
  }
  local gc : colnames G
  local k = 0
  local cEST = 0
  local cREL = 0
  local cWATT = 0
  local cWIT = 0
  local cWI = 0
  local cWITIME = 0
  foreach c of local gc {
    local k = `k' + 1
    if "`c'" == "estimate"    local cEST = `k'
    if "`c'" == "relyr"       local cREL = `k'
    if "`c'" == "w_att"       local cWATT = `k'
    if "`c'" == "w_att_it"    local cWIT = `k'
    if "`c'" == "w_att_i"     local cWI = `k'
    if "`c'" == "w_att_itime" local cWITIME = `k'
  }
  if `cEST'==0 | `cWATT'==0 {
    dis "GTWT|`tag'|LAYOUT(`gc')"
    exit
  }
  local RG = rowsof(G)
  matrix M = e(wdidmainresults)
  local R = rowsof(M)
  local rn : rownames M
  local re : roweq M
  local out ""
  local nneg 0
  foreach e in att att_it att_i att_itime {
    local cw = `cWATT'
    if "`e'" == "att_it"    local cw = `cWIT'
    if "`e'" == "att_i"     local cw = `cWI'
    if "`e'" == "att_itime" local cw = `cWITIME'
    if `cw' == 0 continue
    * the reported value
    local rep = .
    forvalues r = 1/`R' {
      local ee : word `r' of `re'
      local nn : word `r' of `rn'
      if "`ee'" == "main" & "`nn'" == "`e'" local rep = M[`r',1]
    }
    if missing(`rep') continue
    * reconstruct over POST cells and over ALL cells
    local sP = 0
    local wP = 0
    local sA = 0
    local wA = 0
    local lo = .
    local hi = .
    forvalues r = 1/`RG' {
      local w = G[`r',`cw']
      local v = G[`r',`cEST']
      if missing(`w') | missing(`v') continue
      if `w' < 0 local nneg = `nneg' + 1
      local sA = `sA' + `w'*`v'
      local wA = `wA' + `w'
      local rel = 0
      if `cREL' > 0 local rel = G[`r',`cREL']
      if `rel' >= 0 {
        local sP = `sP' + `w'*`v'
        local wP = `wP' + `w'
        if `w' > 0 {
          if missing(`lo') | `v' < `lo' local lo = `v'
          if missing(`hi') | `v' > `hi' local hi = `v'
        }
      }
    }
    local pP = .
    local pA = .
    if `wP' != 0 local pP = `sP'/`wP'
    if `wA' != 0 local pA = `sA'/`wA'
    local dP = 99
    local dA = 99
    if !missing(`pP') local dP = abs(`pP'-`rep')/max(1,abs(`rep'))
    if !missing(`pA') local dA = abs(`pA'-`rep')/max(1,abs(`rep'))
    local best = min(`dP',`dA')
    local hull "in"
    if !missing(`lo') & `rep' < `lo' - 1e-8*max(1,abs(`hi'-`lo')) local hull "BELOW"
    if !missing(`hi') & `rep' > `hi' + 1e-8*max(1,abs(`hi'-`lo')) local hull "ABOVE"
    local out "`out' `e'(rep=`rep' post=`pP' all=`pA' bestrel=`best' hull=`hull')"
  }
  dis "GTWT|`tag'|negwts=`nneg'|`out'"
end

foreach ds in canon unbal rotate rcs biennial fewG10 pw_survey miss_y nevergrp longT hetero_unit onsetgap clustonset sg_gaps {
  foreach ex in "PLAIN" "robust" "cluster(unit)" "AWEIGHT" "esprelength(2) espostlength(2)" "controls(x1 x2)" "fe(region)" "IFCOND" "subgroups(sg2)" "nofastcollapse" {
    local wc ""
    local ifc ""
    local oo "`ex'"
    if "`ex'" == "PLAIN"   local oo ""
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
    if _rc == 0 wdfz_gtwt , tag(`ds'/`ex')
    else dis "GTWT|`ds'/`ex'|RC=`=_rc'"
  }
}
dis "GTWT_DONE"
