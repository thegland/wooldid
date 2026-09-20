*=====================================================================
* p2_contreat_and_label.do
*  (A) redo the contreat-family orphan cases (p1 failed them because
*      contreat() itself demands contreatcontrols())
*  (B) DEEP DIVE: contreatelasticitytype() WITHOUT contreat() is not
*      inert -- it relabels the results table.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set maxvar 20000
which wooldid

cap prog drop edump
prog def edump
  args fn
  tempname fh
  file open `fh' using "`fn'", write replace text
  local ss : e(scalars)
  foreach s of local ss {
    local v = string(e(`s'), "%21x")
    file write `fh' "S `s' = `v'" _n
  }
  local ms : e(macros)
  foreach m of local ms {
    file write `fh' `"M `m' = `e(`m')'"' _n
  }
  local xs : e(matrices)
  foreach x of local xs {
    tempname M
    matrix `M' = e(`x')
    local nr = rowsof(`M')
    local nc = colsof(`M')
    local rn : rownames `M'
    file write `fh' "X `x' dim `nr' `nc' rows(`rn')" _n
    forval r = 1/`nr' {
      local line
      forval c = 1/`nc' {
        local line `line' `= string(`M'[`r',`c'], "%21x")'
      }
      file write `fh' "X `x' `r':`line'" _n
    }
  }
  file close `fh'
end

cap prog drop runcase
prog def runcase
  syntax , nm(string) opts(string asis) [yv(string) wt(string asis)]
  if "`yv'" == "" local yv ync
  cap log close cx
  qui log using "L2_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "E2_`nm'.txt"
  dis "RCLINE `nm' = `rc'"
end

local gN      80
local gT      10
local gt0     2000
local gonsets 2004 2006 2008
local gysd    0.30
local gdose   lattice
local gwtype  unit
local gwspread 2
local gpois   1
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
qui egen latvar = cut(dose), group(4)
qui replace latvar = 0 if latvar >= .
qui save "pdat2.dta", replace

*---------------------------------------------------------------------
* (A) contreat-family orphans, now with the mandatory contreatcontrols()
*---------------------------------------------------------------------
local CC contreat(dose) contreatcontrols(byTreated)

local d1_lab  cfxplottypes_noparent
local d1_b    att robust `CC'
local d1_o    att robust `CC' cfxplottypes(post 2)

local d2_lab  cfxplottypes_GARBAGE_noparent
local d2_b    att robust `CC'
local d2_o    att robust `CC' cfxplottypes(psot 1.5 nosuchvar)

local d3_lab  cfxplotnose_noparent
local d3_b    att robust `CC'
local d3_o    att robust `CC' cfxplotnose

local d4_lab  makecfxplotsbysg_noparent
local d4_b    att robust `CC' subgroups(sg2)
local d4_o    att robust `CC' subgroups(sg2) makecfxplotsbysg

local d5_lab  latticeignoreweights_nolattice
local d5_b    att robust `CC'
local d5_o    att robust `CC' latticeignoreweights

local d6_lab  latticeignoreweights_premadelattice
local d6_b    att robust `CC' lattice(latvar)
local d6_o    att robust `CC' lattice(latvar) latticeignoreweights

local d7_lab  latticeignoreweights_CONTROL_realLattice
local d7_b    att robust `CC' lattice(byTreated 4)
local d7_o    att robust `CC' lattice(byTreated 4) latticeignoreweights

local d8_lab  ccnoabsorb_WITH_contreat_CONTROL
local d8_b    att robust `CC'
local d8_o    att robust `CC' ccnoabsorb

local d9_lab  contreatpoly_WITH_contreat_CONTROL
local d9_b    att robust `CC'
local d9_o    att robust `CC' contreatpoly(3)

forval k = 1/9 {
  use "pdat2.dta", clear
  dis "==== D`k' : `d`k'_lab' ===="
  runcase , nm("d`k'_b") opts(`d`k'_b') wt([aw=w])
  use "pdat2.dta", clear
  runcase , nm("d`k'_o") opts(`d`k'_o') wt([aw=w])
  dis "DLAB `k' = `d`k'_lab'"
}

*---------------------------------------------------------------------
* (B) the elasticity RELABEL on a run with NO continuous treatment
*---------------------------------------------------------------------
local e1_lab  base_binary
local e1_o    att att_i robust
local e2_lab  cet_eyex_noparent
local e2_o    att att_i robust contreatelasticitytype(eyex)
local e3_lab  cet_eydx_noparent
local e3_o    att att_i robust contreatelasticitytype(eydx)
local e4_lab  cet_dyex_noparent
local e4_o    att att_i robust contreatelasticitytype(dyex)
local e5_lab  semielasticity_binary_CONTROL
local e5_o    att att_i robust semielasticity
local e6_lab  cet_eyex_noparent_withES
local e6_o    att att_i robust esprelength(2) espostlength(2) contreatelasticitytype(eyex)
local e7_lab  base_withES
local e7_o    att att_i robust esprelength(2) espostlength(2)

forval k = 1/7 {
  use "pdat2.dta", clear
  dis "==== E`k' : `e`k'_lab' ===="
  runcase , nm("e`k'") opts(`e`k'_o') wt([aw=w])
  dis "ELAB `k' = `e`k'_lab'"
}

* poisson leg: no contreat, elasticity type requested
use "pdat2.dta", clear
dis "==== E8 : poisson_cet_eyex_noparent ===="
runcase , nm("e8") opts(att robust poisson contreatelasticitytype(eyex)) yv(ycount)
use "pdat2.dta", clear
dis "==== E9 : poisson_base ===="
runcase , nm("e9") opts(att robust poisson) yv(ycount)

* --- show the visible header lines side by side, unwrapped as far as possible
set linesize 250
dis "%%%% HEADER COMPARISON %%%%"
use "pdat2.dta", clear
cap noi wooldid ync unit t ttre [aw=w], att att_i robust
dis "HDR_RC_base=" _rc
use "pdat2.dta", clear
cap noi wooldid ync unit t ttre [aw=w], att att_i robust contreatelasticitytype(eyex)
dis "HDR_RC_eyex=" _rc

dis "PROBE2_COMPLETE"
