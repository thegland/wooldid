*=====================================================================
* p1_orphans.do -- does an option WITHOUT its parent option change
*                  anything, and does it say anything?
* Self-contained.  Generates its own data (gen_core.do), runs each
* option twice (base / base+orphan), captures a NAMED log per run and
* an exact-hex dump of the FULL e() per run.  The shell diffs them.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set maxvar 20000

* ---- HAZARD (b): assert we are on the pinned build, not the v1.3 in PLUS ----
which wooldid
dis "HEXFMT=" string(1.5,"%21x")

* ---- exact full-e() dumper -----------------------------------------
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
    local cn : colnames `M'
    file write `fh' "X `x' dim `nr' `nc' rows(`rn') cols(`cn')" _n
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

* ---- one run: named log + e() dump ---------------------------------
cap prog drop runcase
prog def runcase
  syntax , nm(string) opts(string asis) [yv(string) wt(string asis)]
  if "`yv'" == "" local yv ync
  cap log close cx
  qui log using "L_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "E_`nm'.txt"
  dis "RCLINE `nm' = `rc'"
end

* ---- data ----------------------------------------------------------
local gN      80
local gT      10
local gt0     2000
local gonsets 2004 2006 2008
local gysd    0.30
local gdose   lattice
local gwtype  unit
local gwspread 2
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"

* a pre-made lattice variable for the 1-token lattice() form
qui egen latvar = cut(dose), group(4)
qui replace latvar = 0 if latvar >= .
qui count
dis "NOBS=" r(N)
qui levelsof unit, local(uu)
dis "NUNITS=" wordcount("`uu'")
qui save "pdat.dta", replace

* ---- clear away any stale plot files from a previous run ------------
foreach z in zz1 zz2 zz3 zz4 {
  cap erase "`z'_ES.dta"
}

*=====================================================================
* THE PAIRS.  cN_b = base, cN_o = base + the orphan option
*=====================================================================
local NC 16

local c1_lab  saveplots_bare
local c1_b    att robust
local c1_o    att robust saveplots(zz1)

local c2_lab  saveplots_makeplots_nohorizons
local c2_b    att robust makeplots
local c2_o    att robust makeplots saveplots(zz2)

local c3_lab  cfxplottypes
local c3_b    att robust contreat(dose)
local c3_o    att robust contreat(dose) cfxplottypes(post 2)

local c4_lab  cfxplottypes_garbage
local c4_b    att robust contreat(dose)
local c4_o    att robust contreat(dose) cfxplottypes(psot 1.5 nosuchvar)

local c5_lab  cfxplotnose
local c5_b    att robust contreat(dose)
local c5_o    att robust contreat(dose) cfxplotnose

local c6_lab  makecfxplotsbysg
local c6_b    att robust contreat(dose) subgroups(sg2)
local c6_o    att robust contreat(dose) subgroups(sg2) makecfxplotsbysg

local c7_lab  latticeignoreweights_nolattice
local c7_b    att robust contreat(dose)
local c7_o    att robust contreat(dose) latticeignoreweights

local c8_lab  latticeignoreweights_premadelattice
local c8_b    att robust contreat(dose) lattice(latvar)
local c8_o    att robust contreat(dose) lattice(latvar) latticeignoreweights

local c9_lab  ccnoabsorb_nocontreat
local c9_b    att robust
local c9_o    att robust ccnoabsorb

local c10_lab contreatpoly_nocontreat
local c10_b   att robust
local c10_o   att robust contreatpoly(3)

local c11_lab contreatelasticitytype_nocontreat
local c11_b   att robust
local c11_o   att robust contreatelasticitytype(eyex)

local c12_lab contreatcontrols_nocontreat
local c12_b   att robust
local c12_o   att robust contreatcontrols(byTreated)

local c13_lab contreatcontrols_garbage_nocontreat
local c13_b   att robust
local c13_o   att robust contreatcontrols(zzzNotAThing)

local c14_lab suppresssg_nosubgroups
local c14_b   att robust
local c14_o   att robust suppresssgprimaryeffects

local c15_lab simultaneousbands_nohorizons
local c15_b   att robust
local c15_o   att robust simultaneousbands

local c16_lab histogram_att_only
local c16_b   att robust
local c16_o   att robust histogramcohorteffects

forval k = 1/`NC' {
  use "pdat.dta", clear
  dis "==== CASE `k' : `c`k'_lab' ===="
  runcase , nm("c`k'_b") opts(`c`k'_b') wt([aw=w])
  use "pdat.dta", clear
  runcase , nm("c`k'_o") opts(`c`k'_o') wt([aw=w])
  dis "CASELAB `k' = `c`k'_lab'"
}

*=====================================================================
* Control: cfxplottypes(garbage) WITH its parent must error (proving
* the validation exists and is skipped only because the parent is absent)
*=====================================================================
use "pdat.dta", clear
cap noi wooldid ync unit t ttre [aw=w], att robust contreat(dose) makecfxplots cfxplottypes(psot)
dis "CTRL_cfxgarbage_withparent_RC=" _rc

use "pdat.dta", clear
cap noi wooldid ync unit t ttre [aw=w], att robust contreat(dose) contreatcontrols(zzzNotAThing)
dis "CTRL_ccgarbage_withparent_RC=" _rc

*=====================================================================
* FILE EVIDENCE for saveplots without makeplots
*=====================================================================
foreach z in zz1 zz2 {
  foreach sfx in _ES.dta _cfx.dta _att.gph _hist_att_i.gph {
    cap confirm file "`z'`sfx'"
    dis "FILE `z'`sfx' RC=" _rc
  }
}

dis "PROBE1_COMPLETE"
