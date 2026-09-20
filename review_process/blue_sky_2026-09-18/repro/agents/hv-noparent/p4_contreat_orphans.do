*=====================================================================
* p4_contreat_orphans.do
*  (A) the contreat-family orphan options, using contreatcontrols(
*      untreatedZero) which the generator's dose satisfies exactly
*  (B) minimal hand-built repro of the semielasticity failure, and a
*      test for whether it can fail SILENTLY (partial negativity)
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
  qui log using "L4_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "E4_`nm'.txt"
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
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
qui egen latvar = cut(dose), group(4)
qui replace latvar = 0 if latvar >= .
qui save "pdat4.dta", replace
qui sum dose if ttre < . & t >= ttre
dis "DOSE treated mean=" r(mean) " min=" r(min) " max=" r(max)
qui count if dose != 0 & (ttre >= . | t < ttre)
dis "NONZERO_DOSE_OFFTREAT=" r(N)

local CC contreat(dose) contreatcontrols(untreatedZero)

local f1_lab  cfxplottypes_noparent
local f1_b    att robust `CC'
local f1_o    att robust `CC' cfxplottypes(post 2)

local f2_lab  cfxplottypes_GARBAGE_noparent
local f2_b    att robust `CC'
local f2_o    att robust `CC' cfxplottypes(psot 1.5 nosuchvar)

local f3_lab  cfxplotnose_noparent
local f3_b    att robust `CC'
local f3_o    att robust `CC' cfxplotnose

local f4_lab  makecfxplotsbysg_noparent
local f4_b    att robust `CC' subgroups(sg2)
local f4_o    att robust `CC' subgroups(sg2) makecfxplotsbysg

local f5_lab  latticeignoreweights_NOlattice
local f5_b    att robust `CC'
local f5_o    att robust `CC' latticeignoreweights

local f6_lab  latticeignoreweights_PREMADElattice
local f6_b    att robust `CC' lattice(latvar)
local f6_o    att robust `CC' lattice(latvar) latticeignoreweights

local f7_lab  latticeignoreweights_CONTROL_generatedLattice
local f7_b    att robust `CC' lattice(byTreated 4)
local f7_o    att robust `CC' lattice(byTreated 4) latticeignoreweights

local f8_lab  saveplots_with_makecfxplots_CONTROL
local f8_b    att robust `CC' makecfxplots
local f8_o    att robust `CC' makecfxplots saveplots(qq1)

forval k = 1/8 {
  use "pdat4.dta", clear
  dis "==== F`k' : `f`k'_lab' ===="
  runcase , nm("f`k'_b") opts(`f`k'_b') wt([aw=w])
  use "pdat4.dta", clear
  runcase , nm("f`k'_o") opts(`f`k'_o') wt([aw=w])
  dis "FLAB `k' = `f`k'_lab'"
}
foreach sfx in _cfx.dta _cfxatt_cfxpost.gph {
  cap confirm file "qq1`sfx'"
  dis "FILE qq1`sfx' RC=" _rc
}

*=====================================================================
* (B) semielasticity: minimal hand-built repro + silent-failure test
*=====================================================================
cap prog drop tryit
prog def tryit
  syntax , nm(string) opts(string asis) [yv(string) wt(string asis)]
  if "`yv'" == "" local yv yy
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  dis "SENT `nm' RC=" _rc
  if _rc == 0 dis "SENT `nm' ATT=" string(el(e(wdidmainresults),1,1),"%21x")
end

clear
set obs 40
gen long unit = _n
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2003 if unit <= 12
qui replace ttre = 2005 if unit > 12 & unit <= 24
set seed 987654
gen double aw = 0.5 + runiform()
* yy is strictly positive; yneg crosses zero; ymostlypos is >0 except a few units
gen double base = 0.1*unit + 0.05*(t-1999) + 0.8*(ttre < . & t >= ttre) + rnormal()*0.2
gen double yy   = base + 5
gen double yneg = base - 2
gen double ymostlypos = base + 5
qui replace ymostlypos = ymostlypos - 9 if inlist(unit, 3, 30)
qui sum yy
dis "YY  min=" r(min)
qui sum yneg
dis "YNEG min=" r(min)
qui count if ymostlypos <= 0
dis "YMOSTLYPOS nonpositive rows=" r(N)
qui save "pdat4b.dta", replace

use "pdat4b.dta", clear
tryit , nm(g1_pos_noweight)  opts(att robust semielasticity) yv(yy)
use "pdat4b.dta", clear
tryit , nm(g2_pos_aw)        opts(att robust semielasticity) yv(yy) wt([aw=aw])
use "pdat4b.dta", clear
tryit , nm(g3_neg_noweight)  opts(att robust semielasticity) yv(yneg)
use "pdat4b.dta", clear
tryit , nm(g4_neg_aw)        opts(att robust semielasticity) yv(yneg) wt([aw=aw])
use "pdat4b.dta", clear
tryit , nm(g5_neg_pw)        opts(att robust semielasticity) yv(yneg) wt([pw=aw])
use "pdat4b.dta", clear
qui gen int fw = 1 + mod(unit,3)
tryit , nm(g6_neg_fw)        opts(att robust semielasticity) yv(yneg) wt([fw=fw])
use "pdat4b.dta", clear
tryit , nm(g7_mostlypos_noweight) opts(att robust semielasticity) yv(ymostlypos)
use "pdat4b.dta", clear
tryit , nm(g8_mostlypos_aw)  opts(att robust semielasticity) yv(ymostlypos) wt([aw=aw])
* controls: same runs without semielasticity
use "pdat4b.dta", clear
tryit , nm(g9_neg_aw_control) opts(att robust) yv(yneg) wt([aw=aw])
use "pdat4b.dta", clear
tryit , nm(g10_mostlypos_control) opts(att robust) yv(ymostlypos)
* and the contreat elasticity twin, which the help puts in the same family
use "pdat4b.dta", clear
qui gen double dose2 = cond(ttre < . & t >= ttre, 1 + mod(unit,4), 0)
tryit , nm(g11_cet_eyex_neg_aw) opts(att robust contreat(dose2) contreatcontrols(untreatedZero) contreatelasticitytype(eyex)) yv(yneg) wt([aw=aw])

dis "PROBE4_COMPLETE"
