*=====================================================================
* p5_cfx_and_latticewt.do
*  (A) the cfx-family orphan options on a contreat base that estimates
*  (B) EFFICACY of latticeignoreweights WITH its parent lattice(type k)
*      on a genuinely continuous dose whose weighted and unweighted
*      quartiles are known to differ
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
  qui log using "L5_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "E5_`nm'.txt"
  dis "RCLINE `nm' = `rc'"
end

*---------------------------------------------------------------------
* data: CONTINUOUS dose, zero off treatment, weights strongly
*       correlated with the dose so weighted != unweighted quantiles
*---------------------------------------------------------------------
clear
set obs 90
gen long unit = _n
expand 10
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2003 if unit <= 20
qui replace ttre = 2005 if unit > 20 & unit <= 40
qui replace ttre = 2007 if unit > 40 & unit <= 55
set seed 13579
by unit: gen double u_dose = runiform()*10 if _n == 1
by unit: replace u_dose = u_dose[1]
gen byte post = (ttre < . & t >= ttre)
gen double dose = cond(post, u_dose, 0)
* weight correlated with dose: heavy weight on the HIGH-dose units only
gen double wgt = cond(u_dose > 6, 40, 1)
gen byte sg2 = mod(unit, 2)
gen double ync = 0.08*unit + 0.04*(t-1999) + (0.5 + 0.15*u_dose)*post + rnormal()*0.25
qui save "pdat5.dta", replace

* ---- oracle: the quartile cut points of dose on the baseTreated domain,
*      weighted and unweighted.  If these differ, latticeignoreweights MUST
*      change the lattice, hence the cell definitions, hence the results.
qui _pctile dose if ttre < ., p(25 50 75)
dis "UNW_P25=" r(r1) " P50=" r(r2) " P75=" r(r3)
qui _pctile dose [aw=wgt] if ttre < ., p(25 50 75)
dis "WTD_P25=" r(r1) " P50=" r(r2) " P75=" r(r3)

local CC contreat(dose) contreatcontrols(untreatedZero)
local LT lattice(byTreated 4)

*---------------------------------------------------------------------
* (B) EFFICACY of latticeignoreweights WITH lattice(type k)
*---------------------------------------------------------------------
use "pdat5.dta", clear
runcase , nm(h1_lattice_wtd)    opts(att robust `CC' `LT') wt([aw=wgt])
use "pdat5.dta", clear
runcase , nm(h2_lattice_ignwt)  opts(att robust `CC' `LT' latticeignoreweights) wt([aw=wgt])
* control: an unweighted run's lattice must match the ignoreweights lattice
use "pdat5.dta", clear
runcase , nm(h3_lattice_noweight) opts(att robust `CC' `LT')

*---------------------------------------------------------------------
* (A) cfx-family orphans on a base that estimates
*---------------------------------------------------------------------
local BASE att robust `CC' `LT'

local k1_lab cfxplottypes_noparent
local k1_o   `BASE' cfxplottypes(post 2)
local k2_lab cfxplottypes_GARBAGE_noparent
local k2_o   `BASE' cfxplottypes(psot 1.5 nosuchvar)
local k3_lab cfxplotnose_noparent
local k3_o   `BASE' cfxplotnose
local k4_lab makecfxplotsbysg_noparent
local k4_o   `BASE' subgroups(sg2) makecfxplotsbysg
local k5_lab base_subgroups_forK4
local k5_o   `BASE' subgroups(sg2)

use "pdat5.dta", clear
runcase , nm(k0_base) opts(`BASE') wt([aw=wgt])
forval k = 1/5 {
  use "pdat5.dta", clear
  dis "==== K`k' : `k`k'_lab' ===="
  runcase , nm("k`k'") opts(`k`k'_o') wt([aw=wgt])
  dis "KLAB `k' = `k`k'_lab'"
}

* CONTROLS: the same tokens WITH makecfxplots
use "pdat5.dta", clear
runcase , nm(k6_makecfx_ok)      opts(`BASE' makecfxplots cfxplottypes(post 2))
use "pdat5.dta", clear
runcase , nm(k7_makecfx_garbage) opts(`BASE' makecfxplots cfxplottypes(psot 1.5 nosuchvar))
use "pdat5.dta", clear
runcase , nm(k8_makecfx_nose)    opts(`BASE' makecfxplots cfxplotnose)

* saveplots WITH its parent, to show the files DO appear when the parent is there
foreach z in ww1 ww2 {
  foreach sfx in _ES.dta _cfx.dta {
    cap erase "`z'`sfx'"
  }
}
use "pdat5.dta", clear
runcase , nm(k9_saveplots_withparent) opts(att robust esprelength(2) espostlength(2) makeplots saveplots(ww1)) wt([aw=wgt])
use "pdat5.dta", clear
runcase , nm(k10_saveplots_noparent)  opts(att robust esprelength(2) espostlength(2) saveplots(ww2)) wt([aw=wgt])
foreach z in ww1 ww2 {
  foreach sfx in _ES.dta _ES_att.gph {
    cap confirm file "`z'`sfx'"
    dis "FILE `z'`sfx' RC=" _rc
  }
}
cap ls
dis "PROBE5_COMPLETE"
