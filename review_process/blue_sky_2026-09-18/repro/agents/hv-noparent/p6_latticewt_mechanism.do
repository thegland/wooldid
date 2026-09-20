*=====================================================================
* p6_latticewt_mechanism.do -- WHY did latticeignoreweights change
* nothing on lattice(byTreated 4) with aweights, when the weighted and
* unweighted quartiles of the dose demonstrably differ?
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
    file write `fh' "X `x' dim `nr' `nc'" _n
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
  qui log using "L6_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "E6_`nm'.txt"
  dis "RCLINE `nm' = `rc'"
end

use "pdat5.dta", clear
gen byte evertreated = (ttre != .)

*---- 1. does gquantiles honour the weight AT ALL on this data? -------
qui gquantiles dose, _pctile nquantiles(4)
dis "GQ_UNW  p25=" r(r1) " p50=" r(r2) " p75=" r(r3)
qui gquantiles dose [aw=wgt], _pctile nquantiles(4)
dis "GQ_WTD  p25=" r(r1) " p50=" r(r2) " p75=" r(r3)

*---- 2. xtile WITHOUT by() --------------------------------------------
qui gquantiles xa = dose, xtile nquantiles(4)
qui gquantiles xb = dose [aw=wgt], xtile nquantiles(4)
qui count if xa != xb
dis "XTILE_NOBY_DIFFROWS=" r(N)

*---- 3. xtile WITH by(evertreated) -- exactly what wooldid does for
*        lattice(byTreated k) ------------------------------------------
qui gquantiles ya = dose, xtile nquantiles(4) by(evertreated)
qui gquantiles yb = dose [aw=wgt], xtile nquantiles(4) by(evertreated)
qui count if ya != yb
dis "XTILE_BY_DIFFROWS=" r(N)
qui tab ya yb

*---- 4. the "base" branch (_pctile on the baseTreated domain) ---------
qui gquantiles dose if evertreated == 1, _pctile nquantiles(4)
dis "BASE_UNW p25=" r(r1) " p50=" r(r2) " p75=" r(r3)
qui gquantiles dose [aw=wgt] if evertreated == 1, _pctile nquantiles(4)
dis "BASE_WTD p25=" r(r1) " p50=" r(r2) " p75=" r(r3)

*---- 5. wooldid with the OTHER lattice types, +/- latticeignoreweights
local CC contreat(dose) contreatcontrols(untreatedZero)
foreach lt in all baseTreated byTime byCohortByTime byTreatedByTime {
  use "pdat5.dta", clear
  runcase , nm(m_`lt'_w)  opts(att robust `CC' lattice(`lt' 4)) wt([aw=wgt])
  use "pdat5.dta", clear
  runcase , nm(m_`lt'_i)  opts(att robust `CC' lattice(`lt' 4) latticeignoreweights) wt([aw=wgt])
}
dis "PROBE6_COMPLETE"
