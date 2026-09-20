*=====================================================================
* p7_latticewt_controls.do -- does latticeignoreweights do ANYTHING?
* Uses a repeated-cross-section shape (4 rows per unit-period cell) so
* the lattice genuinely subdivides cells and cannot be absorbed by the
* cell dummies.  CONTROLS included that prove the lattice itself bites.
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
    if "`m'" != "cmdline" file write `fh' `"M `m' = `e(`m')'"' _n
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
  if "`yv'" == "" local yv y
  cap log close cx
  qui log using "L7_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "E7_`nm'.txt"
  dis "RCLINE `nm' = `rc'"
end

*---------------------------------------------------------------------
* repeated cross-section: 4 individuals per unit-period cell, dose
* varying WITHIN the cell, weights strongly correlated with the dose
*---------------------------------------------------------------------
clear
set obs 60
gen long unit = _n
expand 10
bysort unit: gen int t = 1999 + _n
expand 4
gen double ttre = .
qui replace ttre = 2003 if unit <= 18
qui replace ttre = 2005 if unit > 18 & unit <= 36
set seed 2468
gen byte post = (ttre < . & t >= ttre)
gen double d    = runiform()*10
gen double dose = cond(post, d, 0)
gen double wgt  = cond(d > 6, 40, 1)
gen double y = 0.07*unit + 0.03*(t-1999) + (0.5 + 0.15*d)*post + rnormal()*0.3
gen byte evertreated = (ttre != .)
qui save "pdat7.dta", replace

* ---- hand-built oracle: the byTreated lattice, weighted vs unweighted
qui gquantiles latw = dose [aw=wgt], xtile nquantiles(4) by(evertreated)
qui gquantiles latu = dose,          xtile nquantiles(4) by(evertreated)
qui count if latw != latu
dis "LATTICE_DIFFROWS=" r(N) " of " _N
qui count
dis "NOBS=" r(N)

local CC contreat(dose) contreatcontrols(untreatedZero)

* ---- A: weighted lattice (the default)
use "pdat7.dta", clear
runcase , nm(a_lat4_wtd)   opts(att robust `CC' lattice(byTreated 4)) wt([aw=wgt])
* ---- B: same + latticeignoreweights  <- THE TEST
use "pdat7.dta", clear
runcase , nm(b_lat4_ignw)  opts(att robust `CC' lattice(byTreated 4) latticeignoreweights) wt([aw=wgt])
* ---- C: CONTROL, coarser lattice: must differ from A if the lattice bites
use "pdat7.dta", clear
runcase , nm(c_lat2_wtd)   opts(att robust `CC' lattice(byTreated 2)) wt([aw=wgt])
* ---- D: CONTROL, no lattice at all
use "pdat7.dta", clear
runcase , nm(d_nolattice)  opts(att robust `CC') wt([aw=wgt])
* ---- E: CONTROL, the hand-built UNWEIGHTED lattice passed as a premade variable
*        (this is what latticeignoreweights is supposed to reproduce)
use "pdat7.dta", clear
qui gquantiles latu = dose, xtile nquantiles(4) by(evertreated)
qui replace latu = 0 if latu >= .
runcase , nm(e_premade_unw) opts(att robust `CC' lattice(latu)) wt([aw=wgt])
* ---- F: CONTROL, the hand-built WEIGHTED lattice passed as a premade variable
use "pdat7.dta", clear
qui gquantiles latw = dose [aw=wgt], xtile nquantiles(4) by(evertreated)
qui replace latw = 0 if latw >= .
runcase , nm(f_premade_wtd) opts(att robust `CC' lattice(latw)) wt([aw=wgt])
* ---- G: pweights instead of aweights
use "pdat7.dta", clear
runcase , nm(g_lat4_pw)     opts(att robust `CC' lattice(byTreated 4)) wt([pw=wgt])
use "pdat7.dta", clear
runcase , nm(h_lat4_pw_ignw) opts(att robust `CC' lattice(byTreated 4) latticeignoreweights) wt([pw=wgt])
* ---- I: the 'all' lattice type on this data
use "pdat7.dta", clear
runcase , nm(i_all4_wtd)     opts(att robust `CC' lattice(all 4)) wt([aw=wgt])
use "pdat7.dta", clear
runcase , nm(j_all4_ignw)    opts(att robust `CC' lattice(all 4) latticeignoreweights) wt([aw=wgt])

dis "PROBE7_COMPLETE"
