* LENS 2 part 5, two impact facts left open:
*  (a) on the draws where the DEFAULT (reghdfe) path does NOT zero out under
*      lattice(reg3), what ame SE does it PRINT, and is that number right?
*      (reghdfe-robust and regress-robust are the same sandwich on the same fit,
*       so the prediction is ~0.156 against a realized sampling SD of ~0.554)
*  (b) does plain bootstrap(standard) -- the remedy wooldid's own caution names --
*      actually deliver numbers on a zero-V spec?
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333

capture program drop mkdata
program define mkdata
  clear
  set obs 12
  gen int st = _n
  gen int ncty = 5 + mod(_n*7, 8)
  expand ncty
  gen long cty = _n
  gen byte reg3 = 1 + mod(_n,3)
  expand 10
  bysort cty: gen int year = _n
  gen double ttre = .
  qui replace ttre = 4 if st <= 3
  qui replace ttre = 6 if st > 3 & st <= 6
  qui replace ttre = 8 if st > 6 & st <= 8
  gen byte post = (ttre < . & year >= ttre)
  bysort cty (year): gen double dbase = 2*runiform() if _n == 1
  bysort cty (year): replace dbase = dbase[1]
  gen double dose = dbase + rnormal(0, 0.5)
  gen double aunit = 0
  bysort cty (year): replace aunit = 3*runiform() if _n == 1
  bysort cty (year): replace aunit = aunit[1]
  gen double y = 5 + aunit + 0.25*year + 0.7*dose + 0.4*dose*post + 0.5*post + rnormal(0,1)
end

dis as result "=== (a) what the DEFAULT path PRINTS when it does not zero out ==="
tempname P
postfile `P' double(ameb amese attb attse fired) using L5mc , replace
forval r = 1/14 {
  mkdata
  local ab = .
  local as = .
  local tb = .
  local ts = .
  local fired = .
  cap qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(reg3) robust
  if _rc == 0 {
    mat M = e(wdidmainresults)
    local tb = M[1,1]
    local ts = M[1,2]
    local ab = M[2,1]
    local as = M[2,2]
    local fired = (M[1,2] >= .)
  }
  post `P' (`ab') (`as') (`tb') (`ts') (`fired')
  dis as result "L5 draw`r' fired=`fired'  DEFAULT att b=" %10.5f `tb' " se=" %10.5f `ts' "   ame b=" %10.5f `ab' " se=" %10.5f `as'
}
postclose `P'
use L5mc , clear
qui sum ameb
dis as result "L5 ame across all draws: mean=" %10.5f r(mean) " realized MC SD=" %10.5f r(sd) "   (TRUTH 0.400)"
qui sum amese if fired == 0
dis as result "L5 ame SE the DEFAULT path PRINTED on the non-zeroing draws: n=" r(N) " mean=" %10.5f r(mean) " min=" %10.5f r(min) " max=" %10.5f r(max)
qui sum attse if fired == 0
dis as result "L5 att SE the DEFAULT path PRINTED on the non-zeroing draws: n=" r(N) " mean=" %10.5f r(mean)
qui sum attb
dis as result "L5 att across all draws: mean=" %10.5f r(mean) " realized MC SD=" %10.5f r(sd)
qui count if fired == 1
dis as result "L5 zero-V incidence lattice(reg3): " r(N) "/14"
qui count if fired == 0 & abs(ameb - 0.4) <= 1.96*amese
qui count if fired == 0
local nf = r(N)
qui count if fired == 0 & abs(ameb - 0.4) <= 1.96*amese
dis as result "L5 cover95 of the PRINTED default-path ame interval: " r(N) "/" `nf'

dis as result "=== (b) plain bootstrap(standard): does the named remedy deliver? ==="
mkdata
cap noi wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust bootstrap(standard) reps(25) rseed(99)
dis "L5b_RC=" _rc
dis "DONE_L5"
