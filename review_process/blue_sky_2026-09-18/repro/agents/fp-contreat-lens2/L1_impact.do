* ============================================================================
* LENS 2 (realism & impact) on fp-contreat-B1.
* Question NOT "does the engine split exist" (prober+skeptic both measured it)
* but: WHAT DOES THE USER LOSE, and is the proposed remedy (adopt the regress
* covariance) safe?
*
* Three arms on the SAME draw, all customvce(robust) (= the regress engine, so
* an SE always exists), 12 draws:
*   A1  lattice(all 3)   - dose-quantile lattice  (prober/skeptic design)
*   B1  lattice(reg3)    - substantive lattice INDEPENDENT of dose, wide
*                          within-bin dose variation  (the most natural applied
*                          use of lattice(): "let the dose effect vary by region")
*   N1  no lattice       - baseline calibration arm
* For A1 and B1 we ALSO run the wooldid DEFAULT (reghdfe) path and record
* whether the SE came back missing, so the zero-V incidence is measured on both.
* We track the att row and the ame row SEPARATELY: they are not the same casualty.
* Truth is computed PER DRAW from the DGP.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do L1_impact.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

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

tempname MC
postfile `MC' double(attA seA ameA seameA defmissA) ///
               double(attB seB ameB seameB defmissB) ///
               double(attN seN ameN seameN) ///
               double(truatt truame) using L1mc , replace

forval rep = 1/12 {
  mkdata
  tempfile F
  qui save `F'
  * ---- per-draw DGP truth over treated-post rows ----
  qui sum dose if post == 1
  local truame = 0.4
  qui gen double tefx = 0.4*dose + 0.5 if post == 1
  qui sum tefx if post == 1
  local truatt = r(mean)
  drop tefx
  qui save `F' , replace

  foreach arm in A B N {
    if "`arm'" == "A" local lat "lattice(all 3)"
    if "`arm'" == "B" local lat "lattice(reg3)"
    if "`arm'" == "N" local lat ""
    * ---- regress engine (customvce): always yields an SE ----
    use `F', clear
    local b1 = .
    local s1 = .
    local b2 = .
    local s2 = .
    cap qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) `lat' customvce(robust)
    if _rc == 0 {
      mat M = e(wdidmainresults)
      local b1 = M[1,1]
      local s1 = M[1,2]
      local b2 = M[2,1]
      local s2 = M[2,2]
    }
    local att`arm' = `b1'
    local se`arm'  = `s1'
    local ame`arm' = `b2'
    local seame`arm' = `s2'
    * ---- default (reghdfe) engine: did the SE vanish? ----
    if "`arm'" != "N" {
      use `F', clear
      local dm = .
      cap qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) `lat' robust
      if _rc == 0 {
        mat M2 = e(wdidmainresults)
        local dm = (M2[1,2] >= .)
        local pb1 = M2[1,1]
        local pb2 = M2[2,1]
        dis as text "  rep`rep' arm`arm' default-path att b=" %14.8f `pb1' "  cvce att b=" %14.8f `b1' "  attSEmissing=`dm'"
      }
      local defmiss`arm' = `dm'
    }
  }
  post `MC' (`attA') (`seA') (`ameA') (`seameA') (`defmissA') ///
            (`attB') (`seB') (`ameB') (`seameB') (`defmissB') ///
            (`attN') (`seN') (`ameN') (`seameN') ///
            (`truatt') (`truame')
  dis as result "REP`rep' truatt=" %9.5f `truatt' " | A att=" %9.5f `attA' " se=" %8.5f `seA' " ame=" %10.5f `ameA' " se=" %8.5f `seameA' " defmiss=`defmissA'"
  dis as result "         " _col(22) " | B att=" %9.5f `attB' " se=" %8.5f `seB' " ame=" %10.5f `ameB' " se=" %8.5f `seameB' " defmiss=`defmissB'"
  dis as result "         " _col(22) " | N att=" %9.5f `attN' " se=" %8.5f `seN' " ame=" %10.5f `ameN' " se=" %8.5f `seameN'
}
postclose `MC'

use L1mc , clear
dis as result "================ SUMMARY: SE vs realized sampling SD ================"
foreach v in attA attB attN ameA ameB ameN {
  qui sum `v'
  local sd_`v' = r(sd)
  local mn_`v' = r(mean)
}
foreach v in seA seB seN seameA seameB seameN {
  qui sum `v'
  local ms_`v' = r(mean)
}
dis as result "att  A(all 3)   : mean=" %10.5f `mn_attA' " MCSD=" %10.5f `sd_attA' " meanSE=" %10.5f `ms_seA' "  SE/SD=" %8.3f `ms_seA'/`sd_attA'
dis as result "att  B(reg3)    : mean=" %10.5f `mn_attB' " MCSD=" %10.5f `sd_attB' " meanSE=" %10.5f `ms_seB' "  SE/SD=" %8.3f `ms_seB'/`sd_attB'
dis as result "att  N(nolat)   : mean=" %10.5f `mn_attN' " MCSD=" %10.5f `sd_attN' " meanSE=" %10.5f `ms_seN' "  SE/SD=" %8.3f `ms_seN'/`sd_attN'
qui sum truatt
dis as result "   (true att per-draw mean = " %10.5f r(mean) ", sd " %10.5f r(sd) ")"
dis as result "ame  A(all 3)   : mean=" %10.5f `mn_ameA' " MCSD=" %10.5f `sd_ameA' " meanSE=" %10.5f `ms_seameA' "  SE/SD=" %8.3f `ms_seameA'/`sd_ameA'
dis as result "ame  B(reg3)    : mean=" %10.5f `mn_ameB' " MCSD=" %10.5f `sd_ameB' " meanSE=" %10.5f `ms_seameB' "  SE/SD=" %8.3f `ms_seameB'/`sd_ameB'
dis as result "ame  N(nolat)   : mean=" %10.5f `mn_ameN' " MCSD=" %10.5f `sd_ameN' " meanSE=" %10.5f `ms_seameN' "  SE/SD=" %8.3f `ms_seameN'/`sd_ameN'
dis as result "   (true ame = 0.400)"
qui count if defmissA == 1
dis as result "zero-V incidence: A(all 3) " r(N) "/12"
qui count if defmissB == 1
dis as result "zero-V incidence: B(reg3)  " r(N) "/12"
* coverage of the nominal 95% CI on the recovered SEs
foreach a in A B N {
  qui count if abs(att`a' - truatt) <= 1.96*se`a'
  local ca = r(N)
  qui count if abs(ame`a' - 0.4) <= 1.96*seame`a'
  local cm = r(N)
  dis as result "cover95 arm `a' : att " `ca' "/12   ame " `cm' "/12"
}
dis "DONE_L1"
