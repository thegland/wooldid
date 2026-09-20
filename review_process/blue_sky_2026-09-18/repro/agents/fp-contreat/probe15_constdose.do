clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150

cap prog drop mk
prog def mk
  args mode
  clear
  set obs 288                    /* 6 groups x 6 units x 8 periods */
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/6)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  if "`mode'" == "global"  qui replace dose = 2.5
  if "`mode'" == "bycell"  {
    tempvar cd
    qui gegen double `cd' = mean(dose), by(i t)
    qui replace dose = `cd'
  }
  if "`mode'" == "byunit"  {
    tempvar ud
    qui gegen double `ud' = mean(dose), by(unit)
    qui replace dose = `ud'
  }
  gen double eff = 0
  qui replace eff = 1.0 + 0.6*dose if ttre<. & t>=ttre
  gen double y = 5 + 0.7*i + 0.4*t + 0.9*dose + eff + rnormal(0,1)
end

foreach mode in vary global bycell byunit {
  dis "########## DOSE MODE = `mode'"
  qui mk `mode'
  qui sum dose
  dis "dose: min=" r(min) " max=" r(max) " sd=" r(sd)
  cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) robust regtol(14)
  dis "CD|`mode'|RC=" _rc
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis "CD|`mode'|att=" %16.8f M[1,1] "|se=" %14.8g M[1,2] "|ame=" %16.8f M[2,1] "|se=" %14.8g M[2,2] "|r2=" %8.6f e(r2)
  }
}
dis "DONE"
