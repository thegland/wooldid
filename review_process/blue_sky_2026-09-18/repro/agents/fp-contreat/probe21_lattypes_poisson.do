clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606060
set sortseed 606060
which wooldid
which reghdfe
which ppmlhdfe
which gtools

* ============ (5) lattice TYPE invariance, exact fit, no bin structure =====
cap prog drop mklt
prog def mklt
  clear
  set obs 384                        /* 4 groups x 12 units x 8 periods */
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/12)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  gen double w    = 0.5 + runiform()
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre<. & t>=ttre
  qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre<. & t>=ttre
  gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose + k0 + b1*dose
end

foreach lt in all byTime byTreated byTtre byCohort byTreatedByTime byTtreByTime byCohortByTime baseTreated baseTreatedRef baseTreatedPost baseTreatedPre {
  foreach extra in "" "latticeignoreweights" {
    set seed 606060
    qui mklt
    qui sum b1 if ttre<. & t>=ttre
    local tru = r(mean)
    local wgt
    if "`extra'" != "" local wgt [aw=w]
    cap qui wooldid y i t ttre `wgt' , att contreat(dose) contreatcontrols(basic) lattice(`lt' 3) `extra' robust regtol(14)
    local rc = _rc
    if `rc' == 0 {
      * with weights the truth is the weighted mean
      if "`extra'" != "" {
        qui sum b1 [aw=w] if ttre<. & t>=ttre
        local tru = r(mean)
      }
      mat M = e(wdidmainresults)
      dis "LT|`lt'|`extra'|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|se=" %12.6g M[2,2] "|r2=" %14.12f e(r2)
    }
    else dis "LT|`lt'|`extra'|RC=`rc'"
  }
}

* ============ (6) poisson contreat, exact-fit log-linear ==================
cap prog drop mkpois
prog def mkpois
  clear
  set obs 384
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/12)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.10 + 0.03*i + 0.02*(t-ttre) if ttre<. & t>=ttre
  qui replace b1 = 0.20 + 0.05*i - 0.03*(t-ttre) if ttre<. & t>=ttre
  gen double xb = 1 + 0.15*i + 0.05*t + 0.10*dose + k0 + b1*dose
  gen double y = exp(xb)
end
foreach et in dydx dyex eyex eydx {
  set seed 606060
  qui mkpois
  qui sum b1 if ttre<. & t>=ttre
  local T_dydx = r(mean)
  local T_eydx = r(mean)
  tempvar q
  qui gen double `q' = b1*dose
  qui sum `q' if ttre<. & t>=ttre
  local T_dyex = r(mean)
  local T_eyex = r(mean)
  cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) poisson ///
       contreatelasticitytype(`et') robust regtol(14)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "POIS|`et'|ame=" %16.10f M[2,1] "|truth=" %16.10f `T_`et'' "|reld=" %10.3e abs(M[2,1]-`T_`et'')/abs(`T_`et'')
  }
  else dis "POIS|`et'|RC=" _rc
}
dis "DONE21"
