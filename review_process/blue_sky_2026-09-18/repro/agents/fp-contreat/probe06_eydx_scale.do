clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* self-contained generator with a tunable fixed-effect scale `s'
cap prog drop mkdat
prog def mkdat
  args s
  clear
  set obs 256
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/8)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  gen double kappatrue = 0
  gen double betatrue  = 0
  qui replace kappatrue = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre < . & t >= ttre
  qui replace betatrue  = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre < . & t >= ttre
  gen double y = 6 + `s'*(0.6*i + 0.35*t + 0.9*dose) + kappatrue + betatrue*dose
end

foreach s in 0 0.1 1 3 {
  foreach et in eydx eyex {
    foreach path in "robust" "customvce(robust)" {
      set seed 771133
      qui mkdat `s'
      tempvar q
      if "`et'"=="eydx" qui gen double `q' = betatrue/y
      if "`et'"=="eyex" qui gen double `q' = betatrue*dose/y
      qui sum `q' if ttre<. & t>=ttre
      local tru = r(mean)
      cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') `path' regtol(14)
      if _rc==0 {
        mat M = e(wdidmainresults)
        dis "SC|s=`s'|`et'|`path'|ame=" %18.12f M[2,1] "|truth=" %18.12f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru')
      }
      else dis "SC|s=`s'|`et'|`path'|RC=" _rc
    }
  }
}

* --- same mechanism on the BINARY side: semielasticity (eydx of targetflag) ---
foreach s in 0 1 3 {
  foreach path in "robust" "customvce(robust)" {
    set seed 771133
    qui mkdat `s'
    * binary truth: mean over treated post rows of (kappa + beta*dose)/y  -- the
    * effect of switching treatment off is exactly -(kappa+beta*dose) under exact fit,
    * and margins' eydx reports d ln yhat/dD = (kappa+beta*dose)/yhat
    tempvar q
    qui gen double `q' = (kappatrue + betatrue*dose)/y
    qui sum `q' if ttre<. & t>=ttre
    local tru = r(mean)
    cap qui wooldid y i t ttre , att semielasticity `path' regtol(14)
    if _rc==0 {
      mat M = e(wdidmainresults)
      dis "SEMI|s=`s'|`path'|att=" %18.12f M[1,1] "|truth=" %18.12f `tru' "|reld=" %10.3e abs(M[1,1]-`tru')/abs(`tru')
    }
    else dis "SEMI|s=`s'|`path'|RC=" _rc
  }
}
dis "DONE"
