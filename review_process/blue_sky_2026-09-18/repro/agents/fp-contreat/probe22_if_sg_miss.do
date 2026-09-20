clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242

cap prog drop mk22
prog def mk22
  clear
  set obs 512                       /* 8 groups x 8 units x 8 periods */
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/8)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 5 if i == 2
  qui replace ttre = 6 if i == 3
  qui replace ttre = 7 if i == 4
  gen byte sg = mod(i,2)            /* constant within i, hence within (i,t) */
  gen double dose = 1 + 2*runiform()
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre<. & t>=ttre
  qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre<. & t>=ttre
  gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose + k0 + b1*dose
end

* ---- (v) contreat() without contreatcontrols() should be refused ---------
qui mk22
cap noisily wooldid y i t ttre , att contreat(dose) robust regtol(14)
dis "NOCC_RC=" _rc

* ---- (ii) if restriction -------------------------------------------------
foreach ifc in "t <= 6" "i != 2" "dose > 1.5" {
  qui mk22
  qui sum b1 if (`ifc') & ttre<. & t>=ttre
  local tru = r(mean)
  local nn  = r(N)
  cap qui wooldid y i t ttre if `ifc' , att contreat(dose) contreatcontrols(basic) robust regtol(14)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "IF|`ifc'|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|r2=" %14.12f e(r2) "|N=" e(N)
  }
  else dis "IF|`ifc'|RC=" _rc
}

* ---- (iii) subgroups -----------------------------------------------------
qui mk22
forval s = 0/1 {
  qui sum b1 if sg==`s' & ttre<. & t>=ttre
  local T`s' = r(mean)
}
cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) subgroups(sg) robust regtol(14)
dis "SG_RC=" _rc
if _rc==0 {
  dis "SG r2=" e(r2)
  forval s = 0/1 {
    cap mat S = e(wdidsgresults_`s')
    if _rc == 0 {
      local rn : rownames S
      local re : roweq S
      dis "SGROWNAMES`s'=`rn'"
      dis "SGROWEQ`s'=`re'"
      local nr = rowsof(S)
      forval r = 1/`nr' {
        local nm : word `r' of `rn'
        local eq : word `r' of `re'
        if "`eq'" == "main" & substr("`nm'",1,3) == "ame" {
          dis "SG|`s'|ame=" %16.10f S[`r',1] "|truth=" %16.10f `T`s'' "|reld=" %10.3e abs(S[`r',1]-`T`s'')/abs(`T`s'')
        }
      }
    }
    else dis "SG|`s'|no matrix e(wdidsgresults_`s')"
  }
}

* ---- (iv) missing dose on 10% of rows ------------------------------------
qui mk22
qui replace dose = . if runiform() < 0.10
qui replace y = . if dose >= .
qui replace y = 2 + 0.6*i + 0.35*t + 0.9*dose + k0 + b1*dose if dose < .
qui sum b1 if ttre<. & t>=ttre & dose<.
local tru = r(mean)
local nn = r(N)
cap noisily qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) robust regtol(14)
if _rc==0 {
  mat M = e(wdidmainresults)
  dis "MISS|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|r2=" %14.12f e(r2) "|N=" e(N) "|Ndrop=" e(N_dropped)
}
else dis "MISS|RC=" _rc
dis "DONE22"
