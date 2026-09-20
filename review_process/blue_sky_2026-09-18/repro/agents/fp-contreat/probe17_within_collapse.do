clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 27182
set sortseed 27182

* ================= (d) contreatwithin: every unit eventually treated =======
cap prog drop mkw
prog def mkw
  clear
  set obs 320
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/10)
  gen double ttre = 3 + i           /* cohorts treated at 4,5,6,7 - NO never-treated */
  gen double dose = 1 + 2*runiform()
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.07*(t-ttre) if t>=ttre
  qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if t>=ttre
  gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose + k0 + b1*dose
end
foreach path in "robust" "customvce(robust)" {
  qui mkw
  qui sum b1 if t>=ttre
  local tru = r(mean)
  cap noisily qui wooldid y i t ttre , att contreat(dose) contreatwithin contreatcontrols(basic) `path' regtol(14)
  if _rc==0 {
    mat M = e(wdidmainresults)
    dis "CW|`path'|ame=" %16.10f M[2,1] "|truth=" %16.10f `tru' "|reld=" %10.3e abs(M[2,1]-`tru')/abs(`tru') "|r2=" %14.12f e(r2)
  }
  else dis "CW|`path'|RC=" _rc
}

* ================= (e) cell collapse vs nofastcollapse =====================
* noisy data, all four estimands, weighted, dydx and dyex: must be bit-identical
cap prog drop mkn
prog def mkn
  clear
  set obs 384
  gen long unit = ceil(_n/8)
  gen int  t    = _n - (unit-1)*8
  gen byte i    = ceil(unit/8)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 5 if i == 2
  qui replace ttre = 6 if i == 3
  gen double dose = 1 + 2*runiform()
  gen double w = 0.3 + 2*runiform()
  gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose + rnormal(0,1)
  qui replace y = y + 1.0 + 0.7*dose if ttre<. & t>=ttre
end
foreach est in att att_it att_i att_itime {
 foreach et in dydx dyex {
   local v1 .
   local v2 .
   foreach fc in "" "nofastcollapse" {
     set seed 27182
     qui mkn
     cap qui wooldid y i t ttre [aw=w] , `est' contreat(dose) contreatcontrols(basic) ///
        contreatelasticitytype(`et') robust esprelength(1) espostlength(2) regtol(14) `fc'
     if _rc==0 {
       mat M = e(wdidmainresults)
       if "`fc'" == "" local v1 = M[3,1]
       else local v2 = M[3,1]
     }
   }
   dis "FC|`est'|`et'|collapse=" %21x `v1' "|full=" %21x `v2' "|equal=" (`v1'==`v2')
 }
}
dis "DONE"
