capture program drop mkfix
program define mkfix
  args sd withlatu
  clear
  set obs 12
  gen int st = _n
  gen int nc = 5 + mod(_n*7,8)
  gen double ttre = .
  qui replace ttre =  8 if inlist(st,1,2,3)
  qui replace ttre = 12 if inlist(st,4,5,6)
  qui replace ttre = 15 if inlist(st,7,8)
  expand nc
  gen long unit = _n
  expand 20
  bysort unit: gen int year = _n
  bysort unit (year): gen double ubase = exp(rnormal(0,0.6)) if _n==1
  bysort unit (year): replace ubase = ubase[1]
  gen double dose_it = ubase*(1 + 0.10*rnormal())
  bysort unit (year): gen double ufe = rnormal(0,1) if _n==1
  bysort unit (year): replace ufe = ufe[1]
  gen double y_dose = 3 + ufe + 0.2*year + 0.5*dose_it + rnormal(0,1)
  qui replace y_dose = y_dose + 0.4*dose_it + 0.3 if ttre<. & year>=ttre
  if "`withlatu'"=="1" qui gquantiles latu = dose_it , xtile nquantiles(3)
end
