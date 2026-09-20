clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 99001
set sortseed 99001

cap prog drop truths
prog def truths, rclass
  syntax varname , IFC(string) [WVAR(string)]
  tempvar q ww wperit wperi tperi tag
  if "`wvar'" == "" qui gen double `ww' = 1
  else qui gen double `ww' = `wvar'
  qui sum `varlist' [aw=`ww'] if `ifc'
  return scalar att = r(mean)
  qui gegen double `wperit' = total(`ww') if `ifc', by(i t)
  qui gen double `q' = `ww'/`wperit' if `ifc'
  qui sum `varlist' [aw=`q'] if `ifc'
  return scalar att_it = r(mean)
  qui gegen double `wperi' = total(`ww') if `ifc', by(i)
  qui replace `q' = `ww'/`wperi' if `ifc'
  qui sum `varlist' [aw=`q'] if `ifc'
  return scalar att_i = r(mean)
  qui gegen byte `tag' = tag(i t) if `ifc'
  qui gegen double `tperi' = total(`tag') if `ifc', by(i)
  qui replace `q' = `ww'/(`wperit'*`tperi') if `ifc'
  qui sum `varlist' [aw=`q'] if `ifc'
  return scalar att_itime = r(mean)
end

* UNBALANCED, messy-but-plausible panel: 6 groups x 10 units x 10 periods,
* 22% of unit-periods missing at random, exact-fit outcome.
cap prog drop mkunb
prog def mkunb
  clear
  set obs 600
  gen long unit = ceil(_n/10)
  gen int  t    = _n - (unit-1)*10
  gen byte i    = ceil(unit/10)
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  qui replace ttre = 8 if i == 3
  gen double relyr = t - ttre
  gen double dose = 1 + 2*runiform()
  gen double w    = 0.4 + 1.6*runiform()
  gen double b1 = 0
  gen double k0 = 0
  qui replace k0 = 0.3 + 0.11*i + 0.07*relyr if relyr>=0 & relyr<.
  qui replace b1 = 0.8 + 0.25*i - 0.13*relyr if relyr>=0 & relyr<.
  gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose + k0 + b1*dose
  qui drop if runiform() < 0.22
end

foreach wt in "none" "aw" {
 foreach est in att att_it att_i att_itime {
   set seed 99001
   qui mkunb
   local wgt
   local wv
   if "`wt'" == "aw" {
     local wgt [aw=w]
     local wv w
   }
   truths b1 , ifc(ttre<. & t>=ttre) wvar(`wv')
   local T = r(`est')
   cap qui wooldid y i t ttre `wgt' , `est' contreat(dose) contreatcontrols(basic) robust regtol(14)
   if _rc == 0 {
     mat M = e(wdidmainresults)
     dis "UNB|`wt'|`est'|ame=" %16.10f M[2,1] "|truth=" %16.10f `T' "|reld=" %10.3e abs(M[2,1]-`T')/abs(`T') "|r2=" %14.12f e(r2) "|N=" e(N)
   }
   else dis "UNB|`wt'|`est'|RC=" _rc
 }
}
dis "DONE"
