clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 40404
set sortseed 40404
clear
set obs 60
gen unit = _n
gen cg = ceil(_n/15)
gen ttre = .
replace ttre = 5 if cg==1
replace ttre = 7 if cg==2
replace ttre = 9 if cg==3
gen double ui = rnormal(0,0.80)
expand 12
bysort unit: gen t = _n
gen double te = 0.05*t
gen double eps = rnormal(0,0.60)
gen byte post = !missing(ttre) & t>=ttre
gen double dose = 0
replace dose = runiform(0.2,2.0) if post
gen double ycont = 2 + ui + te + post*(0.45*dose) + eps
foreach cc in untreatedZero {
  dis "=== ame cluster(unit), contreatcontrols(`cc') ==="
  capture noisily wooldid ycont cg t ttre, att cluster(unit) impvalues contreat(dose) contreatcontrols(`cc')
  dis "U_RC=" _rc
  capture matrix CU = e(wdidmainresults)
  capture noisily matrix list CU
  dis "=== ame cluster(cg) ==="
  capture noisily wooldid ycont cg t ttre, att cluster(cg) impvalues contreat(dose) contreatcontrols(`cc')
  dis "G_RC=" _rc
  capture matrix CG = e(wdidmainresults)
  capture noisily matrix list CG
  capture {
    local nc = rowsof(CU)
    forvalues j = 1/`nc' {
      local rn : word `j' of `: rownames CU'
      dis "CROW `j' `rn' : se_u=" %10.7f CU[`j',2] " se_g=" %10.7f CG[`j',2] " ratio=" %8.5f CU[`j',2]/CG[`j',2]
    }
  }
}
dis "PRED= " sqrt(3*15*2/(44*45))
dis "DONE_CONTREAT"
