* R8 -- (D) calibration arm for SKEP-01: SAME spec as arm A (i=unit, cluster(unit)) but a
*           SMALL control pool (20 treated of 60 units).  Mechanism prediction: SE/SD is
*           driven by the control-to-treated ratio, so this arm should land near the
*           documented regime while the 20-of-300 arm does not.
*       (E) does the si8 preflight refuse the very remedy the help recommends for i=unit?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
which wooldid

local R = 60
mat D = J(`R',3,.)
forvalues r = 1/`R' {
  qui {
    clear
    set obs 60
    gen long unit = _n
    gen byte tg = (unit <= 20)
    expand 8
    bysort unit: gen int t = 1999 + _n
    gen double ttre = .
    replace ttre = 2002 + mod(unit,3) if tg == 1
    gen double ui = rnormal()
    bysort unit (t): replace ui = ui[1]
    gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
  }
  cap qui wooldid y unit t ttre, att cluster(unit)
  if _rc==0 {
    mat bb=e(b)
    mat vv=e(V)
    mat D[`r',1]=bb[1,1]
    mat D[`r',2]=sqrt(vv[1,1])
    mat D[`r',3]= abs(bb[1,1]-1.5) <= invttail(e(df_r),0.025)*sqrt(vv[1,1])
  }
}
clear
svmat D
qui sum D1
local sd=r(sd)
local n=r(N)
qui sum D2
local ms=r(mean)
qui sum D3
local cv=r(mean)
dis "ARMD_20of60_NREP=" `n' " SD_B=" %9.6f `sd' " MEAN_SE=" %9.6f `ms' " SE_OVER_SD=" %9.6f `ms'/`sd' " COV95=" %6.4f `cv'

* ---- (E) does si8 refuse the help's own remedy for i=unit? ----
clear
set obs 4800
gen long unit = _n
gen byte tg = (unit <= 20)
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 + mod(unit,3) if tg == 1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
dis "--- E1 bootstrap(standard) reps(2) cluster(unit), default safety, maxvar " c(maxvar)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(2) rseed(1)
dis "E1_RC=" _rc
dis "R8_DONE"
