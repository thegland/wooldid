* R7 -- control arms for the R6 measurement.
* A: i=unit, cluster(unit)   (cluster AT the level of i -- the documented bad case)
* B: i=coarse cohort id + fe(unit), cluster(unit)  (the help's INTENDED practice: cluster below i)
* C: i=unit, robust
* Same DGP in all three.  If B lands at SE/SD ~1 and coverage ~.95, the harness is sound
* and arm A's number is wooldid's, not mine.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
which wooldid

local R = 60
mat A = J(`R',3,.)
mat B = J(`R',3,.)
mat C = J(`R',3,.)
forvalues r = 1/`R' {
  qui {
    clear
    set obs 300
    gen long unit = _n
    gen byte tg = (unit <= 20)
    gen long cohid = cond(tg==1, mod(unit,3)+1, 9)
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
    mat A[`r',1]=bb[1,1]
    mat A[`r',2]=sqrt(vv[1,1])
    mat A[`r',3]= abs(bb[1,1]-1.5) <= invttail(e(df_r),0.025)*sqrt(vv[1,1])
  }
  cap qui wooldid y cohid t ttre, att cluster(unit) fe(unit)
  if _rc==0 {
    mat bb=e(b)
    mat vv=e(V)
    mat B[`r',1]=bb[1,1]
    mat B[`r',2]=sqrt(vv[1,1])
    mat B[`r',3]= abs(bb[1,1]-1.5) <= invttail(e(df_r),0.025)*sqrt(vv[1,1])
  }
  cap qui wooldid y unit t ttre, att robust
  if _rc==0 {
    mat bb=e(b)
    mat vv=e(V)
    mat C[`r',1]=bb[1,1]
    mat C[`r',2]=sqrt(vv[1,1])
    mat C[`r',3]= abs(bb[1,1]-1.5) <= invttail(e(df_r),0.025)*sqrt(vv[1,1])
  }
}
foreach arm in A B C {
  clear
  svmat `arm'
  qui sum `arm'1
  local sd = r(sd)
  local n  = r(N)
  local mb = r(mean)
  qui sum `arm'2
  local ms = r(mean)
  qui sum `arm'3
  local cv = r(mean)
  dis "ARM`arm'_NREP=" `n' "  MEANB=" %9.6f `mb' "  SD_B=" %9.6f `sd' "  MEAN_SE=" %9.6f `ms' "  SE_OVER_SD=" %9.6f `ms'/`sd' "  COV95=" %6.4f `cv'
}
dis "R7_DONE"
