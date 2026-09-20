*======================================================================
* s00_smoke -- assert build, see which cautions fire at G=10/G=9,
*              and time the analytic + bootstrap arms.
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 771001
set sortseed 771001

* ---- my own generator: G states, T periods, NC cohorts of KP states each ----
cap prog drop mkpan
prog def mkpan
  syntax , G(integer) T(integer) SEEDX(integer) [KPER(integer 2) RHO(real 0.4) NTRTC(integer 3)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `g'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  gen double ttre = .
  forvalues c = 1/`ntrtc' {
    local lo = (`c'-1)*`kper' + 1
    local hi = `c'*`kper'
    local on = 2003 + 2*`c'
    replace ttre = `on' if rk >= `lo' & rk <= `hi'
  }
  gen double ai = rnormal()*1.2
  drop u rk
  expand `t'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.6*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal() if _n==1
  by st: replace e = `rho'*e[_n-1] + rnormal()*sqrt(1-`rho'^2) if _n>1
  gen double y = ai + ft + e
  keep st t ttre y
  order st t ttre y
end

* ===== A. which cautions fire?  G=10 (6 treated) vs G=9 (6 treated) vs G=10 (4 treated)
dis "===== CAUTIONCHECK G10T6 ====="
mkpan , g(10) t(10) seedx(9001) kper(2) ntrtc(3)
qui tab ttre, missing
cap noi wooldid y st t ttre, cluster(st) att bootstrap(standard) reps(19) rseed(11)
dis "G10T6_RC=" _rc

dis "===== CAUTIONCHECK G9T6 ====="
mkpan , g(9) t(10) seedx(9002) kper(2) ntrtc(3)
cap noi wooldid y st t ttre, cluster(st) att bootstrap(standard) reps(19) rseed(11)
dis "G9T6_RC=" _rc

dis "===== CAUTIONCHECK G12T6 ====="
mkpan , g(12) t(10) seedx(9003) kper(2) ntrtc(3)
cap noi wooldid y st t ttre, cluster(st) att bootstrap(standard) reps(19) rseed(11)
dis "G12T6_RC=" _rc

* ===== B. timing: 5 draws x 4 analytic specs
timer clear
timer on 1
forvalues d = 1/5 {
  mkpan , g(10) t(10) seedx(`=770000+`d'') kper(2) ntrtc(3)
  tempfile dd
  qui save "`dd'", replace
  qui use "`dd'", clear
  cap noi qui wooldid y st t ttre, cluster(st) att
  qui use "`dd'", clear
  cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3)
  qui use "`dd'", clear
  cap noi qui wooldid y st t ttre, cluster(st) att esfixedbaseperiod
  qui use "`dd'", clear
  cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3) esfixedbaseperiod
}
timer off 1

* ===== C. timing: 2 draws bootstrap(standard) reps 199
timer on 2
forvalues d = 1/2 {
  mkpan , g(10) t(10) seedx(`=778000+`d'') kper(2) ntrtc(3)
  cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(199) rseed(`=1000+`d'')
  dis "BSTD_RC=" _rc
}
timer off 2
timer on 3
forvalues d = 1/2 {
  mkpan , g(10) t(10) seedx(`=778000+`d'') kper(2) ntrtc(3)
  cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(bayes) reps(199) rseed(`=1000+`d'')
  dis "BBAY_RC=" _rc
}
timer off 3
timer on 4
mkpan , g(10) t(10) seedx(778001) kper(2) ntrtc(3)
cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) boottstat reps(199) rseed(1001)
dis "BTST_RC=" _rc
timer off 4
timer list
dis "SMOKE_DONE"
