*======================================================================
* r04_s7geom -- SELF-CONTAINED.  THE experiment MCBP-S7COVER says would settle
* it: bootstrap(standard) coverage of att_i at cluster()==i, holding G = 10
* states fixed and varying ONLY the number of treated states per onset cohort.
*   A = 3 cohorts x 1 state  (3 treated, 7 never)   <- the §7 geometry's feature
*   B = 3 cohorts x 2 states (6 treated, 4 never)   <- the seeded finding's
* Sharp null (truth = 0), reps(99), 45 draws, identical data seeds across arms.
* Analytic arms reported beside them for reference.
* Run:  stata-mp -e do r04_s7geom.do           (~35 min)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpan
prog def mkpan
  syntax , G(integer) T(integer) SEEDX(integer) [KPER(integer 2) RHO(real 0.4) NTRTC(integer 3) NSUB(integer 1)]
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
  if `nsub' > 1 {
    expand `nsub'
    bysort st t: gen int sub = _n
    gen double v = rnormal()*0.8
  }
  else gen double v = 0
  gen double y = ai + ft + e + v
  keep st t ttre y
  order st t ttre y
end

* pull one main-table row out of e(wdidmainresults): b se stat p lo hi
cap prog drop getrow
prog def getrow, rclass
  syntax , ROWN(string)
  return scalar b  = .
  return scalar se = .
  return scalar p  = .
  return scalar lo = .
  return scalar hi = .
  cap mat TT = e(wdidmainresults)
  cap confirm matrix TT
  if _rc exit
  local rn : rownames TT
  local re : roweq TT
  local nr : word count `rn'
  forvalues z = 1/`nr' {
    local nz : word `z' of `rn'
    local ez : word `z' of `re'
    if "`nz'"=="`rown'" & "`ez'"=="main" {
      return scalar b  = TT[`z',1]
      return scalar se = TT[`z',2]
      return scalar p  = TT[`z',4]
      return scalar lo = TT[`z',5]
      return scalar hi = TT[`z',6]
      continue, break
    }
  }
end

tempname P
postfile `P' int draw str6 arm double b double se double p double lo double hi ///
   double bi double sei double pi double loi double hii double ncomp int rc using r04_out.dta, replace

local ND 45
forvalues d = 1/`ND' {
  local sx = 720000 + `d'
  local rs = 950000 + 17*`d'
  * ---- geometry A: one treated state per cohort
  mkpan , g(10) t(10) seedx(`sx') kper(1) ntrtc(3) nsub(1)
  tempfile da
  qui save "`da'", replace
  * ---- geometry B: two treated states per cohort, same data seed
  mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(1)
  tempfile db
  qui save "`db'", replace
  foreach a in A_BS B_BS A_AN B_AN {
    local f = cond(substr("`a'",1,1)=="A", "`da'", "`db'")
    qui use "`f'", clear
    ereturn clear
    if substr("`a'",3,2)=="BS" cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(99) rseed(`rs')
    else                        cap noi qui wooldid y st t ttre, cluster(st) att att_i
    local rc = _rc
    local nc = .
    if `rc'==0 cap local nc = e(boot_completed)
    qui getrow , rown(att)
    local b = r(b)
    local se = r(se)
    local p = r(p)
    local lo = r(lo)
    local hi = r(hi)
    qui getrow , rown(att_i)
    post `P' (`d') ("`a'") (`b') (`se') (`p') (`lo') (`hi') (r(b)) (r(se)) (r(p)) (r(lo)) (r(hi)) (`nc') (`rc')
  }
  if mod(`d',5)==0 dis "... r04 draw `d'"
}
postclose `P'

use r04_out.dta, clear
gen byte cov  = (lo  <= 0 & hi  >= 0) if lo<.
gen byte covi = (loi <= 0 & hii >= 0) if loi<.
gen double widi = hii - loi
dis "=== R04  att_i coverage at cluster(st)==i, sharp null, G=10 ==="
foreach a in A_BS B_BS A_AN B_AN {
  qui sum covi if arm=="`a'"
  local cvi = r(mean)
  local ni = r(N)
  qui sum cov if arm=="`a'"
  local cv = r(mean)
  qui sum bi if arm=="`a'"
  local sdi = r(sd)
  qui sum sei if arm=="`a'"
  local msi = r(mean)
  qui sum widi if arm=="`a'"
  dis "R04 `a': n=`ni' cov(att_i)=" %6.4f `cvi' "  cov(att)=" %6.4f `cv' ///
      "  sd(att_i)=" %7.4f `sdi' "  mean se=" %7.4f `msi' "  se/sd=" %6.3f `=`msi'/`sdi'' "  mean width=" %7.4f r(mean)
}
count if rc != 0
dis "R04 NONZERO_RC=" r(N)
tab arm rc
dis "R04_DONE"
