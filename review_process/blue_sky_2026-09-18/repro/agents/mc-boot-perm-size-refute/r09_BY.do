*======================================================================
* r09_BY -- SELF-CONTAINED.  Independent size test of bootstrap(BY) at
* G=10 with cluster()==i under a TRUE SHARP NULL.  130 fresh draws, reps(79)
* (the seeded findings' own replication count).  MCBP-BAYES claims .163-.204
* rejection for bayes; MCBP-BOOTSTD claims .122 for standard.
* Run:  stata-mp -e do r09_BY.do          (~40 min)
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
postfile `P' int draw double b double se double p double lo double hi double bi double sei double pi double loi double hii double ncomp int rc using r09_BY_out.dta, replace
forvalues d = 1/130 {
  local sx = 1310000 + `d'
  local rs = 1440000 + 19*`d'
  mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(1)
  ereturn clear
  if "BY"=="BY" cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(bayes) reps(79) rseed(`rs')
  if "BY"=="BS" cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs')
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
  post `P' (`d') (`b') (`se') (`p') (`lo') (`hi') (r(b)) (r(se)) (r(p)) (r(lo)) (r(hi)) (`nc') (`rc')
  if mod(`d',20)==0 dis "... r09_BY draw `d'"
}
postclose `P'
use r09_BY_out.dta, clear
gen byte rej = (p <= .05+1e-9) if p<.
gen byte cov = (lo <= 0 & hi >= 0) if lo<.
gen byte reji = (pi <= .05+1e-9) if pi<.
gen byte covi = (loi <= 0 & hii >= 0) if loi<.
qui sum b
local sd = r(sd)
local n = r(N)
qui sum se
local ms = r(mean)
qui sum rej
local rj = r(mean)
local nr = r(sum)
qui sum cov
dis "R09_BY att: n=`n' sd(b)=" %7.4f `sd' " mse=" %7.4f `ms' " se/sd=" %6.3f `=`ms'/`sd'' " rej@5%=" %6.3f `rj' " (`nr' of `n') cov=" %6.3f r(mean)
qui sum reji
if r(N) > 0 {
  local rji = r(mean)
  local nri = r(sum)
  qui sum bi
  local sdi = r(sd)
  qui sum sei
  local msi = r(mean)
  qui sum covi
  dis "R09_BY att_i: sd(b)=" %7.4f `sdi' " mse=" %7.4f `msi' " se/sd=" %6.3f `=`msi'/`sdi'' " rej@5%=" %6.3f `rji' " (`nri') cov=" %6.3f r(mean)
}
foreach x in .05 .10 .20 .30 {
  qui count if p <= `x'+1e-9 & p<.
  local c1 = r(N)
  qui count if p<.
  dis "R09_BY P(p<=`x') = " %6.3f `=`c1'/r(N)'
}
qui sum ncomp
dis "R09_BY completions min=" r(min)
count if rc != 0
dis "R09_BY NONZERO_RC=" r(N)
dis "R09_BY_DONE"
