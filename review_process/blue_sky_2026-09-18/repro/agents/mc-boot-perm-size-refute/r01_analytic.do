*======================================================================
* r01_analytic -- SELF-CONTAINED refutation probe.
*  Q1 (MCBP-ESFBP): at cluster()==i on a G=10 sharp-null staggered panel,
*      is the DEFAULT att analytic SE understated, and does esfixedbaseperiod
*      "fix" it?  Measured over 300 draws (their evidence was 25 and ~49),
*      with e(df_r), e(N) and K captured so a small-sample-df explanation can
*      be tested rather than guessed.
*  Q1b: the same on a panel with 5 sub-units per state-year (i still = state,
*      cluster(st) still == i).  If the "fix" is a df artifact it dies here.
*  Q2 (MCBP-JOINTSIZE): size of the ANALYTIC joint pre-trend test e(joint_horiz)
*      and of e(adversarial_pval_att), same draws.
* Run:  stata-mp -e do r01_analytic.do            (~15 min)
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

* pull the main-table att row out of e(wdidmainresults)
cap prog drop getatt
prog def getatt, rclass
  return scalar b  = .
  return scalar se = .
  return scalar p  = .
  cap mat TT = e(wdidmainresults)
  cap confirm matrix TT
  if _rc exit
  local rn : rownames TT
  local re : roweq TT
  local nr : word count `rn'
  forvalues z = 1/`nr' {
    local nz : word `z' of `rn'
    local ez : word `z' of `re'
    if "`nz'"=="att" & "`ez'"=="main" {
      return scalar b  = TT[`z',1]
      return scalar se = TT[`z',2]
      return scalar p  = TT[`z',4]
      continue, break
    }
  }
end

tempname P
postfile `P' int draw str4 spec double b double se double p double dfr double nn double kk ///
    double pjoint double padv int rc using r01_out.dta, replace

local NDRAW 300
forvalues d = 1/`NDRAW' {
  local sx = 430000 + `d'
  mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(1)
  tempfile dd
  qui save "`dd'", replace
  foreach s in S0 S1 S2 S3 S5 {
    qui use "`dd'", clear
    ereturn clear
    if "`s'"=="S0" cap noi qui wooldid y st t ttre, cluster(st) att
    if "`s'"=="S1" cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3)
    if "`s'"=="S2" cap noi qui wooldid y st t ttre, cluster(st) att esfixedbaseperiod
    if "`s'"=="S3" cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3) esfixedbaseperiod
    if "`s'"=="S5" cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial
    local rc = _rc
    local dfr = .
    local nn = .
    local kk = .
    local pj = .
    local pa = .
    if `rc'==0 {
      cap local dfr = e(df_r)
      cap local nn  = e(N)
      cap local kk  = colsof(e(b))
      cap local pj  = e(joint_horiz)
      cap local pa  = e(adversarial_pval_att)
    }
    qui getatt
    post `P' (`d') ("`s'") (r(b)) (r(se)) (r(p)) (`dfr') (`nn') (`kk') (`pj') (`pa') (`rc')
  }
}
postclose `P'

* ---- Q1b: 5 sub-units per state-year, i still = state, cluster(st) still == i
tempname Q
postfile `Q' int draw str4 spec double b double se double p double dfr double nn double kk int rc ///
    using r01b_out.dta, replace
forvalues d = 1/120 {
  local sx = 560000 + `d'
  mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(5)
  tempfile dd
  qui save "`dd'", replace
  foreach s in M0 M2 {
    qui use "`dd'", clear
    ereturn clear
    if "`s'"=="M0" cap noi qui wooldid y st t ttre, cluster(st) att
    if "`s'"=="M2" cap noi qui wooldid y st t ttre, cluster(st) att esfixedbaseperiod
    local rc = _rc
    local dfr = .
    local nn = .
    local kk = .
    if `rc'==0 {
      cap local dfr = e(df_r)
      cap local nn  = e(N)
      cap local kk  = colsof(e(b))
    }
    qui getatt
    post `Q' (`d') ("`s'") (r(b)) (r(se)) (r(p)) (`dfr') (`nn') (`kk') (`rc')
  }
}
postclose `Q'

*---------------------------- report --------------------------------
use r01_out.dta, clear
gen byte rej = (p <= .05 + 1e-9) if p < .
dis "=== R01 MAIN (G=10, 1 row per state-year, cluster(st)==i, sharp null) ==="
foreach s in S0 S1 S2 S3 S5 {
  qui sum b if spec=="`s'"
  local sd = r(sd)
  local nb = r(N)
  qui sum se if spec=="`s'"
  local mse = r(mean)
  qui sum rej if spec=="`s'"
  local rj = r(mean)
  qui sum dfr if spec=="`s'"
  local md = r(mean)
  qui sum kk if spec=="`s'"
  local mk = r(mean)
  qui sum nn if spec=="`s'"
  local mn = r(mean)
  dis "R01 `s': n=`nb' sd(b)=" %7.4f `sd' " mse=" %7.4f `mse' " se/sd=" %6.3f `=`mse'/`sd'' ///
      " rej=" %5.3f `rj' " df_r=" %7.2f `md' " N=" %6.1f `mn' " K=" %5.1f `mk'
}
qui corr b if spec=="S0"
dis "--- per-draw correlations of att across specs ---"
qui keep draw spec b se
qui reshape wide b se, i(draw) j(spec) string
corr bS0 bS1 bS2 bS3
dis "R01 sqrt((N-K_S0)/(N-K_S2)) test: see K/N above"

use r01_out.dta, clear
gen byte rj = (pjoint <= .05+1e-9) if pjoint<.
gen byte ra = (padv   <= .05+1e-9) if padv<.
qui sum rj
dis "R01 JOINT pre-trend rejection @5% = " %6.4f r(mean) "  n=" r(N)
qui sum ra
dis "R01 ADVERSARIAL rejection @5% = " %6.4f r(mean) "  n=" r(N)
foreach a in .01 .05 .10 .20 .30 .50 {
  qui count if pjoint <= `a' + 1e-9 & pjoint <.
  local c1 = r(N)
  qui count if pjoint <.
  dis "R01 joint ECDF at `a' = " %6.4f `=`c1'/r(N)'
}

use r01b_out.dta, clear
gen byte rej = (p <= .05 + 1e-9) if p < .
dis "=== R01b MULTI-UNIT (5 sub-units per state-year, i=state, cluster(st)==i) ==="
foreach s in M0 M2 {
  qui sum b if spec=="`s'"
  local sd = r(sd)
  local nb = r(N)
  qui sum se if spec=="`s'"
  local mse = r(mean)
  qui sum rej if spec=="`s'"
  local rj = r(mean)
  qui sum dfr if spec=="`s'"
  local md = r(mean)
  qui sum kk if spec=="`s'"
  local mk = r(mean)
  qui sum nn if spec=="`s'"
  local mn = r(mean)
  dis "R01b `s': n=`nb' sd(b)=" %7.4f `sd' " mse=" %7.4f `mse' " se/sd=" %6.3f `=`mse'/`sd'' ///
      " rej=" %5.3f `rj' " df_r=" %7.2f `md' " N=" %6.1f `mn' " K=" %5.1f `mk'
}
dis "R01_DONE"
