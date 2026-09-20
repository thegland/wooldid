clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 98221
set sortseed 98221

which wooldid

local R = 120
tempname pf
capture erase mcres.dta
postfile `pf' int rep double bA seA lA uA double bB seB lB uB double bC seC lC uC ///
    using mcres, replace

forvalues r = 1/`R' {
  quietly {
    clear
    set obs 60
    gen unit = _n
    gen cg = ceil(_n/15)
    gen ttre = .
    replace ttre = 5 if cg==1
    replace ttre = 7 if cg==2
    replace ttre = 9 if cg==3
    gen double ui = rnormal(0,0.80)
    gen double uhet = rnormal(0,0.30)
    expand 12
    bysort unit: gen t = _n
    gen double te = 0.05*t
    gen double eps = rnormal(0,0.60)
    gen double effc = 0
    replace effc = 1.00 if cg==1
    replace effc = 1.20 if cg==2
    replace effc = 0.80 if cg==3
    gen byte post = !missing(ttre) & t>=ttre
    gen double y = 2 + ui + te + post*(effc + uhet) + eps
  }

  foreach arm in A B C {
    scalar b`arm' = .
    scalar s`arm' = .
    scalar l`arm' = .
    scalar u`arm' = .
  }

  capture quietly wooldid y cg t ttre, att cluster(unit) impvalues
  if _rc==0 {
    matrix MM = e(wdidmainresults)
    scalar bA = MM[1,1]
    scalar sA = MM[1,2]
    scalar lA = MM[1,5]
    scalar uA = MM[1,6]
  }
  capture quietly wooldid y cg t ttre, att cluster(cg) impvalues
  if _rc==0 {
    matrix MM = e(wdidmainresults)
    scalar bB = MM[1,1]
    scalar sB = MM[1,2]
    scalar lB = MM[1,5]
    scalar uB = MM[1,6]
  }
  capture quietly wooldid y cg t ttre, att cluster(unit)
  if _rc==0 {
    matrix MM = e(wdidmainresults)
    scalar bC = MM[1,1]
    scalar sC = MM[1,2]
    scalar lC = MM[1,5]
    scalar uC = MM[1,6]
  }

  post `pf' (`r') (bA) (sA) (lA) (uA) (bB) (sB) (lB) (uB) (bC) (sC) (lC) (uC)
  if mod(`r',10)==0 dis "PROGRESS rep=`r'"
}
postclose `pf'

use mcres, clear
count
* targets
scalar tgt_IM  = (1.00+1.20+0.80)/3
scalar tgt_att = (8*1.00 + 6*1.20 + 4*0.80)/18
dis "TARGET_IM  = " tgt_IM
dis "TARGET_ATT = " tgt_att

foreach arm in A B C {
  gen byte cov_own`arm' = .
  gen byte cov_tgt`arm' = .
}
foreach arm in A B {
  quietly sum b`arm'
  scalar mn`arm' = r(mean)
  quietly replace cov_own`arm' = (l`arm'<=mn`arm' & u`arm'>=mn`arm') if !missing(l`arm')
  quietly replace cov_tgt`arm' = (l`arm'<=tgt_IM & u`arm'>=tgt_IM) if !missing(l`arm')
}
quietly sum bC
scalar mnC = r(mean)
quietly replace cov_ownC = (lC<=mnC & uC>=mnC) if !missing(lC)
quietly replace cov_tgtC = (lC<=tgt_att & uC>=tgt_att) if !missing(lC)

dis "================ RESULTS ================"
foreach arm in A B C {
  quietly count if !missing(s`arm')
  local nn = r(N)
  quietly sum b`arm'
  local mb = r(mean)
  local sdb = r(sd)
  quietly sum s`arm'
  local ms = r(mean)
  quietly sum cov_own`arm'
  local co = r(mean)
  quietly sum cov_tgt`arm'
  local ct = r(mean)
  dis "ARM `arm': n=`nn'  mean(b)=" %9.5f `mb' "  sd(b)=" %9.5f `sdb' ///
      "  mean(se)=" %9.5f `ms' "  se/sd=" %7.4f `ms'/`sdb' ///
      "  cov_ownmean=" %6.3f `co' "  cov_target=" %6.3f `ct'
}
dis "3*binomSE at .95, R=" _N " : " 3*sqrt(.95*.05/_N)
dis "DONE_MC"
