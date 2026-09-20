clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 880001
set sortseed 880001
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/oob-impvalues/mkdgp.doi"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/oob-impvalues"

di as txt "SENTINEL_START"

*---- targets: the estimand each arm actually reports, from NOISELESS data ----
mkD 1 0 0 1
replace y = base
tempfile nf
qui save `nf'
foreach cv in unit cg reg super {
  use `nf', clear
  cap noi qui wooldid y cg t ttre, att cluster(`cv') impvalues level(95)
  mat Q = e(wdidmainresults)
  local TIM_`cv' = Q[1,1]
  use `nf', clear
  cap noi qui wooldid y cg t ttre, att cluster(`cv') level(95)
  mat Q = e(wdidmainresults)
  local TAN_`cv' = Q[1,1]
  di as txt "TARGET|`cv'|IM=" %18.12f `TIM_`cv'' "|AN=" %18.12f `TAN_`cv''
}
di as txt "SENTINEL_TARGETS_DONE"

tempname pf
postfile `pf' int rep str6 geom str3 arm double(b se ll ul) int rc using "`WD'/mcres.dta", replace

local R 200
timer clear 1
timer on 1
forvalues r = 1/`R' {
  mkD `=880000+`r'' 0.60 0.15 1
  tempfile dr
  qui save `dr'
  foreach cv in unit cg reg super {
    use `dr', clear
    cap qui wooldid y cg t ttre, att cluster(`cv') impvalues level(95)
    if _rc post `pf' (`r') ("`cv'") ("IM") (.) (.) (.) (.) (_rc)
    else {
      mat Q = e(wdidmainresults)
      post `pf' (`r') ("`cv'") ("IM") (Q[1,1]) (Q[1,2]) (Q[1,5]) (Q[1,6]) (0)
    }
    use `dr', clear
    cap qui wooldid y cg t ttre, att cluster(`cv') level(95)
    if _rc post `pf' (`r') ("`cv'") ("AN") (.) (.) (.) (.) (_rc)
    else {
      mat Q = e(wdidmainresults)
      post `pf' (`r') ("`cv'") ("AN") (Q[1,1]) (Q[1,2]) (Q[1,5]) (Q[1,6]) (0)
    }
  }
  if mod(`r',20)==0 {
    timer off 1
    qui timer list 1
    di as txt "PROGRESS rep `r' of `R'   elapsed=" r(t1)
    timer on 1
  }
}
timer off 1
postclose `pf'
di as txt "SENTINEL_MC_DONE"

*---- summarise ----
use "`WD'/mcres.dta", clear
gen double target = .
foreach cv in unit cg reg super {
  replace target = `TIM_`cv'' if geom=="`cv'" & arm=="IM"
  replace target = `TAN_`cv'' if geom=="`cv'" & arm=="AN"
}
gen byte cov = (ll<=target & target<=ul) if !missing(ll,ul)
foreach cv in unit cg reg super {
  foreach a in IM AN {
    qui count if geom=="`cv'" & arm=="`a'" & rc==0
    local nn = r(N)
    qui su b if geom=="`cv'" & arm=="`a'" & rc==0
    local mb = r(mean)
    local sdb = r(sd)
    qui su se if geom=="`cv'" & arm=="`a'" & rc==0
    local mse = r(mean)
    qui su cov if geom=="`cv'" & arm=="`a'" & rc==0
    local cv_ = r(mean)
    local secov = sqrt(`cv_'*(1-`cv_')/`nn')
    qui count if geom=="`cv'" & arm=="`a'" & rc!=0
    local nfail = r(N)
    di as txt "RES|`cv'|`a'|n=" `nn' "|fail=" `nfail' "|target=" %10.6f `=cond("`a'"=="IM",`TIM_`cv'',`TAN_`cv'')' "|meanb=" %10.6f `mb' "|sd_b=" %10.6f `sdb' "|mean_se=" %10.6f `mse' "|se_over_sd=" %8.4f `=`mse'/`sdb'' "|cov=" %6.3f `cv_' "|covSE=" %6.4f `secov'
  }
}
di as txt "SENTINEL_ALLDONE"
