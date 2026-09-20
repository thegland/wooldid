clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 660221
set sortseed 660221
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/oob-impvalues/mkdgp.doi"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/oob-impvalues"
di as txt "SENTINEL_START"

*========== A. NEAR-HOMOGENEOUS effects: does the IM interval get better or worse? ==========
* hetm = 0.15 shrinks the cross-cohort spread to ~0.045..0.24.
mkD 1 0 0 0.15
replace y = base
foreach cv in unit cg {
  cap noi qui wooldid y cg t ttre, att cluster(`cv') impvalues level(95)
  mat Q = e(wdidmainresults)
  local TIM_`cv' = Q[1,1]
  cap noi qui wooldid y cg t ttre, att cluster(`cv') level(95)
  mat Q = e(wdidmainresults)
  local TAN_`cv' = Q[1,1]
  di as txt "HOMTARGET|`cv'|IM=" %18.12f `TIM_`cv'' "|AN=" %18.12f `TAN_`cv''
}
tempname pf
postfile `pf' int rep str6 geom str3 arm double(b se ll ul) int rc using "`WD'/mcres_hom.dta", replace
local R 120
forvalues r = 1/`R' {
  mkD `=770000+`r'' 0.60 0.15 0.15
  tempfile dr
  qui save `dr'
  foreach cv in unit cg {
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
  if mod(`r',30)==0 di as txt "HOMPROGRESS `r'"
}
postclose `pf'
use "`WD'/mcres_hom.dta", clear
gen double target = .
foreach cv in unit cg {
  replace target = `TIM_`cv'' if geom=="`cv'" & arm=="IM"
  replace target = `TAN_`cv'' if geom=="`cv'" & arm=="AN"
}
gen byte cov = (ll<=target & target<=ul) if !missing(ll,ul)
foreach cv in unit cg {
  foreach a in IM AN {
    qui count if geom=="`cv'" & arm=="`a'" & rc==0
    local nn = r(N)
    qui su b if geom=="`cv'" & arm=="`a'" & rc==0
    local mb = r(mean)
    local sdb = r(sd)
    qui su se if geom=="`cv'" & arm=="`a'" & rc==0
    local mse = r(mean)
    qui su cov if geom=="`cv'" & arm=="`a'" & rc==0
    local cc = r(mean)
    di as txt "HOMRES|`cv'|`a'|n=" `nn' "|meanb=" %10.6f `mb' "|sd_b=" %10.6f `sdb' "|mean_se=" %10.6f `mse' "|se_over_sd=" %8.4f `=`mse'/`sdb'' "|cov=" %6.3f `cc' "|covSE=" %6.4f `=sqrt(`cc'*(1-`cc')/`nn')'
  }
}
di as txt "SENTINEL_A_DONE"

*========== B. identically-composed crossing clusters on NOISY data ==========
* 13 i-levels (12 treated over 4 onsets + 1 never-treated), 6 units each,
* region r = exactly one unit from each i-level -> every region has the SAME
* composition -> every per-cluster read is numerically identical.
clear
set obs 78
gen long unit = _n
gen byte  cg  = floor((_n-1)/6)+1
gen int   reg = mod(_n-1,6)+1
gen double ai = 0.8*mod(_n*23,31)/31 - 0.4
expand 10
bys unit: gen int t = 2000 + _n
gen double ttre = .
replace ttre = 2004 if inrange(cg,1,3)
replace ttre = 2005 if inrange(cg,4,6)
replace ttre = 2006 if inrange(cg,7,9)
replace ttre = 2007 if inrange(cg,10,12)
gen double truetau = 0
replace truetau = 0.30 + 0.12*(cg-1) if cg<=12 & t>=ttre
set seed 33127
set sortseed 33127
by unit: gen double us = rnormal(0,0.25) if _n==1
by unit: replace us = us[1]
gen double y = ai + 0.35*(t-2000) + truetau + us*(t-2005.5) + rnormal(0,0.90)
qui su y
di as txt "sd(y) = " as res r(sd) as txt "  (NOISY)"
di as txt "VVVVV BEGIN VERBATIM degenerate-crossing impvalues VVVVV"
cap noi wooldid y cg t ttre, att cluster(reg) impvalues level(95)
di as txt "^^^^^ END VERBATIM ^^^^^"
di as txt "VVVVV BEGIN VERBATIM same data, analytic cluster(reg) VVVVV"
cap noi wooldid y cg t ttre, att cluster(reg) level(95)
di as txt "^^^^^ END VERBATIM ^^^^^"
di as txt "SENTINEL_B_DONE"


*========== C. the honest analogue: fit the per-cluster estimates SEPARATELY ==========
* C1. cluster(unit) -- can a single unit support its own DiD estimate at all?
mkD 4101 0.60 0.15 1
tempfile full
qui save `full'
di as txt "VVVVV C1: wooldid on ONE unit's data (unit 3, a treated unit) VVVVV"
use `full', clear
keep if unit==3
cap noi wooldid y cg t ttre, att cluster(unit) level(95)
di as txt "C1_RC=" _rc
di as txt "^^^^^ END C1 ^^^^^"
di as txt "VVVVV C1b: one treated unit + all never-treated units VVVVV"
use `full', clear
keep if unit==3 | cg>12
cap noi wooldid y cg t ttre, att cluster(unit) level(95)
di as txt "C1b_RC=" _rc
if _rc==0 {
  mat Q = e(wdidmainresults)
  di as txt "C1b|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2]
}
di as txt "^^^^^ END C1b ^^^^^"

* C2. genuine per-COHORT separate fits (cohort c + the never-treated pool),
*     over 3 datasets, versus wooldid's own per-cohort reads from the pooled fit.
foreach s in 4101 4102 4103 {
  mkD `s' 0.60 0.15 1
  tempfile fx
  qui save `fx'
  * wooldid's internal read dispersion, recovered from the IM se: sd = se*sqrt(G)
  use `fx', clear
  cap qui wooldid y cg t ttre, att cluster(cg) impvalues level(95)
  mat Q = e(wdidmainresults)
  local wse = Q[1,2]
  local wb  = Q[1,1]
  local Geff = e(df_r)+1
  local wsd = `wse'*sqrt(`Geff')
  * genuine separate fits
  local ssum 0
  local ssq  0
  local nok  0
  forvalues c = 1/12 {
    use `fx', clear
    keep if cg==`c' | cg>12
    cap qui wooldid y cg t ttre, att cluster(unit) level(95)
    if _rc==0 {
      mat Q2 = e(wdidmainresults)
      local ssum = `ssum' + Q2[1,1]
      local ssq  = `ssq'  + Q2[1,1]^2
      local ++nok
    }
  }
  local sepmean = `ssum'/`nok'
  local sepsd = sqrt((`ssq' - `nok'*`sepmean'^2)/(`nok'-1))
  di as txt "SEP|`s'|nok=" `nok' "|wooldid_read_mean=" %12.8f `wb' "|wooldid_read_sd=" %12.8f `wsd' "|separatefit_mean=" %12.8f `sepmean' "|separatefit_sd=" %12.8f `sepsd' "|IMse_wooldid=" %12.8f `wse' "|IMse_if_separate=" %12.8f `=`sepsd'/sqrt(`nok')'
}
di as txt "SENTINEL_C_DONE"
di as txt "SENTINEL_ALLDONE"
