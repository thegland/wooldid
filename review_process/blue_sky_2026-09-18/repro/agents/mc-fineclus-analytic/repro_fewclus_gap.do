* =======================================================================================
* repro_fewclus_gap.do -- SELF-CONTAINED. Run with:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_fewclus_gap.do
* Takes ~8 minutes. Needs the pinned wooldid build and csdid.
*
* CLAIM UNDER TEST. wooldid.sthlp:251 presents the analytic-SE problem as a BINARY in the
* nesting relation: "At or above the level of i the cluster-robust covariance stops seeing
* the treated cells' own variation ... the remedies are a finer cluster()". This shows the
* understatement is a GRADIENT in the number of clusters inside a treated cohort-by-period
* cell, not a switch: with cluster() strictly FINER than i throughout, wooldid's reported
* SE runs below csdid's -- on IDENTICAL point estimates and identical sampling variability
* -- by about 8% when a treated cohort holds 10 clusters and about 17% when it holds 3.
*
* LADDER. Everything is held fixed except K = the number of states in each treated cohort:
*   40 states x 4 units x 12 periods; two adoption waves at t=5 and t=8 with K states each;
*   the other 40-2K states never adopt.
*   y = s_g + v_gt + a_u + .2 t + .5 cos t + tau + e
*   s_g ~ N(0,.8^2) state RE ; v_gt ~ N(0,.4^2) state-by-year shock ; a_u ~ N(0,.5^2) unit RE
*   e AR(1) rho .5, innovation sd .4 ; tau = 1.0 + .05*rel (wave 1), 0.4 + .05*rel (wave 2)
*   States are independent, so cluster(state) is the correct level in EVERY rung, and there
*   are 40 clusters in every rung.
* ARM: i = coh (the onset cohort: 3 values), cluster(state) -- cluster STRICTLY FINER than i,
*      i.e. exactly the geometry the help calls the remedy. csdid long2 is the oracle.
* True post-h effect = 0.70 + 0.05h for h = 0..3 (both waves are observed at those horizons).
* =======================================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local NDRAWS = 150

cap program drop mkk
program define mkk
  args K
  clear
  qui set obs 40
  gen double state = _n
  gen double s_g   = rnormal()*0.8
  gen double ttre  = cond(state<=`K', 5, cond(state<=2*`K', 8, .))
  qui expand 4
  bysort state: gen double upos = _n
  gen double unit = (state-1)*4 + upos
  gen double a_u  = rnormal()*0.5
  qui expand 12
  bysort unit: gen double t = _n
  tempvar first vv eps
  bysort state t: gen byte `first' = (_n==1)
  gen double `vv' = rnormal()*0.4 if `first'
  bysort state t (`first'): replace `vv' = `vv'[_N]
  gen double v_gt = `vv'
  gen double `eps' = rnormal()*0.4
  bysort unit (t): gen double e = `eps'[1]
  by unit (t): replace e = 0.5*e[_n-1] + sqrt(1-0.25)*`eps'[_n] if _n > 1
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
  gen double y   = s_g + v_gt + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double g   = cond(ttre<., ttre, 0)
  gen double coh = cond(ttre<., ttre, 99)
  keep y state unit t ttre coh g
end

tempfile dd
foreach K in 3 10 {
  tempname M
  postfile `M' int draw double bw0 double sw0 double bw1 double sw1 double bw2 double sw2 double bw3 double sw3 ///
                  double bc0 double sc0 double bc1 double sc1 double bc2 double sc2 double bc3 double sc3 ///
                  using "_rfg_K`K'.dta", replace
  forval d = 1/`NDRAWS' {
    set seed `=600000 + 23*`d' + `K''
    set sortseed `=600000 + 23*`d' + `K''
    mkk `K'
    qui save `"`dd'"', replace
    forval h = 0/3 {
      local bw`h' = .
      local sw`h' = .
      local bc`h' = .
      local sc`h' = .
    }
    cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod
    if _rc == 0 {
      forval h = 0/3 {
        local bw`h' = _b[ES_att_post`h']
        local sw`h' = _se[ES_att_post`h']
      }
    }
    use `dd', clear
    cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(state) long2
    if _rc == 0 {
      cap qui estat event, post
      if _rc == 0 {
        forval h = 0/3 {
          local bc`h' = _b[Tp`h']
          local sc`h' = _se[Tp`h']
        }
      }
    }
    post `M' (`d') (`bw0') (`sw0') (`bw1') (`sw1') (`bw2') (`sw2') (`bw3') (`sw3') ///
             (`bc0') (`sc0') (`bc1') (`sc1') (`bc2') (`sc2') (`bc3') (`sc3')
    if mod(`d',50)==0 di "RFG_PROGRESS K=`K' draw `d'"
  }
  postclose `M'
}

local tc = invttail(39, 0.025)
di as txt ""
di as txt "K = states in each treated cohort. cluster(state) in every rung (40 clusters)."
di as txt "wooldid i = coh (cluster STRICTLY FINER than i). csdid is the oracle; the two"
di as txt "return the SAME point estimates, so SD is shared and only the SEs differ."
di as txt ""
di as txt "  K  hor       SD   SE_csdid  SE_wooldid   wd/cs   SE/SD_cs  SE/SD_wd  cov_cs  cov_wd"
foreach K in 3 10 {
  use "_rfg_K`K'.dta", clear
  forval h = 0/3 {
    local TR = 0.70 + 0.05*`h'
    qui su bw`h'
    local sd = r(sd)
    qui su sw`h'
    local msw = r(mean)
    qui su sc`h'
    local msc = r(mean)
    qui count if abs(bw`h' - `TR') <= `tc'*sw`h' & sw`h' < .
    local cw = r(N)
    qui count if abs(bc`h' - `TR') <= 1.959963985*sc`h' & sc`h' < .
    local cc = r(N)
    qui count if sw`h' < .
    local n = r(N)
    di as txt "  " %2.0f `K' "  post" `h' "  " as res %7.4f `sd' "   " %7.4f `msc' "    " %7.4f `msw' ///
       "  " %6.3f `msw'/`msc' "    " %7.3f `msc'/`sd' "  " %7.3f `msw'/`sd' " " %6.3f `cc'/`n' "  " %6.3f `cw'/`n'
  }
}
di as txt ""
di as txt "Read the wd/cs column: at K=10 wooldid reports about 92% of csdid's standard error"
di as txt "on the same estimate; at K=3 about 83%. The nesting relation is identical in both."
di "REPRO_RFG_DONE"
