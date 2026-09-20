* =======================================================================================
* repro_remedy_works.do  -- SELF-CONTAINED. Run with:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_remedy_works.do
* Takes ~6 minutes. Needs only the pinned wooldid build and csdid.
*
* QUESTION. wooldid's help (wooldid.sthlp, "Standard Errors" / "What level to cluster at")
* says the analytic clustered SEs are understated when cluster() sits AT OR ABOVE the level
* of `i', "measured at roughly 0.66 to 0.75 of the reference", and names the remedy: "the
* remedies are a finer cluster(), or bootstrap() or permute()". Does the finer-cluster()
* remedy actually deliver correct coverage on the analytic path?
*
* DESIGN. 40 states x 5 units x 12 periods. Treatment onset is a STATE-level policy:
* states 1-10 adopt at t=5, states 11-20 adopt at t=8, states 21-40 never adopt.
*   y = s_g + v_gt + a_u + .2 t + .5 cos t + tau + e
*   s_g  ~ N(0,.8^2)  state random effect
*   v_gt ~ N(0,.4^2)  state-by-year shock  (the classic Bertrand-Duflo-Mullainathan term)
*   a_u  ~ N(0,.5^2)  unit random effect ;  e AR(1) rho .5, innovation sd .4, within unit
*   tau  = 1.0 + .05*rel for the t=5 adopters, 0.4 + .05*rel for the t=8 adopters
* Units are independent ACROSS states, so cluster(state) is exactly the right level.
* The true att (the cell-weighted ATT on the treated) is 0.9153844 by construction:
*   (8*(1.0+.05*3.5) + 5*(0.4+.05*2)) / 13 .
*
* THREE ARMS, on the SAME draws, all with cluster(state) and all with 40 clusters:
*   AT_I     i = state  -> cluster() sits AT the level of i  (the documented-bad regime)
*   FINER    i = coh    -> cluster() is strictly FINER than i (the documented REMEDY:
*                          each treated cohort-by-period cell holds 10 states)
*   CSDID    csdid long2 + estat event, cluster(state)      (external control)
* The point estimates are identical in the two wooldid arms, so only the variance differs.
* =======================================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local NDRAWS = 150
local TRUTH  = 0.9153843986204159844   // the cell-weighted att, exact

cap program drop mkd
program define mkd
  clear
  qui set obs 40
  gen double state = _n
  gen double s_g   = rnormal()*0.8
  gen double ttre  = cond(state<=10, 5, cond(state<=20, 8, .))
  qui expand 5
  bysort state: gen double unit = (state-1)*5 + _n
  gen double a_u   = rnormal()*0.5
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
  gen double coh = cond(ttre<., ttre, 99)
  gen double g   = cond(ttre<., ttre, 0)
  keep y state unit t ttre coh g
end

tempname M
postfile `M' int draw double b_ati double se_ati double b_fin double se_fin ///
              double b_cs double se_cs double bdiff using "_repro_rrw.dta", replace

tempfile dd
forval d = 1/`NDRAWS' {
  set seed `=900000 + 17*`d''
  set sortseed `=900000 + 17*`d''
  mkd
  qui save `"`dd'"', replace

  local b1 = .
  local s1 = .
  cap qui wooldid y state t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod
  if _rc == 0 {
    local b1 = _b[main_att]
    local s1 = _se[main_att]
  }

  use `dd', clear
  local b2 = .
  local s2 = .
  cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod
  if _rc == 0 {
    local b2 = _b[main_att]
    local s2 = _se[main_att]
  }

  use `dd', clear
  local b3 = .
  local s3 = .
  cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(state) long2
  if _rc == 0 {
    cap qui estat event, post
    if _rc == 0 {
      local b3 = _b[Tp0]
      local s3 = _se[Tp0]
    }
  }
  post `M' (`d') (`b1') (`s1') (`b2') (`s2') (`b3') (`s3') (abs(`b1'-`b2'))
  if mod(`d',25)==0 di "RRW_PROGRESS draw `d'"
}
postclose `M'

use "_repro_rrw.dta", clear
qui su bdiff
di as txt "-------------------------------------------------------------------"
di as txt "Largest |b(i=state) - b(i=coh)| over `NDRAWS' draws: " as res %12.0g r(max)
di as txt "  (cluster() and the choice of i must not move the point estimate)"
di as txt "-------------------------------------------------------------------"

* wooldid reports CIs on t(N_clust-1) = t(39); csdid on the normal.
local tc = invttail(39, 0.025)
foreach a in ati fin {
  qui su b_`a'
  local sd = r(sd)
  qui su se_`a'
  local mse = r(mean)
  qui count if abs(b_`a' - `TRUTH') <= `tc'*se_`a' & se_`a' < .
  local ncov = r(N)
  qui count if se_`a' < .
  local n = r(N)
  local nm "i = state, cluster(state)   [cluster AT i, the documented-bad regime]"
  if "`a'" == "fin" local nm "i = coh,   cluster(state)   [cluster FINER than i, THE REMEDY]"
  di as txt "`nm'"
  di as txt "   n = " as res `n' as txt "   empirical SD = " as res %6.4f `sd' ///
     as txt "   mean reported SE = " as res %6.4f `mse'
  di as txt "   mean(SE)/SD = " as res %6.3f `mse'/`sd' ///
     as txt "     95% coverage of the true att = " as res %6.3f `ncov'/`n'
}
qui su b_cs
local sdc = r(sd)
qui su se_cs
local msec = r(mean)
di as txt "csdid Tp0 (external control, cluster(state)):  mean(SE)/SD = " as res %6.3f `msec'/`sdc'
di as txt ""
di as txt "PASS RULE: mean(SE)/SD >= 0.90 and coverage >= 0.93."
di "REPRO_RRW_DONE"
