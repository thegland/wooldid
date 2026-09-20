*==============================================================================
* LENS 0 (independent reproduction) of FP-SEMIELAS-FE
* Built from the CLAIM ALONE. Own DGP: an administrative-shaped county panel.
* Blocks:
*   A  assert the pinned build
*   B  main design, 4 arms: {default, customvce(robust)} x {dydx, eydx}
*   C  homogeneous-level control arm (claim says the gap must vanish)
*   D  mechanism WITHOUT wooldid + hand oracles
*   E  NOVEL: within-DEFAULT-ENGINE reparameterization, fe(sector) vs controls(dummies)
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
set linesize 200

*------------------------------------------------------------------ A: build
capture noisily which wooldid
dis "A_RC=" _rc

*------------------------------------------------------------------ DGP
* 60 counties x 14 years. Heterogeneous unit levels (the ordinary reason one
* includes unit FE at all), a gentle calendar trend, a time-varying 3-level
* sector with real teeth, three adoption waves, 15 never-treated counties.
capture noisily {
  clear
  set obs 60
  gen i = _n
  gen double alpha = 380 + 120*invnormal(uniform())       // county base level
  replace alpha = 140 if alpha < 140
  gen byte wave = mod(_n,4)                                // 0 = never treated
  gen ttre = .
  replace ttre = 2012 if wave==1
  replace ttre = 2015 if wave==2
  replace ttre = 2018 if wave==3
  expand 14
  bys i: gen t = 2007 + _n
  gen double gamma = 6*(t-2008)                            // calendar trend
  gen byte sector = 1 + mod(i + 2*t, 3)                    // varies WITHIN unit
  gen double seceff = cond(sector==2, 60, cond(sector==3, -40, 0))
  gen byte sec2 = sector==2
  gen byte sec3 = sector==3
  gen byte post = (ttre<.) & (t>=ttre)
  gen double y = alpha + gamma + seceff + 25*post + 10*invnormal(uniform())
  qui su y
  dis "DGP y: min=" r(min) " mean=" r(mean) " max=" r(max)
  qui count if y<=0
  dis "DGP nonpositive y = " r(N)
  save "main.dta", replace
}
dis "DGP_RC=" _rc

*------------------------------------------------------------------ B: 4 arms
capture noisily {
  use "main.dta", clear
  qui wooldid y i t ttre, att
  local lev_def = e(b)[1,1]
  qui wooldid y i t ttre, att customvce(robust)
  local lev_exp = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity
  local sem_def = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity customvce(robust)
  local sem_exp = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity unconditionalse safety(off)
  local sem_unc = e(b)[1,1]
  dis "B_LEV_DEF = " %21x `lev_def'
  dis "B_LEV_EXP = " %21x `lev_exp'
  dis "B_LEV_RATIO = " %20.12f `lev_def'/`lev_exp'
  dis "B_SEM_DEF = " %20.12f `sem_def'
  dis "B_SEM_EXP = " %20.12f `sem_exp'
  dis "B_SEM_UNC = " %20.12f `sem_unc'
  dis "B_SEM_RATIO = " %20.12f `sem_def'/`sem_exp'
}
dis "B_RC=" _rc

*------------------------------------------------------------------ B2: what the user SEES
* Does anything at all warn?  Print the full default semielasticity run noisily.
capture noisily {
  use "main.dta", clear
  noi wooldid y i t ttre, att semielasticity
}
dis "B2_RC=" _rc

*------------------------------------------------------------------ C: control arm
* Homogeneous unit levels, no calendar trend, no sector effect: the absorbed FE
* are then ~0 and xb ~= yhat.  Claim predicts the gap COLLAPSES here.
capture noisily {
  use "main.dta", clear
  replace y = 400 + 25*post + 10*invnormal(uniform())
  qui wooldid y i t ttre, att semielasticity
  local c_def = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity customvce(robust)
  local c_exp = e(b)[1,1]
  dis "C_SEM_DEF = " %20.12f `c_def'
  dis "C_SEM_EXP = " %20.12f `c_exp'
  dis "C_SEM_RATIO = " %20.12f `c_def'/`c_exp'
}
dis "C_RC=" _rc

*------------------------------------------------------------------ D: mechanism, no wooldid
capture noisily {
  use "main.dta", clear
  qui reghdfe y post, absorb(i.i i.t) tol(1e-9)
  local b_h = _b[post]
  qui margins, dydx(post) nose
  local d_h = r(b)[1,1]
  qui margins, eydx(post) nose
  local e_h = r(b)[1,1]
  qui predict double xb_h, xb
  qui reg y i.i i.t post
  local b_f = _b[post]
  qui margins, dydx(post) nose
  local d_f = r(b)[1,1]
  qui margins, eydx(post) nose
  local e_f = r(b)[1,1]
  qui predict double xb_f, xb
  dis "D_B_HDFE = " %21x `b_h'
  dis "D_B_FULL = " %21x `b_f'
  dis "D_DYDX_HDFE = " %20.12f `d_h'
  dis "D_DYDX_FULL = " %20.12f `d_f'
  dis "D_EYDX_HDFE = " %20.12f `e_h'
  dis "D_EYDX_FULL = " %20.12f `e_f'
  * hand oracles: eydx(x) = mean( (dyhat/dx) / yhat )
  qui gen double o_h = `b_h'/xb_h
  qui gen double o_f = `b_f'/xb_f
  qui su o_h, meanonly
  dis "D_ORACLE_HDFE = " %20.12f r(mean)
  qui su o_f, meanonly
  dis "D_ORACLE_FULL = " %20.12f r(mean)
  qui corr xb_h xb_f
  dis "D_CORR_XB = " %20.12f r(rho)
  qui su xb_h
  dis "D_SD_XB_HDFE = " %20.12f r(sd) "   mean=" %20.12f r(mean)
  qui su xb_f
  dis "D_SD_XB_FULL = " %20.12f r(sd) "   mean=" %20.12f r(mean)
  qui su y, meanonly
  dis "D_MEAN_Y = " %20.12f r(mean)
  * arithmetic check: is the hdfe denominator just the grand mean?
  qui su post, meanonly
  local pbar = r(mean)
  qui su y, meanonly
  dis "D_GRANDMEAN_CONS = " %20.12f r(mean)-`b_h'*`pbar'
}
dis "D_RC=" _rc

*------------------------------------------------------------------ E: NOVEL
* Same DEFAULT reghdfe engine, same column space, two parameterizations of the
* sector effect.  fe(sector) absorbs it; controls(sec2 sec3) keeps it in e(b).
* Identical model => every reported quantity must be identical.
capture noisily {
  use "main.dta", clear
  * sanity: confirm the two fits ARE the same model
  qui reghdfe y post sec2 sec3, absorb(i.i i.t) tol(1e-9)
  dis "E_CHK_CTRL_b_post = " %21x _b[post] "  rank=" e(rank) "  r2=" %20.14f e(r2)
  qui reghdfe y post, absorb(i.i i.t sector) tol(1e-9)
  dis "E_CHK_FE_b_post   = " %21x _b[post] "  rank=" e(rank) "  r2=" %20.14f e(r2)

  qui wooldid y i t ttre, att fe(sector)
  local e_lev_fe = e(b)[1,1]
  qui wooldid y i t ttre, att controls(sec2 sec3)
  local e_lev_ct = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity fe(sector)
  local e_sem_fe = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity controls(sec2 sec3)
  local e_sem_ct = e(b)[1,1]
  dis "E_LEV_FE   = " %21x `e_lev_fe'
  dis "E_LEV_CTRL = " %21x `e_lev_ct'
  dis "E_LEV_RATIO = " %20.12f `e_lev_fe'/`e_lev_ct'
  dis "E_SEM_FE   = " %20.12f `e_sem_fe'
  dis "E_SEM_CTRL = " %20.12f `e_sem_ct'
  dis "E_SEM_RATIO = " %20.12f `e_sem_fe'/`e_sem_ct'
  * and the fully-correct reference for this model
  qui wooldid y i t ttre, att semielasticity controls(sec2 sec3) customvce(robust)
  dis "E_SEM_CTRL_EXPLICIT = " %20.12f e(b)[1,1]
}
dis "E_RC=" _rc

*------------------------------------------------------------------ F: reach
* Does it move an event-study horizon and a pretreatment coefficient too?
capture noisily {
  use "main.dta", clear
  qui wooldid y i t ttre, att semielasticity esprelength(3) espostlength(3)
  matrix Fd = e(b)
  qui wooldid y i t ttre, att semielasticity esprelength(3) espostlength(3) customvce(robust)
  matrix Fe = e(b)
  local nc = colsof(Fd)
  dis "F_NCOLS=" `nc'
  forvalues k = 1/`nc' {
    local nm : word `k' of `: colnames Fd'
    dis "F_ROW " %-28s "`nm'" "  def=" %16.9f Fd[1,`k'] "  exp=" %16.9f Fe[1,`k'] "  ratio=" %12.6f Fd[1,`k']/Fe[1,`k']
  }
}
dis "F_RC=" _rc

dis "LENS0_DONE"
