*======================================================================
* v01_validate.do -- harness calibration BEFORE the big MC.
*  (1) 60 draws of the PLAIN analytic run at G=10 -> known-bad analytic
*      joint / adversarial oversize, and analytic att coverage
*  (2) e(perm_npossible) sanity: 10!/(2!2!2!4!) = 18900
*  (3) permdist bit-identity across p-value conventions at one rseed
*  (4) timing of the att+att_i bootstrap (reid regime) and permtstat
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

*---------------- DGP: SHARP NULL state-year panel (no treatment effect) -----
* ngrp   states (= i = cluster);  nper periods 2000..2000+nper-1
* onsets 2004/2005/2006 to 20% of states each; 40% never treated
* y = a_i + f_t + e_it ,  e_it AR(1) rho=.5 within state (Bertrand-Duflo-
* Mullainathan serial correlation).  hetv scales the error sd of TREATED
* states only (hetv=1 -> homogeneous).
cap prog drop mkdgp
prog def mkdgp
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1) HETV(real 1) NPERCELL(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `ngrp'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  local nc = round(`ngrp'*0.2)
  gen double ttre = .
  replace ttre = 2004 if rk <= `nc'
  replace ttre = 2005 if rk >    `nc' & rk <= 2*`nc'
  replace ttre = 2006 if rk >  2*`nc' & rk <= 3*`nc'
  gen double ai = rnormal()
  gen double esd = cond(ttre < ., `hetv', 1)
  drop u rk
  expand `nper'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()*esd if _n==1
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(1-0.25)*esd if _n>1
  gen double y = ai + ft + e
  if `npercell' > 1 {
    expand `npercell'
    bysort st t: gen int sub = _n
    gen long unit = st*1000 + sub
    replace y = y + rnormal()*0.5
  }
  keep st t ttre y
  order st t ttre y
end

*---------------- (1) analytic positive controls ----------------
tempname P
postfile `P' int draw double p_att double lo_att double hi_att ///
   double p_atti double lo_atti double hi_atti double p_joint double p_adv int rc ///
   using v01_analytic.dta, replace

forvalues d = 1/60 {
  local sd = 70000 + `d'
  mkdgp , ngrp(10) nper(8) seedx(`sd')
  ereturn clear
  cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) jointtests adversarial
  local rc = _rc
  local pa = .
  local la = .
  local ha = .
  local pi = .
  local li = .
  local hi = .
  local pj = .
  local pv = .
  if `rc' == 0 {
    mat M = e(wdidmainresults)
    local rn : rownames M
    local k = 0
    foreach r of local rn {
      local ++k
      if "`r'" == "att" {
        local pa = M[`k',4]
        local la = M[`k',5]
        local ha = M[`k',6]
      }
      if "`r'" == "att_i" {
        local pi = M[`k',4]
        local li = M[`k',5]
        local hi = M[`k',6]
      }
    }
    local pj = e(joint_horiz)
    local pv = e(adversarial_pval_att)
  }
  post `P' (`d') (`pa') (`la') (`ha') (`pi') (`li') (`hi') (`pj') (`pv') (`rc')
  if `d' == 1 {
    dis "FIRSTDRAW rc=`rc' p_att=`pa' p_joint=`pj' p_adv=`pv'"
    mat list M
    local rnn : rownames M
    dis "ROWNAMES=`rnn'"
    dis "NCLUST=" e(N_clust) " N=" e(N)
  }
}
postclose `P'

use v01_analytic.dta, clear
gen byte rej_att   = (p_att   <= 0.05 + 1e-9) if p_att   < .
gen byte rej_atti  = (p_atti  <= 0.05 + 1e-9) if p_atti  < .
gen byte rej_joint = (p_joint <= 0.05 + 1e-9) if p_joint < .
gen byte rej_adv   = (p_adv   <= 0.05 + 1e-9) if p_adv   < .
gen byte cov_att   = (lo_att  <= 0 & hi_att  >= 0) if lo_att < .
gen byte cov_atti  = (lo_atti <= 0 & hi_atti >= 0) if lo_atti < .
dis "=== ANALYTIC CONTROLS G=10, 60 draws ==="
sum rej_att rej_atti rej_joint rej_adv cov_att cov_atti
count if rc != 0
dis "NRC_NONZERO=" r(N)

*---------------- (2)+(3) permute sanity ----------------
mkdgp , ngrp(10) nper(8) seedx(70001)
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) permute(shuffle) permreps(79) rseed(555)
dis "PERM_RC=" _rc
dis "NPOSSIBLE=" e(perm_npossible) "  str=" e(perm_npossible_str)
dis "COMPLETED=" e(perm_completed) " FAILED=" e(perm_failed) " BLANK=" e(perm_blank)
dis "P_DEFAULT=" %21x e(perm_p_att)
dis "PERM_LB=" e(perm_lb_att) " PERM_UB=" e(perm_ub_att)
mat D0 = e(permdist_b)
dis "DIM_permdist_b=" rowsof(D0) "x" colsof(D0)
mat C0 = e(perm_comparator_b)
mat list C0

ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) permute(shuffle) permreps(79) rseed(555) permabsp
dis "PERMABSP_RC=" _rc
dis "P_ABSP=" %21x e(perm_p_att)
mat D1 = e(permdist_b)

ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) permute(shuffle) permreps(79) rseed(555) permrawp
dis "PERMRAWP_RC=" _rc
dis "P_RAWP=" %21x e(perm_p_att)
mat D2 = e(permdist_b)

mata: st_numscalar("dd1", max(abs(st_matrix("D0")-st_matrix("D1"))))
mata: st_numscalar("dd2", max(abs(st_matrix("D0")-st_matrix("D2"))))
dis "MAXABSDIFF_D0_D1=" %21x dd1
dis "MAXABSDIFF_D0_D2=" %21x dd2

*---------------- (4) timing: att+att_i bootstrap (reid regime), permtstat ----
mkdgp , ngrp(10) nper(8) seedx(70001)
timer clear
timer on 1
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(99) rseed(555)
timer off 1
dis "BOOT_ATTI_RC=" _rc " completed=" e(boot_completed) " failed=" e(boot_failed)
mat list e(boot_p_percentile)
mat list e(boot_ci_percentile)
mat list e(boot_p_bc)
mat list e(boot_ci_normal)

timer on 2
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(555) permtstat
timer off 2
dis "PERMT_RC=" _rc " completed=" e(perm_completed) " failed=" e(perm_failed) " blank=" e(perm_blank)
dis "PERMT_p_att=" e(perm_p_att) " p_att_i=" e(perm_p_att_i)

timer on 3
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(bayes) reps(99) rseed(555)
timer off 3
dis "BAYES_RC=" _rc " completed=" e(boot_completed)

timer on 4
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(99) rseed(555) boottstat
timer off 4
dis "BOOTT_RC=" _rc " completed=" e(boot_completed)
mat list e(boot_tpctl_mainresults)

timer list
