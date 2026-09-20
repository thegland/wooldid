* ===== two controls the advisor flagged =====
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid

dis "########## 1. EC-2 mechanism: how does a Stata `real' render? ##########"
foreach L in 90 95 99 97.5 99.5 68.5 99.9 90.1 99.99 99.95 68.3 92.5 {
  local x = `L'
  dis "LVLFMT|`L'|macro=`x'|%18.0g=" %18.0g `L' "|%12.0g=" %12.0g `L' "|hex=" %21x `L'
  cap set level `x'
  dis "   set level `x' (via a macro-expanded real) -> rc " _rc
  set level 95
}

dis "########## 2. EC-4 control: the design the help says SATISFIES its conditions ##########"
* the help: "a coarser i (cohorts, or treatment-onset groups) with unit fixed effects
* supplied through fe() is the design that satisfies them"
clear
set obs 120
gen long id = _n
gen double ttre = .
replace ttre = 4  if inrange(id,1,30)
replace ttre = 7  if inrange(id,31,60)
replace ttre = 10 if inrange(id,61,80)
gen long coh = 1 + (ttre==7) + 2*(ttre==10)
replace coh = 4 if missing(ttre)
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()
save S3.dta, replace

cap prog drop SH
prog def SH
  args tag
  cap confirm matrix e(V_es_att)
  if _rc {
    dis "SH|`tag'|no V_es_att"
    exit
  }
  local rk = .
  local ev = ""
  local md = .
  local ch = .
  cap mata: st_local("rk", strofreal(rank(st_matrix("e(V_es_att)"))))
  cap mata: st_local("ev", invtokens(strofreal(sort(Re(eigenvalues(st_matrix("e(V_es_att)")))',1)',"%9.3e")))
  cap mata: st_local("md", strofreal(min(diagonal(st_matrix("e(V_es_att)"))),"%9.3e"))
  cap mata: st_local("ch", strofreal(sum(cholesky(makesymmetric(st_matrix("e(V_es_att)"))) :>= .)))
  local nd = .
  cap local nd = e(es_ndegen_att)
  local rkv = .
  cap mata: st_local("rkv", strofreal(rank(st_matrix("e(V)"))))
  dis "SH|`tag'|H=" colsof(e(V_es_att)) "|rank=`rk'|mindiag=`md'|cholmiss=`ch'|es_ndegen=`nd'|K=" colsof(e(b)) "|rank(e(V))=`rkv'"
  dis "SHEIG|`tag'|`ev'"
end

* i = unit  (the geometry the help flags)
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) cluster(clu)
dis "U_RC=" _rc
SH i_is_unit_cluster
* i = cohort + fe(id)   (the geometry the help endorses), clustered
cap noi wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) fe(id) cluster(clu)
dis "V_RC=" _rc
SH i_is_cohort_fe_cluster
* i = cohort + fe(id), robust
cap noi wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) fe(id) robust
dis "W_RC=" _rc
SH i_is_cohort_fe_robust
* i = cohort + fe(id), iid
cap noi wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) fe(id)
dis "X_RC=" _rc
SH i_is_cohort_fe_iid
* i = cluster (a coarser i that is exactly the cluster), clustered
cap noi wooldid y clu t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) fe(id) cluster(clu)
dis "Y_RC=" _rc
SH i_is_cluster

dis "########## 3. EC-6: are pweight SEs actually robust? ##########"
use S3.dta, clear
gen double pw = 10 + 40*runiform()
cap prog drop SE1
prog def SE1
  args tag
  tempname M
  cap confirm matrix e(wdidmainresults)
  if _rc exit
  matrix `M' = e(wdidmainresults)
  dis "SE|`tag'|att est=" %12.9f `M'[1,1] " se=" %12.9f `M'[1,2] "|vce=`e(vce)'|vcetype=`e(vcetype)'|wtype=`e(wtype)'|N=" e(N)
end
cap noi wooldid y id t ttre [pw=pw], att esfixedbaseperiod esprelength(3) espostlength(3)
dis "PW_RC=" _rc
SE1 pweight
cap noi wooldid y id t ttre [aw=pw], att esfixedbaseperiod esprelength(3) espostlength(3)
dis "AW_RC=" _rc
SE1 aweight
cap noi wooldid y id t ttre [aw=pw], att esfixedbaseperiod esprelength(3) espostlength(3) robust
dis "AWR_RC=" _rc
SE1 aweight_robust
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3)
dis "NW_RC=" _rc
SE1 noweight

dis "CONTROLS_DONE"
