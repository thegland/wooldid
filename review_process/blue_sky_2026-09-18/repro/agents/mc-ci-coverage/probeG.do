clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
include "$WD/shapelib.doi"
use "$WD/fixA.dta", clear

* ---- permute: shift interval (perm_ci) and permtci, level honouring & shape ----
foreach L in 90 95 99 87.5 {
  local Lt = subinstr("`L'",".","p",.)
  cap noi wooldid y cg t ttre, att att_i cluster(unit) permute(shuffle) permreps(199) ///
        rseed(4711) permtstat permtci level(`L') esprelength(2) espostlength(2) dots(0)
  dis "RC|G_`Lt'|" _rc
  if _rc == 0 {
    dis "PERMLVL|`Lt'|elevel=" e(level) "|perm_reps=" e(perm_completed)
    wdck G_`Lt' e(wdidmainresults)
    mat G`Lt' = e(wdidmainresults)
    dis "PERMCOLS|`Lt'|" `"`: colnames e(wdidmainresults)'"'
    foreach e in att att_i {
      dis "PERMSHIFT|`Lt'|`e'|b=" %12.8f e(perm_obs_`e') "|shift=[" %12.8f e(perm_lb_`e') "," %12.8f e(perm_ub_`e') "]|tci=[" %12.8f e(perm_tlb_`e') "," %12.8f e(perm_tub_`e') "]|se=" %12.8f e(perm_se_`e')
    }
  }
}
wdnest PERM_90in99 G90 G99

* ---- permute without permtstat (pure shift interval) ----
cap noi wooldid y cg t ttre, att att_i cluster(unit) permute(shuffle) permreps(199) ///
      rseed(4711) level(95) esprelength(2) espostlength(2) dots(0)
dis "RC|G_plain|" _rc
if _rc==0 {
  wdck G_plain e(wdidmainresults)
  foreach e in att att_i {
    dis "PERMSHIFT|plain|`e'|b=" %12.8f e(perm_obs_`e') "|shift=[" %12.8f e(perm_lb_`e') "," %12.8f e(perm_ub_`e') "]"
  }
}
dis "PROBEG_DONE"

* ---- extra: impvalues at the help-RECOMMENDED pairing i=cohort, cluster=unit ----
use "$WD/fixA.dta", clear
cap noi wooldid y cg t ttre, att cluster(unit) impvalues level(95)
dis "RC|G_imp_iunitclus|" _rc
if _rc==0 {
  mat Q = e(wdidmainresults)
  dis "IMPX|icg_clusunit|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2] "|lb=" %14.9f Q[1,5] "|ub=" %14.9f Q[1,6] "|df=" e(df_r)
}
* same data, homogeneous effects: rebuild truetau constant
replace y = y - truetau
replace truetau = cond(ttre<. & t>=ttre, 1, 0)
replace y = y + truetau
cap noi wooldid y cg t ttre, att cluster(unit) impvalues level(95)
dis "RC|G_imp_homog|" _rc
if _rc==0 {
  mat Q = e(wdidmainresults)
  dis "IMPX|homog_icg_clusunit|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2] "|lb=" %14.9f Q[1,5] "|ub=" %14.9f Q[1,6] "|df=" e(df_r)
}
cap noi wooldid y cg t ttre, att cluster(unit) level(95)
if _rc==0 {
  mat Q = e(wdidmainresults)
  dis "IMPX|homog_analytic|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2]
}
dis "PROBEG2_DONE"
