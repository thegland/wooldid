clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
include "$WD/shapelib.doi"
use "$WD/fixA.dta", clear

* ---- H1: permute with cluster AT i (cluster nests i trivially) ----
foreach L in 90 95 99 87.5 {
  local Lt = subinstr("`L'",".","p",.)
  cap noi wooldid y cg t ttre, att att_i cluster(cg) permute(shuffle) permreps(199) ///
        rseed(4711) permtstat permtci level(`L') esprelength(2) espostlength(2) dots(0)
  dis "RC|H1_`Lt'|" _rc
  if _rc == 0 {
    dis "PERMLVL|H1_`Lt'|elevel=" e(level) "|completed=" e(perm_completed)
    wdck H1_`Lt' e(wdidmainresults)
    mat H`Lt' = e(wdidmainresults)
    dis "PERMCOLS|H1_`Lt'|" `"`: colnames e(wdidmainresults)'"'
    foreach e in att att_i {
      dis "PERMSHIFT|H1_`Lt'|`e'|b=" %12.8f e(perm_obs_`e') "|shift=[" %12.8f e(perm_lb_`e') "," %12.8f e(perm_ub_`e') "]|tci=[" %12.8f e(perm_tlb_`e') "," %12.8f e(perm_tub_`e') "]|se=" %12.8f e(perm_se_`e')
    }
  }
}
wdnest H1_90in99 H90 H99

* ---- H2: permonsetcohorts (finer-than-i cluster, the 2026-09-16 feature) ----
foreach L in 90 95 99 {
  cap noi wooldid y cg t ttre, att att_i cluster(unit) permute(shuffle) permreps(199) ///
        rseed(4711) permtstat permtci permonsetcohorts level(`L') esprelength(2) espostlength(2) dots(0)
  dis "RC|H2_`L'|" _rc
  if _rc == 0 {
    dis "PERMLVL|H2_`L'|elevel=" e(level) "|completed=" e(perm_completed) "|irebuilt=" e(perm_irebuilt)
    wdck H2_`L' e(wdidmainresults)
    mat J`L' = e(wdidmainresults)
    foreach e in att att_i {
      dis "PERMSHIFT|H2_`L'|`e'|b=" %12.8f e(perm_obs_`e') "|shift=[" %12.8f e(perm_lb_`e') "," %12.8f e(perm_ub_`e') "]|tci=[" %12.8f e(perm_tlb_`e') "," %12.8f e(perm_tub_`e') "]|se=" %12.8f e(perm_se_`e')
    }
  }
}
wdnest H2_90in99 J90 J99
dis "PROBEH_DONE"
