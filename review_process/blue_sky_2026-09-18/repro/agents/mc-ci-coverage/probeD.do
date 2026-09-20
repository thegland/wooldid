clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
set seed 20260916
set sortseed 20260916
include "$WD/shapelib.doi"

local FN 100
local FT 8
local FSEED 31
include "$WD/mkfix.doi"
gen int reg25 = mod(unit-1,25)+1
save "$WD/fixD.dta", replace

* ============ D1: impvalues (IM) at a cross-cutting cluster, G=25 ============
foreach L in 90 95 99 87.5 {
  local Lt = subinstr("`L'",".","p",.)
  cap noi wooldid y cg t ttre, att cluster(reg25) impvalues level(`L') esprelength(2) espostlength(2)
  dis "RC|D1_`Lt'|" _rc
  if _rc == 0 {
    dis "IMP|`Lt'|level=" e(level) "|df_r=" e(df_r) "|Nclust=" e(N_clust)
    wdck D1_`Lt' e(wdidmainresults)
    mat IM`Lt' = e(wdidmainresults)
    mat list IM`Lt', format(%14.9f)
    * implied critical value vs invttail(G-1)
    mat Q = e(wdidmainresults)
    forvalues i = 1/`=rowsof(Q)' {
      local b = Q[`i',1]
      local s = Q[`i',2]
      local lo = Q[`i',5]
      local hi = Q[`i',6]
      if `s'<. & `s'>0 & `lo'<. {
        local rn : rownames Q
        local rr : word `i' of `rn'
        dis "IMPCRIT|`Lt'|`rr'|z=" %12.9f ((`hi'-`lo')/(2*`s')) "|t24=" %12.9f invttail(24,(100-`L')/200) "|tn=" %12.9f invttail(e(df_r),(100-`L')/200)
      }
    }
  }
}
wdnest IMP_90in99 IM90 IM99

* ============ D2: permute + permtci ============
foreach L in 90 95 99 {
  cap noi wooldid y cg t ttre, att att_i cluster(unit) permute(clusters) permreps(199) ///
        rseed(4711) permtstat permtci level(`L') esprelength(2) espostlength(2) dots(0)
  dis "RC|D2_`L'|" _rc
  if _rc == 0 {
    wdck D2_`L' e(wdidmainresults)
    dis "PERMLVL|`L'|e(level)=" e(level)
    foreach e in att att_i {
      dis "PERMCI|`L'|`e'|shift=[" e(perm_lb_`e') "," e(perm_ub_`e') "]|tci=[" e(perm_tlb_`e') "," e(perm_tub_`e') "]|b=" e(perm_obs_`e')
    }
    cap confirm matrix e(perm_ci)
    if _rc==0 {
      mat list e(perm_ci), format(%14.9f)
    }
    if `L'==90 {
      mat P90 = e(wdidmainresults)
    }
    if `L'==99 {
      mat P99 = e(wdidmainresults)
    }
  }
}
wdnest PERM_90in99 P90 P99

* ============ D3: poisson / poisexpresults / semielasticity shape ============
gen double ycount = rpoisson(exp(0.5 + 0.2*(t-1999) + 0.4*(t>=ttre & ttre<.)))
cap noi wooldid ycount cg t ttre, att poisson cluster(unit) level(95) esprelength(2) espostlength(2)
dis "RC|D3_pois|" _rc
if _rc==0 {
  wdck D3_pois e(wdidmainresults)
}
cap noi wooldid ycount cg t ttre, att poisson poisexpresults cluster(unit) level(95) esprelength(2) espostlength(2)
dis "RC|D3_poisexp|" _rc
if _rc==0 {
  wdck D3_poisexp e(wdidmainresults)
  mat list e(wdidmainresults), format(%14.9f)
}
cap noi wooldid y cg t ttre, att semielasticity cluster(unit) level(95) esprelength(2) espostlength(2)
dis "RC|D3_semi|" _rc
if _rc==0 {
  wdck D3_semi e(wdidmainresults)
}
dis "PROBED_DONE"
