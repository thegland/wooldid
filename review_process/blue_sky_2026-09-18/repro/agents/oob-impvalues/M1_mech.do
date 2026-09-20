clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771031
set sortseed 771031
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/oob-impvalues/mkdgp.doi"

di as txt "SENTINEL_START"

*================= 0. harness sanity: noiseless target vs hand-computed =========
mkD 1 0 0 1
replace y = base
cap noi wooldid y cg t ttre, att cluster(unit) level(95)
mat Q = e(wdidmainresults)
di as txt "NOISELESS analytic att b = " as res %20.14f Q[1,1]
* hand-computed cell-weighted ATT on the treated: mean of truetau over treated post cells
qui su truetau if cg<=12 & t>=ttre
di as txt "HAND cell-weighted mean truetau (equal cell sizes) = " as res %20.14f r(mean)
cap noi wooldid y cg t ttre, att cluster(unit) impvalues level(95)
mat Q = e(wdidmainresults)
di as txt "NOISELESS impvalues b (cluster unit) = " as res %20.14f Q[1,1]
cap noi wooldid y cg t ttre, att cluster(cg) impvalues level(95)
mat Q = e(wdidmainresults)
di as txt "NOISELESS impvalues b (cluster cg)   = " as res %20.14f Q[1,1]
di as txt "SENTINEL_S0_DONE"

*================= 1. mechanism: algebraic identity across geometries ==========
* K = 12 treated i-levels, n = 5 units each -> predicted se ratio sqrt((K-1)/(nK-1))
local K 12
local n 5
local predratio = sqrt((`K'-1)/(`n'*`K'-1))
di as txt "PREDICTED se_unit/se_cg = " as res %12.9f `predratio'

foreach s in 4101 4102 4103 {
  mkD `s' 0.60 0.15 1
  di as txt "---- dataset seed `s' ----"
  foreach cv in unit cg reg super {
    cap noi wooldid y cg t ttre, att cluster(`cv') impvalues level(95)
    if _rc {
      di as txt "MECH|`s'|`cv'|RC=" _rc
    }
    else {
      mat Q = e(wdidmainresults)
      di as txt "MECH|`s'|`cv'|b=" %20.14f Q[1,1] "|se=" %20.14f Q[1,2] "|dfr=" e(df_r) "|N_clust=" e(N_clust)
    }
    cap noi wooldid y cg t ttre, att cluster(`cv') level(95)
    if !_rc {
      mat Q = e(wdidmainresults)
      di as txt "ANLY|`s'|`cv'|b=" %20.14f Q[1,1] "|se=" %20.14f Q[1,2] "|dfr=" e(df_r)
    }
  }
}
di as txt "SENTINEL_S1_DONE"

*================= 2. noise sweep: does the IM se carry the sampling component? =
foreach ns in 0.20 0.60 2.00 {
  mkD 5200 `ns' 0.15 1
  foreach cv in unit cg {
    cap noi wooldid y cg t ttre, att cluster(`cv') impvalues level(95)
    if !_rc {
      mat Q = e(wdidmainresults)
      di as txt "SWEEP_IM|ns=`ns'|`cv'|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2]
    }
    cap noi wooldid y cg t ttre, att cluster(`cv') level(95)
    if !_rc {
      mat Q = e(wdidmainresults)
      di as txt "SWEEP_AN|ns=`ns'|`cv'|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2]
    }
  }
}
di as txt "SENTINEL_S2_DONE"

*================= 3. what does the user SEE under geometry (a)? ===============
mkD 9001 0.60 0.15 1
di as txt "VVVVV BEGIN VERBATIM impvalues cluster(unit) OUTPUT VVVVV"
noi wooldid y cg t ttre, att cluster(unit) impvalues level(95)
di as txt "^^^^^ END VERBATIM ^^^^^"
di as txt "e(wdid_vcespec) = " as res `"`e(wdid_vcespec)'"'
di as txt "SENTINEL_S3_DONE"
di as txt "SENTINEL_ALLDONE"
