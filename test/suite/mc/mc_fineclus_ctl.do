* mc_fineclus_ctl.do - prereg F4/P6 and P8: the 2 - 1/G NEGATIVE CONTROL for the
* bootstrap re-identification leg, the only readout that shows the leg is working.
*
* A twice-drawn cluster has to count as TWO clusters in the replication's own
* cluster-robust standard error. The leg delivers that by substituting the resampled
* cluster identity __wdbns into the replication's vce; the substitution is a literal
* string edit that can silently no-op, and a no-op returns the full E[m^2]/E[m] =
* 2 - 1/G inflation the leg exists to remove. A run that merely produces plausible
* standard errors does not demonstrate the leg is working; a run that can reproduce
* the defect on demand and then not reproduce it does. WOOLDID_BOOTNOSUB, an internal
* undocumented diagnostic global, is what suppresses it.
*
* Two arms, both on dgpD (i = state, 12 states, 114 units):
*   nest  cluster(state) == i          G = 12,  predicted ratio sqrt(2 - 1/12)
*   fine  cluster(unit), finer than i  G = 114, predicted ratio sqrt(2 - 1/114)  [P8]
* PASS is |ratio/predicted - 1| <= .15, the tolerance PREREG_FEATURE_FINECLUS F0
* records as set without calibration.
*
* P7 rides along: wooldid asserts e(N_clust) against the distinct __wdbns count on
* every replication of the ACTIVE arms, so those arms failing rc != 0 IS the P7 result.
*
* usage: <stata> -e do mc_fineclus_ctl.do <ado>
clear all
set more off
set linesize 200
set graphics off
local ADO `1'
qui do `"`ADO'"'
run mc_fineclus_dgp.doh

cap prog drop fcarm
prog def fcarm, rclass
  * args: label  clustervar  reps  seed
  args lbl cv R sd
  qui dgpD 5150 1 1
  global WOOLDID_BOOTNOSUB
  cap noi wooldid y state t ttre, att cluster(`cv') bootstrap(standard) boottstat reps(`R') rseed(`sd') dots(0) bootsave(a_`lbl'.dta) replace
  local rcA = _rc
  local G = e(N_clust)


  qui dgpD 5150 1 1
  global WOOLDID_BOOTNOSUB 1
  cap noi wooldid y state t ttre, att cluster(`cv') bootstrap(standard) boottstat reps(`R') rseed(`sd') dots(0) bootsave(s_`lbl'.dta) replace
  local rcS = _rc
  global WOOLDID_BOOTNOSUB

  dis as text "FCARM=`lbl' cluster=`cv' rcA=`rcA' rcS=`rcS' G=`G'"
  if `rcA' != 0 | `rcS' != 0 exit

  use a_`lbl'.dta, clear
  qui sum se_main_att if _rep > 0
  local mA = r(mean)
  local nA = r(N)
  qui sum main_att if _rep > 0
  local nbA = r(N)
  use s_`lbl'.dta, clear
  qui sum se_main_att if _rep > 0
  local mS = r(mean)
  local nS = r(N)

  local pred = sqrt(2 - 1/`G')
  local ratio = `mS'/`mA'
  local relerr = `ratio'/`pred' - 1
  dis as text "FINECLUS_CTL|`lbl'|G=`G'|meanse_active=" %12.8f `mA' "|nA=`nA'|meanse_suppressed=" %12.8f `mS' "|nS=`nS'|ratio=" %9.6f `ratio' "|predicted=" %9.6f `pred' "|relerr=" %9.6f `relerr' "|verdict=" cond(abs(`relerr') <= .15, "PASS", "FAIL")
  return scalar ratio = `ratio'
end

* ---- P6: under NESTING (cluster == i == state)
fcarm nest state 199 4242

* ---- P8: under STRICTLY FINER clustering (cluster = unit, inside state)
fcarm fine unit 199 4242

dis as text "MC_FINECLUS_CTL_DONE"
