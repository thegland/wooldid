*! r01a — REFUTE FCW-01 on a REALISTIC panel.
*! Question: does ONE straddling (unit,period) cell really move se(att_i) ~23%,
*! or is 23% an artifact of a near-degenerate DGP (tau = f(cluster) exactly, noise sd .02)?
*! Self-contained: generates its own county x firm panel, no project data, no project suite.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , dseed(int) [taucsd(real .4) tauisd(real .4) esd(real 1)]
  local G 15
  local U 4
  local T 7
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs `G'
  gen long cty = _n
  gen double a_g = rnormal(0,`taucsd')
  gen onset = .
  replace onset = 4 if cty <= 5
  replace onset = 6 if cty > 5 & cty <= 10
  expand `U'
  sort cty
  gen long unit = _n
  gen double b_i = rnormal(0,`tauisd')
  gen double uf  = rnormal(0,1)
  expand `T'
  bysort unit (cty): gen int t = _n
  gen double tau = 1 + a_g + b_i
  gen byte D = (onset < . & t >= onset)
  gen double y = 2 + uf + .1*t + tau*D + rnormal(0,`esd')
  keep y unit cty t onset
end

tempname pf
postfile `pf' dseed rseed str6 arm double(b_att b_atti se_att se_atti) nunits nobs using "r01a.dta", replace

foreach ds in 4101 4102 4103 {
 foreach rs in 31 32 {
  * ---- arm A : cluster(cty) nests unit (every unit wholly inside one county)
  mkpanel , dseed(`ds')
  qui gegen double cmn = min(cty), by(unit)
  qui gegen double cmx = max(cty), by(unit)
  qui count if cmn != cmx
  dis "ARM_A dseed=`ds' straddlers=" r(N)
  drop cmn cmx
  cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(199) rseed(`rs') dots(0)
  dis "A_RC=" _rc
  mat bA = e(b)
  mat VA = e(V)
  post `pf' (`ds') (`rs') ("A") (bA[1,1]) (bA[1,2]) (sqrt(VA[1,1])) (sqrt(VA[2,2])) (60) (e(N))

  * ---- arm B : ONE cell of the cluster variable moved. unit 1 spends its last
  *      period in county 2 (a firm reclassified mid-panel). Nothing else changes.
  mkpanel , dseed(`ds')
  qui replace cty = 2 if unit == 1 & t == 7
  qui gegen double cmn = min(cty), by(unit)
  qui gegen double cmx = max(cty), by(unit)
  qui count if cmn != cmx
  dis "ARM_B dseed=`ds' straddlers=" r(N)
  drop cmn cmx
  cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(199) rseed(`rs') dots(0)
  dis "B_RC=" _rc
  mat bB = e(b)
  mat VB = e(V)
  post `pf' (`ds') (`rs') ("B") (bB[1,1]) (bB[1,2]) (sqrt(VB[1,1])) (sqrt(VB[2,2])) (60) (e(N))
 }
}
postclose `pf'

use "r01a.dta", clear
list, sep(0) noobs
reshape wide b_att b_atti se_att se_atti nunits nobs, i(dseed rseed) j(arm) string
gen double ratio_atti = se_attiB/se_attiA
gen double ratio_att  = se_attB/se_attA
gen double bdiff      = b_attiB - b_attiA
format ratio* %9.4f
list dseed rseed ratio_atti ratio_att bdiff, sep(0) noobs
summ ratio_atti ratio_att
dis "MEAN_RATIO_ATTI=" r(mean)
