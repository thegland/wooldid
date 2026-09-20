*! n02 (NEW PROBE, not seeded) — the THIRD refusal T2 deleted: bootstrap(standard) with a
*! cell-weighted estimand and NO cluster().  Old ado:6353 refused this r198.  It now runs,
*! resampling OBSERVATIONS.  Question: on a panel with ordinary within-unit serial
*! dependence (a unit random effect), is the resulting se(att_i) honest?
*! Oracle: Monte-Carlo SD of att_i over 60 redraws of the noise with the design fixed.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkH
prog def mkH
  syntax , nseed(int)
  clear
  set seed 12345
  set sortseed 12345
  set obs 12
  gen long sys = _n
  gen onset = cond(sys <= 4, 4, cond(sys <= 8, 5, .))
  gen double tau_s = cond(missing(onset), 0, .4 + .3*mod(sys,4))
  expand 5
  sort sys
  gen long hosp = _n
  expand 6
  bysort hosp (sys): gen int t = _n
  gen byte D = (onset < . & t >= onset)
  set seed `nseed'
  set sortseed `nseed'
  qui gen double hre = .
  qui bysort hosp (t): replace hre = rnormal(0,.8) if _n==1
  qui bysort hosp (t): replace hre = hre[1]
  gen double y = 2 + hre + .1*t + tau_s*D + rnormal(0,1)
  keep y sys hosp t onset
end

tempname mc
postfile `mc' draw double(b_att b_atti) using "n02_mc.dta", replace
forvalues d = 1/60 {
  mkH , nseed(`=20000+`d'')
  qui wooldid y sys t onset, att att_i robust
  mat bb = e(b)
  post `mc' (`d') (bb[1,1]) (bb[1,2])
}
postclose `mc'
preserve
use "n02_mc.dta", clear
summ b_atti
local mcsd_i = r(sd)
summ b_att
local mcsd_a = r(sd)
restore
dis "MC_SD_ATTI=" `mcsd_i'
dis "MC_SD_ATT=" `mcsd_a'

tempname pf
postfile `pf' nseed str10 arm double(se_atti se_att) using "n02_boot.dta", replace
foreach ns in 20003 20011 20027 {
  mkH , nseed(`ns')
  cap noi wooldid y sys t onset, att att_i bootstrap(standard) reps(99) rseed(61) dots(0)
  dis "NOCLUS_RC=" _rc " ns=`ns'"
  mat V1 = e(V)
  post `pf' (`ns') ("noclus") (sqrt(V1[2,2])) (sqrt(V1[1,1]))
  mkH , nseed(`ns')
  cap noi wooldid y sys t onset, att att_i cluster(hosp) bootstrap(standard) reps(99) rseed(61) dots(0)
  dis "CLHOSP_RC=" _rc " ns=`ns'"
  mat V2 = e(V)
  post `pf' (`ns') ("cl_hosp") (sqrt(V2[2,2])) (sqrt(V2[1,1]))
}
postclose `pf'
use "n02_boot.dta", clear
format se_* %12.0g
list, sep(2) noobs
bysort arm: summ se_atti se_att
dis "MC_SD_ATTI was " `mcsd_i' "  MC_SD_ATT was " `mcsd_a'
