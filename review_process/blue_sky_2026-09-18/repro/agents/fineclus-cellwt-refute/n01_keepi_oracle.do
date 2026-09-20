*! n01 (NEW PROBE, not seeded) — does the NEW keep-i branch (T2/T3, 2026-09-16) deliver
*! an honest standard error for att_i when cluster() is genuinely FINER than i?
*! That is the configuration T2 was built to allow and the help now recommends
*! ("cluster BELOW the level of i").  Nothing screens it any more, so it must be right.
*!
*! Design: i = hospital SYSTEM (the cohort), cluster = hospital (5 per system), T = 6.
*! Treatment onset is a system-level property; tau_s is FIXED across the Monte Carlo,
*! so the target sampling distribution is "the cohorts in the data", exactly what the
*! keep-i branch claims to bootstrap.  Oracle = Monte-Carlo SD of att_i over 60 redraws
*! of the noise, with everything else held fixed.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkH
prog def mkH
  syntax , nseed(int)
  clear
  * fixed structure and fixed tau_s, set from a constant seed
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
  * only the NOISE is redrawn
  set seed `nseed'
  set sortseed `nseed'
  qui gen double hre = .
  qui bysort hosp (t): replace hre = rnormal(0,.8) if _n==1
  qui bysort hosp (t): replace hre = hre[1]
  gen double y = 2 + hre + .1*t + tau_s*D + rnormal(0,1)
  keep y sys hosp t onset
end

* ---------- oracle: MC sampling SD of att_i (tau_s FIXED) ----------
tempname mc
postfile `mc' draw double(b_att b_atti) using "n01_mc.dta", replace
forvalues d = 1/60 {
  mkH , nseed(`=20000+`d'')
  qui wooldid y sys t onset, att att_i robust
  mat bb = e(b)
  post `mc' (`d') (bb[1,1]) (bb[1,2])
}
postclose `mc'
preserve
use "n01_mc.dta", clear
summ b_atti
local mcsd_i = r(sd)
local mctru_i = r(mean)
summ b_att
local mcsd_a = r(sd)
restore
dis "MC_SD_ATTI=" `mcsd_i'
dis "MC_MEAN_ATTI=" `mctru_i'
dis "MC_SD_ATT=" `mcsd_a'

* ---------- bootstrap arms on 4 of those datasets ----------
tempname pf
postfile `pf' nseed str10 arm double(b_atti se_atti b_att se_att) using "n01_boot.dta", replace
foreach ns in 20003 20011 20027 {
  * keep-i regime: cluster(hosp) is strictly finer than i = sys
  mkH , nseed(`ns')
  cap noi wooldid y sys t onset, att att_i cluster(hosp) bootstrap(standard) reps(99) rseed(61) dots(0)
  dis "KEEPI_RC=" _rc " ns=`ns'"
  mat b1 = e(b)
  mat V1 = e(V)
  post `pf' (`ns') ("keepi") (b1[1,2]) (sqrt(V1[2,2])) (b1[1,1]) (sqrt(V1[1,1]))
  * nesting regime on the SAME data: cluster(sys) == i
  mkH , nseed(`ns')
  cap noi wooldid y sys t onset, att att_i cluster(sys) bootstrap(standard) reps(99) rseed(61) dots(0)
  dis "NEST_RC=" _rc " ns=`ns'"
  mat b2 = e(b)
  mat V2 = e(V)
  post `pf' (`ns') ("nest_eq_i") (b2[1,2]) (sqrt(V2[2,2])) (b2[1,1]) (sqrt(V2[1,1]))
  * analytic cluster(hosp): legitimate reference, since cluster is BELOW i
  mkH , nseed(`ns')
  cap noi wooldid y sys t onset, att att_i cluster(hosp)
  dis "ANAL_RC=" _rc " ns=`ns'"
  mat b3 = e(b)
  mat V3 = e(V)
  post `pf' (`ns') ("analytic") (b3[1,2]) (sqrt(V3[2,2])) (b3[1,1]) (sqrt(V3[1,1]))
}
postclose `pf'
use "n01_boot.dta", clear
format b_* se_* %12.0g
list, sep(3) noobs
bysort arm: summ se_atti se_att
dis "MC_SD_ATTI was " `mcsd_i' "   MC_SD_ATT was " `mcsd_a'
dis "predicted cross-system composition SE if it leaked in = sqrt(Var_s/8) = " sqrt(.1125/8)
