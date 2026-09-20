* p11 (my own angle): the help sells bayes on "no cluster ever leaves the sample,
* and replication estimability is stable" (sthlp:532). With FEW never-treated
* clusters, bootstrap(standard) FAILS the replications that drew none - and SAYS
* SO ("n of reps failed"). bayes cannot fail: a control cluster that draws a
* near-zero Dirichlet mass is still "in the sample" but contributes nothing, and
* the replicate estimate blows up SILENTLY. Measure what that does to e(V).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

clear
set obs 20
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2005 if st <= 6
replace onset = 2008 if st > 6 & st <= 12
replace onset = 2011 if st > 12 & st <= 18
* st 19, 20 never treated  (2 control clusters out of 20 - common in state policy work)
expand 12
bysort st: gen int year = 2000 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.03*(year-2000) + 1.0*post + rnormal()*0.5
save p11_panel, replace

local R 500
foreach md in standard bayes {
  use p11_panel, clear
  cap noi wooldid y st year onset, att cluster(st) bootstrap(`md') reps(`R') rseed(2468) dots(0) bootsave(p11_`md') replace
  dis "`md'_RC=" _rc
  dis "`md' completed=" e(boot_completed) " failed=" e(boot_failed)
  matrix BB = e(b)
  matrix VV = e(V)
  matrix BD = e(bootdist)
  dis "`md' b=" BB[1,1] "  se(e(V))=" sqrt(VV[1,1])
  matrix CN = e(boot_ci_normal)
  matrix CP = e(boot_ci_percentile)
  dis "`md' wald CI = [" CN[1,1] ", " CN[2,1] "]  width=" CN[2,1]-CN[1,1]
  dis "`md' pctl CI = [" CP[1,1] ", " CP[2,1] "]  width=" CP[2,1]-CP[1,1]
  preserve
    clear
    qui svmat double BD
    qui drop if missing(BD1)
    qui su BD1, detail
    dis "`md' draws: n=" r(N) " sd=" r(sd) " skew=" r(skewness) " kurt=" r(kurtosis)
    dis "`md' draws: min=" r(min) " p1=" r(p1) " p50=" r(p50) " p99=" r(p99) " max=" r(max)
    * robust scale for comparison
    qui su BD1, detail
    local iqr = r(p75)-r(p25)
    dis "`md' IQR-based sd = " `iqr'/1.349 "   ratio sd/IQRsd = " r(sd)/(`iqr'/1.349)
  restore
}

*--- mechanism: link the control clusters' Dirichlet mass to the replicate value
use p11_bayes, clear
qui keep _rep w_19 w_20 main_att
qui keep if _rep >= 1
gen double ctrlmass = w_19 + w_20
qui su ctrlmass, detail
dis "control-cluster Dirichlet mass: min=" r(min) " p1=" r(p1) " p50=" r(p50)
qui count if ctrlmass < 0.01
dis "reps with control mass < 0.01 : " r(N)
qui su main_att
local mb = r(mean)
gen double dev = abs(main_att - `mb')
corr dev ctrlmass
dis "corr(|deviation|, control mass) = " r(rho)
qui su dev if ctrlmass < 0.02
dis "mean |dev| when control mass<0.02 : " r(mean) " (n=" r(N) ")"
qui su dev if ctrlmass >= 0.02
dis "mean |dev| otherwise              : " r(mean) " (n=" r(N) ")"
gsort -dev
list _rep main_att ctrlmass in 1/8
dis "PROBE_DONE"
