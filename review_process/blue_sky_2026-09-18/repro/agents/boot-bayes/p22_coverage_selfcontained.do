*=============================================================================
* p22  SELF-CONTAINED REPRO: COVERAGE of bootstrap(bayes) intervals at
*      cluster() == i, nominal 95%, against bootstrap(standard) and the analytic
*      clustered interval on the SAME 30 simulated panels.
*      (identical to p12b_coverage.do + p20_cov_analyze.do, with the str8 -> str16
*       fix that stopped "main_att_i" truncating onto "main_att")
*=============================================================================
* COVERAGE of bootstrap(bayes) intervals at cluster() == i, nominal 95%.
* Context: the project has an OPEN, unruled finding that bootstrap(STANDARD) at
* cluster()==i OVER-covers att_i (149/150 = .9933).  The bayes side of that cell -
* and especially bayes + a cell-weighted estimand under safety(off), the combination
* the package refuses by default - has never been measured.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

local NDRAW 30
local REPS  79
local TRUE  1.0

tempname mem
postfile `mem' int draw str10 arm str16 est double b double se double lo_p double hi_p double lo_w double hi_w using p22_res, replace

forval d = 1/`NDRAW' {
  clear
  set seed `=90000 + `d''
  set sortseed `=90000 + `d''
  qui set obs 30
  qui gen int st = _n
  qui gen double a_i = rnormal()
  qui gen int onset = .
  qui replace onset = 2005 if st <= 8
  qui replace onset = 2008 if st > 8 & st <= 16
  qui replace onset = 2011 if st > 16 & st <= 22
  qui expand 10
  qui bysort st: gen int year = 2001 + _n
  qui gen byte post = (year >= onset) & !missing(onset)
  * homogeneous effect => att and att_i share the same truth
  qui gen double y = a_i + 0.03*(year-2001) + `TRUE'*post + rnormal()*0.7
  tempfile dd
  qui save `dd'

  * arm 1: bayes, att AND att_i (the safety(off) combination)
  use `dd', clear
  cap qui wooldid y st year onset, att att_i cluster(st) bootstrap(bayes) reps(`REPS') rseed(`=7000+`d'') dots(0) safety(off)
  if _rc == 0 {
    matrix BB = e(b)
    matrix VV = e(V)
    matrix CP = e(boot_ci_percentile)
    matrix CN = e(boot_ci_normal)
    local cn : colnames BB
    local k 0
    foreach nm of local cn {
      local ++k
      post `mem' (`d') ("bayes") ("`nm'") (BB[1,`k']) (sqrt(VV[`k',`k'])) (CP[1,`k']) (CP[2,`k']) (CN[1,`k']) (CN[2,`k'])
    }
  }
  * arm 2: standard, att only (att_i there takes the physical re-id leg, ~2x cost)
  use `dd', clear
  cap qui wooldid y st year onset, att cluster(st) bootstrap(standard) reps(`REPS') rseed(`=7000+`d'') dots(0)
  if _rc == 0 {
    matrix BB = e(b)
    matrix VV = e(V)
    matrix CP = e(boot_ci_percentile)
    matrix CN = e(boot_ci_normal)
    post `mem' (`d') ("standard") ("main_att") (BB[1,1]) (sqrt(VV[1,1])) (CP[1,1]) (CP[2,1]) (CN[1,1]) (CN[2,1])
  }
  * arm 3: analytic clustered, for reference
  use `dd', clear
  cap qui wooldid y st year onset, att att_i cluster(st)
  if _rc == 0 {
    matrix MR = e(wdidmainresults)
    local rn : rownames MR
    local k 0
    foreach nm of local rn {
      local ++k
      post `mem' (`d') ("analytic") ("`nm'") (MR[`k',1]) (MR[`k',2]) (.) (.) (MR[`k',6]) (MR[`k',7])
    }
  }
  if mod(`d',10)==0 dis "DRAWDONE `d'"
}
postclose `mem'

*--- analysis ---------------------------------------------------------------
use p22_res, clear
gen long ord = _n
bysort arm draw (ord): gen byte k = _n
gen str12 est2 = cond(k==1, "att", "att_i")
replace est2 = "att" if arm=="standard"
gen byte covp = (lo_p <= `TRUE') & (hi_p >= `TRUE') if !missing(lo_p)
gen byte covw = (lo_w <= `TRUE') & (hi_w >= `TRUE') if !missing(lo_w)
* the analytic arm posted b and se; build its normal interval here
replace covw = (b - 1.959964*se <= `TRUE') & (b + 1.959964*se >= `TRUE') if arm=="analytic"
gen double widp = hi_p - lo_p
collapse (sum) ncovp=covp ncovw=covw (mean) mse=se mwidp=widp (sd) sdb=b (count) n=b, by(arm est2)
gen double ratio = mse/sdb
list arm est2 n ncovp ncovw mse sdb ratio mwidp, noobs sepby(arm) abbrev(12)
dis "ncovp / ncovw are COUNTS covered out of n at nominal 95% (percentile / wald)"
dis "ratio = mean REPORTED se divided by the TRUE sampling sd of the estimate"
dis "PROBE_DONE"
