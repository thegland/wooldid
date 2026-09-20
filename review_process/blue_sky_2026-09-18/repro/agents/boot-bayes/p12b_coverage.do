* p12b: COVERAGE of bootstrap(bayes) intervals at cluster() == i, nominal 95%.
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
postfile `mem' int draw str10 arm str8 est double b double se double lo_p double hi_p double lo_w double hi_w using p12_res, replace

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

use p12_res, clear
gen byte covp = (lo_p <= `TRUE') & (hi_p >= `TRUE') if !missing(lo_p)
gen byte covw = (lo_w <= `TRUE') & (hi_w >= `TRUE') if !missing(lo_w)
gen double widw = hi_w - lo_w
bysort arm est: gen long nn = _N
collapse (mean) covp covw widw se (sd) sdb=b (count) n=b, by(arm est)
list, sepby(arm) noobs
dis "COVERAGE: covp = percentile family, covw = wald/normal family (analytic arm: analytic CI)"
dis "se vs sdb: se is the MEAN reported standard error, sdb the true sampling sd of b"
dis "PROBE_DONE"
