* p12: coverage of bootstrap(bayes) intervals at cluster() == i, the configuration
* where the project has an OPEN unruled finding that bootstrap(standard) OVER-covers
* (att_i 149/150 = .9933, recon known_bydesign.md sec 7). Nominal 95%.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

local NDRAW 100
local REPS  99
local TRUE  1.0

tempname mem
postfile `mem' int draw str8 mode double b double se double lo_p double hi_p double lo_w double hi_w int fail using p12_res, replace

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
  qui gen double y = a_i + 0.03*(year-2001) + `TRUE'*post + rnormal()*0.7
  tempfile dd
  qui save `dd'
  foreach md in bayes standard {
    use `dd', clear
    cap qui wooldid y st year onset, att cluster(st) bootstrap(`md') reps(`REPS') rseed(`=7000+`d'') dots(0)
    if _rc == 0 {
      matrix BB = e(b)
      matrix VV = e(V)
      matrix CP = e(boot_ci_percentile)
      matrix CN = e(boot_ci_normal)
      post `mem' (`d') ("`md'") (BB[1,1]) (sqrt(VV[1,1])) (CP[1,1]) (CP[2,1]) (CN[1,1]) (CN[2,1]) (`=e(boot_failed)')
    }
  }
  if mod(`d',10)==0 dis "draw `d' done"
}
postclose `mem'

use p12_res, clear
gen byte covp = (lo_p <= `TRUE') & (hi_p >= `TRUE')
gen byte covw = (lo_w <= `TRUE') & (hi_w >= `TRUE')
gen double wid_p = hi_p - lo_p
gen double wid_w = hi_w - lo_w
bysort mode: su covp covw wid_p wid_w se b fail
table mode, stat(mean covp covw wid_w se) stat(sd b) stat(n b)
dis "PROBE_DONE"
