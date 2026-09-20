* p7_nest_seeds.do -- is the p6 att_i standard-error collapse systematic, or one seed?
* Same design, four bootstrap seeds, reps(400).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootsample"

clear
set obs 16
gen int state = _n
gen double sfx = rnormal()*0.5
gen int onset = .
replace onset = 2003 + mod(_n-1,4) if _n <= 8
expand 3
bysort state: gen int unit = (state-1)*3 + _n
gen double ufx = rnormal()*0.4
expand 10
bysort unit: gen int year = 1999 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .10*(year-1999) + sfx + ufx + (.3 + .4*sfx)*post + rnormal()*0.6
gen byte stray = 0
set obs `=_N+1'
replace unit  = 19   in l
replace state = 8    in l
replace year  = 2009 in l
replace onset = 2005 in l
replace y     = .    in l
replace stray = 1    in l
tempfile full
save `full', replace

foreach s in 101 202 303 404 {
  use `full', clear
  cap noi qui wooldid y unit year onset, cluster(state) att_i bootstrap(standard) reps(400) rseed(`s') dots(0)
  mat vf = e(V)
  local sef = sqrt(vf[1,1])
  local nf  = e(N)
  use `full', clear
  drop if stray
  cap noi qui wooldid y unit year onset, cluster(state) att_i bootstrap(standard) reps(400) rseed(`s') dots(0)
  mat vo = e(V)
  local seo = sqrt(vo[1,1])
  local no  = e(N)
  dis "SEED`s'_NFULL=" `nf' " NORACLE=" `no' " SE_FULL=" `sef' " SE_ORACLE=" `seo' " RATIO=" `sef'/`seo'
}
* negative control: give the stray row a real outcome, so the straddle is REAL and
* the keep-i regime is the correct one on both sides
use `full', clear
replace y = 1.5 if stray
cap noi qui wooldid y unit year onset, cluster(state) att_i bootstrap(standard) reps(400) rseed(101) dots(0)
mat vc = e(V)
dis "CTRL_REALSTRADDLE_N=" e(N) " SE=" sqrt(vc[1,1])
dis "P7_ALLDONE"
