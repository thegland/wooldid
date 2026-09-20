*! r03 — REFUTE FCW-03: is the bootstrap SE of att_i at cluster()==i dominated by an
*! ADDED cohort-composition term Var_s(tau_s)/S?  Decisive cheap test: identical noise
*! draws, identical truth for att_i, only the SPREAD of tau across treated units differs.
*! HET: tau in {.4,1.2} (Var = .16).  HOM: tau = .8 for every treated unit (Var = 0).
*! If the composition mechanism is right, se(att_i) must be ~sqrt(.16/S) under HET and
*! ~0 under HOM even though the two runs see the same noise.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkE
prog def mkE
  syntax , dseed(int) het(int) [esd(real .02)]
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs 40
  gen long state = _n
  gen onset = cond(mod(state,2)==1, 4, .)
  gen double tau = cond(`het'==1, cond(mod(state,4)==1, .4, 1.2), .8)
  replace tau = 0 if missing(onset)
  gen double sf = rnormal(0,1)
  expand 6
  bysort state (onset): gen int t = _n
  gen byte D = (onset < . & t >= onset)
  gen double y = 2 + sf + .15*t + tau*D + rnormal(0,`esd')
  keep y state t onset
end

tempname pf
postfile `pf' het dseed rseed str8 clusvar double(b_att b_atti se_att se_atti) using "r03.dta", replace

foreach ds in 9001 9002 {
 foreach h in 1 0 {
  foreach rs in 51 52 {
   * cluster() == i : the NESTING regime, __wdbnewi = group(drawid, i)
   mkE , dseed(`ds') het(`h')
   cap noi wooldid y state t onset, att att_i cluster(state) bootstrap(standard) reps(199) rseed(`rs') dots(0)
   dis "EQ_RC=" _rc " het=`h' ds=`ds' rs=`rs'"
   mat b1 = e(b)
   mat V1 = e(V)
   post `pf' (`h') (`ds') (`rs') ("eq_i") (b1[1,1]) (b1[1,2]) (sqrt(V1[1,1])) (sqrt(V1[2,2]))
  }
 }
}
postclose `pf'
use "r03.dta", clear
format b_* se_* %16.0g
list, sep(0) noobs
dis "predicted composition SE under HET = sqrt(.16/20) = " sqrt(.16/20)
