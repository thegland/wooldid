*! r02 — REFUTE FCW-02: can a row the ESTIMATOR DISCARDS (missing outcome) change the
*! bootstrap standard error?  Arms differ ONLY in the cluster code carried by a row
*! whose outcome is missing.  e(N), e(sample) and the point estimates must be identical.
*! Self-contained; two DGPs (realistic, and near-degenerate pure-composition).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , dseed(int) dgp(string)
  local G 15
  local U 4
  local T 7
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs `G'
  gen long cty = _n
  if "`dgp'" == "real"  gen double a_g = rnormal(0,.4)
  if "`dgp'" == "degen" gen double a_g = .3*cty
  gen onset = .
  replace onset = 4 if cty <= 5
  replace onset = 6 if cty > 5 & cty <= 10
  expand `U'
  sort cty
  gen long unit = _n
  if "`dgp'" == "real"  gen double b_i = rnormal(0,.4)
  if "`dgp'" == "degen" gen double b_i = 0
  gen double uf = cond("`dgp'"=="real", rnormal(0,1), 0)
  expand `T'
  bysort unit (cty): gen int t = _n
  gen double tau = 1 + a_g + b_i
  gen byte D = (onset < . & t >= onset)
  local esd = cond("`dgp'"=="real", 1, .02)
  gen double y = 2 + uf + .1*t + tau*D + rnormal(0,`esd')
  * the PHANTOM row: unit 1 at t = 3 has a missing outcome in BOTH arms.
  replace y = . if unit == 1 & t == 3
  keep y unit cty t onset
end

tempname pf
postfile `pf' str6 dgp dseed rseed str6 arm double(b_att b_atti se_att se_atti) nobs using "r02.dta", replace

foreach g in real {
 foreach ds in 7201 7202 {
  foreach rs in 41 42 {
   * arm A: cluster variable untouched. unit 1 wholly inside county 1.
   mkpanel , dseed(`ds') dgp("`g'")
   cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(199) rseed(`rs') dots(0)
   dis "A_RC=" _rc " dgp=`g' ds=`ds' rs=`rs' N=" e(N)
   mat bA = e(b)
   mat VA = e(V)
   post `pf' ("`g'") (`ds') (`rs') ("A") (bA[1,1]) (bA[1,2]) (sqrt(VA[1,1])) (sqrt(VA[2,2])) (e(N))

   * arm D: the SAME missing-outcome row is recorded under county 2. No estimable
   * row moves; the estimation sample is byte-identical to arm A.
   mkpanel , dseed(`ds') dgp("`g'")
   qui replace cty = 2 if unit == 1 & t == 3
   cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(199) rseed(`rs') dots(0)
   dis "D_RC=" _rc " dgp=`g' ds=`ds' rs=`rs' N=" e(N)
   mat bD = e(b)
   mat VD = e(V)
   post `pf' ("`g'") (`ds') (`rs') ("D") (bD[1,1]) (bD[1,2]) (sqrt(VD[1,1])) (sqrt(VD[2,2])) (e(N))
  }
 }
}
postclose `pf'
use "r02.dta", clear
format b_* se_* %16.0g
list, sep(0) noobs
reshape wide b_att b_atti se_att se_atti nobs, i(dgp dseed rseed) j(arm) string
gen double r_atti = se_attiD/se_attiA
gen double r_att  = se_attD/se_attA
gen double dbi    = b_attiD - b_attiA
gen double dN     = nobsD - nobsA
format r_* %9.4f
list dgp dseed rseed r_atti r_att dbi dN, sep(0) noobs
