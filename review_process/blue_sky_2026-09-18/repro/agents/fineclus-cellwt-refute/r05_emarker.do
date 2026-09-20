*! r05 — REFUTE FCW-05: is the bootstrap REGIME really unrecorded?  Instead of grepping,
*! run both regimes on one dataset with the same seed and DIFF every e() name and value,
*! plus the full printed output.  If any e() or any printed line differs other than the
*! numbers themselves, a user CAN tell which rule ran and FCW-05 is refuted.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 8801
set sortseed 8801
set obs 15
gen long cty = _n
gen double a_g = .3*cty
gen onset = .
replace onset = 4 if cty <= 5
replace onset = 6 if cty > 5 & cty <= 10
expand 4
sort cty
gen long unit = _n
expand 7
bysort unit (cty): gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double y = 2 + .1*t + (1+a_g)*D + rnormal(0,.05)
keep y unit cty t onset
save r05dat.dta, replace

cap prog drop dumpe
prog def dumpe
  syntax , tag(string)
  file open fh using "r05_`tag'.txt", write replace
  foreach nm in `:e(scalars)' {
    file write fh "S `nm' = " %21x (e(`nm')) _n
  }
  foreach nm in `:e(macros)' {
    file write fh "M `nm' = " `"`e(`nm')'"' _n
  }
  foreach nm in `:e(matrices)' {
    file write fh "X `nm'" _n
  }
  file close fh
end

use r05dat.dta, clear
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(99) rseed(71) dots(0)
dis "A_RC=" _rc
dumpe , tag(A)
mat VA = e(V)
dis "A se_atti = " %21x sqrt(VA[2,2])

use r05dat.dta, clear
replace cty = 2 if unit == 1 & t == 7
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(99) rseed(71) dots(0)
dis "B_RC=" _rc
dumpe , tag(B)
mat VB = e(V)
dis "B se_atti = " %21x sqrt(VB[2,2])
dis "RATIO = " sqrt(VB[2,2])/sqrt(VA[2,2])
