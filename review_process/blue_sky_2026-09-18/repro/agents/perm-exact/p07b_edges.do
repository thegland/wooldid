* p07_edges.do -- guards, boundary and membership-stress for permexact.
* Every case is wrapped in `capture noisily' and prints a sentinel, so batch Stata does not
* abort at the first legitimate refusal.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

capture program drop mkpanel
program define mkpanel
  syntax , n(integer) [PERcell(integer 3)]
  clear
  set obs `n'
  gen long unit = _n
  gen double ons = .
  expand 10
  bysort unit: gen int t = _n
  expand `percell'
end
capture program drop rpt
program define rpt
  args tag
  dis "RES|`tag'|rc=`=_rc'|NPOSS=" e(perm_npossible) "|REPS=" e(perm_reps) ///
      "|DONE=" e(perm_completed) "|FAIL=" e(perm_failed) "|EXACT=" e(perm_exact) ///
      "|INCOMPLETE=<`e(perm_incomplete)'>|ADDONE=<`e(perm_addone)'>"
  cap mat li e(perm_p)
end

*================================================================ E4/E5 permstrata pathologies
mkpanel, n(6)
replace ons = 4 if unit == 1
replace ons = 7 if unit == 2
gen byte stra = cond(unit <= 3, 1, 2)
replace stra = . if unit == 5
gen double y = 0.6*unit + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact permstrata(stra)
dis "E4_RC=" _rc
gen byte strb = cond(unit <= 3, 1, 2)
replace strb = 2 if unit == 1 & t > 5
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact permstrata(strb)
dis "E5_RC=" _rc

*================================================================ E6 over the 10000 ceiling
mkpanel, n(10) percell(1)
replace ons = 4 if inlist(unit, 1, 2)
replace ons = 6 if inlist(unit, 3, 4)
replace ons = 8 if inlist(unit, 5, 6)
gen double y = 0.6*unit + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact
dis "E6_RC=" _rc
dis "E6_perm_reps=" e(perm_reps) " E6_perm_exact=" e(perm_exact) " E6_cmd=`e(cmd)'"

*================================================================ E7/E8/E9 refusals
mkpanel, n(6)
replace ons = 4 if unit == 1
replace ons = 7 if unit == 2
gen double y = 0.6*unit + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
save p07_base, replace
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact permreps(100)
dis "E7_RC=" _rc
cap noi wooldid y unit t ons, cluster(unit) permute(probabilities) permexact permprobs((1 2 3 4 5 6) 4 .5 . .5)
dis "E8_RC=" _rc
gen byte onestr = unit
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact permstrata(onestr)
dis "E9_RC=" _rc

*================================================================ E10-E13 membership stress
use p07_base, clear
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact customvce(cluster unit)
rpt "E10_customvce"
use p07_base, clear
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact permtstat customvce(cluster unit)
rpt "E11_customvce_tstat"
use p07_base, clear
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact unconditionalse
rpt "E12_unconditionalse"
use p07_base, clear
replace y = exp(y/3)
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact poisson
rpt "E13_poisson"

*================================================================ E14 markout-dropped cluster
use p07_base, clear
replace y = . if unit == 6
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact
rpt "E14_phantomcluster"

*================================================================ E15 weights
use p07_base, clear
gen double w = 1 + mod(unit, 3)
cap noi wooldid y unit t ons [aw=w], cluster(unit) permute(shuffle) permexact
rpt "E15_aweight"

*================================================================ E16 permstepdown
use p07_base, clear
cap noi wooldid y unit t ons, cluster(unit) esprelength(2) espostlength(3) permute(shuffle) permexact permstepdown
rpt "E16_stepdown"
cap mat li e(perm_stepdown)

dis "ALLDONE"
