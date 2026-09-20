*=============================================================================
* p14b  SELF-CONTAINED REPRO (needs only the pinned build; rwgen ships with Stata 19).
*
* CLAIM: bootstrap(bayes) has an EXTERNAL ORACLE after all - Stata's own
* rwgen/bayesboot - and wooldid matches it BITWISE at a common seed:
*   bootalpha(size)  == rwgen bayes ..., id(cluster)
*   bootalpha(equal) == rwgen bayes ..., id(cluster) priorpowers(-n_c), times n_c
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
which wooldid
dis "Stata version = " c(stata_version)
local R 20
local SEED 31415
set seed 20260916
set sortseed 20260916
clear
set obs 12
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2005 if st <= 3
replace onset = 2008 if st > 3 & st <= 6
replace onset = 2011 if st > 6 & st <= 8
expand 14
bysort st: gen int year = 2000 + _n
bysort st: drop if _n > (6 + st/2)
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.03*(year-2000) + 1.0*post + rnormal()*0.5
tempfile dd
qui save `dd'
foreach am in size equal {
  use `dd', clear
  qui wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(`am') reps(`R') rseed(`SEED') bootsave(p14b_`am') replace dots(0)
  use p14b_`am', clear
  qui keep if _rep >= 1
  qui sort _rep
  qui ds w_*
  mkmat `r(varlist)', matrix(WD`am')
}
use `dd', clear
qui bysort st (year): gen long nc = _N
qui gen double negnc = -nc
set seed `SEED'
qui rwgen bayes RS, reps(`R') id(st)
set seed `SEED'
qui rwgen bayes RE, reps(`R') id(st) priorpowers(negnc)
qui bysort st (year): keep if _n == 1
qui sort st
mkmat RS1-RS`R', matrix(RS)
mkmat RE1-RE`R', matrix(RE)
mkmat nc, matrix(NCv)
mata {
  WS = st_matrix("WDsize"); WE = st_matrix("WDequal")
  S = st_matrix("RS")'; E = st_matrix("RE")'; nn = st_matrix("NCv")
  printf("bootalpha(size)  bitwise == rwgen bayes id(st)            : %f   max|diff| = %21x\n", (WS==S), max(abs(WS-S)))
  printf("bootalpha(equal) bitwise == rwgen priorpowers(-n_c) * n_c : %f   max rel diff = %g\n", (WE==(E:*nn')), max(abs((WE-E:*nn'):/(E:*nn'))))
}
dis "PROBE_DONE"
