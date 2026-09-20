* probe10 -- permstrata(): is the randomization design the one the help describes, and is
* e(perm_npossible) the right count?  Checked structurally off the permsave dump.
* (a) sampled: 24 units in 4 strata of 6, each stratum 2/2/2 over {2004, 2006, never}
*     -> npossible must be (6!/(2!2!2!))^4 = 90^4 = 65,610,000, and EVERY replication must
*     preserve each stratum's onset composition.
* (b) permexact + permstrata: 8 units in 2 strata of 4, each 2 treated / 2 never
*     -> M must be C(4,2)^2 = 36 and the enumeration must be exactly the 36 distinct
*     composition-preserving assignments.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"

* ---------------- (a) ----------------
clear
set obs 24
gen long unit = _n
gen long strat = ceil(unit/6)
gen int pos = unit - 6*(strat-1)
gen double ttre = .
replace ttre = 2004 if inlist(pos,1,2)
replace ttre = 2006 if inlist(pos,3,4)
gen double ai = rnormal()
keep unit strat ttre ai
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ai + 0.2*(t-2000) + 1.3*post + 0.7*rnormal()
tab strat ttre, m
tempfile basea
save `basea'
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) permute(shuffle) permstrata(strat) permreps(99) rseed(23) permsave(st_a) replace
dis "RC_STRA=" _rc
dis "npossible = " e(perm_npossible) "  (expected 65610000)"
local htag STRATA
local habsp 0
local hexact 0
include "harvest.do"
preserve
use st_a.dta, clear
qui ds onset_*
local ons `r(varlist)'
gen byte bad = 0
* stratum 1 = units 1-6, stratum 2 = 7-12, stratum 3 = 13-18, stratum 4 = 19-24
forvalues s = 1/4 {
  local lo = 6*(`s'-1)+1
  local hi = 6*`s'
  local l04
  local l06
  local lnv
  forvalues u = `lo'/`hi' {
    local l04 `l04' + (onset_`u' == 2004)
    local l06 `l06' + (onset_`u' == 2006)
    local lnv `lnv' + missing(onset_`u')
  }
  gen byte n04_`s' = 0 `l04'
  gen byte n06_`s' = 0 `l06'
  gen byte nnv_`s' = 0 `lnv'
  qui replace bad = 1 if n04_`s' != 2 | n06_`s' != 2 | nnv_`s' != 2
}
count if bad
dis "STRATA composition violations = " r(N) " of " _N " rows"
count
dis "STRATA rows in dump = " r(N)
restore

* ---------------- (b) ----------------
clear
set obs 8
gen long unit = _n
gen long strat = ceil(unit/4)
gen int pos = unit - 4*(strat-1)
gen double ttre = .
replace ttre = 2004 if inlist(pos,1,2)
gen double ai = rnormal()
keep unit strat ttre ai
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ai + 0.2*(t-2000) + 2.0*post + 0.5*rnormal()
tab strat ttre, m
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) permute(shuffle) permstrata(strat) permexact permsave(st_b) replace
dis "RC_STRB=" _rc
dis "npossible = " e(perm_npossible) "  (expected 36)"
local htag STRATAX
local habsp 0
local hexact 1
include "harvest.do"
matrix list e(perm_p), format(%12.8f)
preserve
use st_b.dta, clear
dis "rows in exact dump = " _N
qui ds onset_*
local ons `r(varlist)'
egen str40 key = concat(`ons'), punct("|")
qui drop if _rep == 0
duplicates report key
bysort key: gen byte dup = _N
count if dup > 1
dis "EXACT duplicate assignments = " r(N)
gen byte bad = 0
forvalues s = 1/2 {
  local lo = 4*(`s'-1)+1
  local hi = 4*`s'
  local l04
  forvalues u = `lo'/`hi' {
    local l04 `l04' + (onset_`u' == 2004)
  }
  gen byte n04_`s' = 0 `l04'
  qui replace bad = 1 if n04_`s' != 2
}
count if bad
dis "EXACT composition violations = " r(N)
restore
dis "PROBE10_DONE"
