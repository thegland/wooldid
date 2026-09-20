*=========================================================================
* repro_S3_esprelength.do -- SELF-CONTAINED, ~1 minute.
*
* CLAIM DEMONSTRATED: on ONE unchanged dataset, adding esprelength() changes
* the reported headline att and widens its standard error, with no message.
* espostlength() on its own does NOT: the mover is esprelength() alone,
* because the reference-normalization pool is {relyr < -esprelength} plus
* each treated unit's forced first in-sample observation, so requested PRE
* horizons are taken out of the reference pool.  e(escomp_refN) records the
* shrinkage (270 -> 180 here) but is only described under Stored Results;
* the esprelength() option entry does not mention the effect.
*
* e(N) and e(df_r) are identical across all runs, which rules out any change
* of sample or degrees of freedom.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_S3_esprelength.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

set seed 4411
set sortseed 4411
clear
set obs 60
gen long unit = _n
gen byte cg = mod(_n-1,4)+1
gen double ttre = .
replace ttre = 2004 if cg==1
replace ttre = 2006 if cg==2
replace ttre = 2008 if cg==3
gen double ui   = rnormal(0,0.70)
gen double hmul = 1 + rnormal(0,0.25)
expand 12
bys unit: gen int t = 2000 + _n - 1
gen double e = .
bys unit (t): replace e = rnormal(0,0.55) if _n==1
bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
gen double truetau = 0
replace truetau = 1.00*hmul if cg==1 & t>=ttre
replace truetau = 1.15*hmul if cg==2 & t>=ttre
replace truetau = 0.85*hmul if cg==3 & t>=ttre
gen double y = ui + 0.25*(t-2000) + truetau + e
tempfile d
qui save "`d'"

di as txt _n "spec                                b            se        N   df_r   refN"
foreach spec in none post2 pre2 both2 both3 {
  use "`d'", clear
  local opt ""
  if "`spec'"=="post2" local opt "espostlength(2)"
  if "`spec'"=="pre2"  local opt "esprelength(2) espostlength(-999)"
  if "`spec'"=="both2" local opt "esprelength(2) espostlength(2)"
  if "`spec'"=="both3" local opt "esprelength(3) espostlength(3)"
  qui wooldid y cg t ttre, att cluster(unit) `opt'
  mat T = e(wdidmainresults)
  local rr = .
  cap local rr = e(escomp_refN)
  di as txt %-12s "`spec'" "  " %14.9f T[1,1] "  " %12.9f T[1,2] "  " %5.0f e(N) "  " %5.0f e(df_r) "   `rr'"
}
di as txt _n "Now the same runs NOISILY, so that any Caution would be visible:"
foreach spec in none post2 pre2 both3 {
  use "`d'", clear
  local opt ""
  if "`spec'"=="post2" local opt "espostlength(2)"
  if "`spec'"=="pre2"  local opt "esprelength(2) espostlength(-999)"
  if "`spec'"=="both3" local opt "esprelength(3) espostlength(3)"
  di as txt _n(2) "=========== spec = `spec' : `opt' ==========="
  wooldid y cg t ttre, att cluster(unit) `opt'
}
di as txt _n "ALLDONE_S3"
