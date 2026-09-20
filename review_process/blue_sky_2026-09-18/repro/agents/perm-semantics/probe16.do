* probe16 -- IS THE SAMPLED SHUFFLE UNIFORM over the assignment set?  Everything the
* p-value's exactness rests on assumes the sampled draws are uniform over the M distinct
* label permutations.  8 clusters, 2 treated -> M = 28.  Draw 1400 replications and test
* the empirical assignment frequencies against 50 each.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
clear
set obs 8
gen long unit = _n
gen double ttre = cond(inlist(unit,1,2), 2004, .)
gen double ai = rnormal()
expand 6
bysort unit: gen int t = 2002 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ai + 0.2*(t-2002) + 1.0*post + 0.5*rnormal()
cap noi wooldid y unit t ttre, cluster(unit) att permute(shuffle) permreps(1400) rseed(97) permsave(u_shuf) replace
dis "RC_U=" _rc
dis "completed=" e(perm_completed) " failed=" e(perm_failed) " npossible=" e(perm_npossible)
use u_shuf.dta, clear
drop if _rep == 0
qui ds onset_*
local ons `r(varlist)'
qui tostring `ons', replace force
egen str100 akey = concat(`ons'), punct("|")
qui levelsof akey, local(keys)
dis "distinct assignments drawn = " `: word count `keys''
preserve
contract akey
sum _freq, detail
gen double dev = (_freq - 50)^2/50
qui sum dev
dis "CHI2 = " r(sum) "  df = " _N - 1
dis "p = " chi2tail(_N-1, r(sum))
list, noobs sepby(akey)
restore
dis "PROBE16_DONE"
