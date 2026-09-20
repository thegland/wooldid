clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"
local OO "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/out"
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"
local BASE wooldid y cohort t ttre

* build a row-shuffled copy of dC and a large-mean-outcome copy
use "`DD'/dC.dta", clear
set seed 4242
gen double _shuf = runiform()
sort _shuf
drop _shuf
save "`DD'/dC_shuf.dta", replace

use "`DD'/dC.dta", clear
gen double ybig = 250000 + 30000*y
label var ybig "outcome with a very large mean (dollars)"
save "`DD'/dC_big.dta", replace

*=== rank-deficient event-study covariance in a NON-subgroup run ===
* cluster(cohort): 4 clusters against 9 horizons -> clustered VCE is rank deficient,
* which is exactly where the adversarial line test's 1e-12 rank truncation is live.
wdcell , dta("`DD'/dC.dta") tag(q01_advrank) base(`BASE') out("`OO'") opts(att esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests)
wdcell , dta("`DD'/dC.dta") tag(q02_advband) base(`BASE') out("`OO'") opts(att att_i esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests simultaneousbands pretrendpower(0.5))
wdcell , dta("`DD'/dC.dta") tag(q03_clunit)  base(`BASE') out("`OO'") opts(att att_i esprelength(3) espostlength(3) cluster(unit) jointtests)

*=== forwarding into permute() and bootstrap() ===
wdcell , dta("`DD'/dC.dta") tag(q04_perm)  base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) permute(shuffle) permreps(40) rseed(7))
wdcell , dta("`DD'/dC.dta") tag(q05_boot)  base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) bootstrap(standard) reps(40) rseed(7) cluster(unit))

*=== the ineligible list: both arms must be identical AND fcol==0 in both ===
wdcell , dta("`DD'/dC.dta") tag(q06_lattice) base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) lattice(all 3))
wdcell , dta("`DD'/dC.dta") tag(q07_inter)   base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) interactivecontrols(x1))

*=== very large outcome mean: the help says no displayed digit is affected ===
wdcell , dta("`DD'/dC_big.dta") tag(q08_big) base(wooldid ybig cohort t ttre) out("`OO'") opts(att att_i att_it att_itime esprelength(3) espostlength(3) jointtests)

*=== row-order calibration: help claims shuffling rows moves things ~19x MORE ===
qui use "`DD'/dC.dta", clear
qui wooldid y cohort t ttre, att att_i att_it att_itime esprelength(3) espostlength(3) jointtests nofastcollapse
wddump , file("`OO'/r_orig_nofc.txt")
qui use "`DD'/dC_shuf.dta", clear
qui wooldid y cohort t ttre, att att_i att_it att_itime esprelength(3) espostlength(3) jointtests nofastcollapse
wddump , file("`OO'/r_shuf_nofc.txt")
qui use "`DD'/dC.dta", clear
qui wooldid y cohort t ttre, att att_i att_it att_itime esprelength(3) espostlength(3) jointtests
wddump , file("`OO'/r_orig_fc.txt")
dis "P05_DONE"
