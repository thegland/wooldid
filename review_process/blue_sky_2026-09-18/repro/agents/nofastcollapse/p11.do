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

* a NOISELESS twin of dCc: every cluster's effect is identical, so the IM
* degeneracy guard (se < 1e-10*|b|) sits at a knife edge.
use "`DD'/dCc.dta", clear
gen double ydet = 2*t + 3*cohort + 1.5*(ttre<. & t>=ttre) + 0.5*x1
save "`DD'/dCdet.dta", replace

local BASE wooldid y cohort t ttre
* --- does a 1e-6 rounding move survive into a RANK-based p-value? ---
wdcell , dta("`DD'/dCc.dta") tag(v01_permex) base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) cluster(cohort) permute(shuffle) permexact controls(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(v02_perm)   base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) cluster(region) permute(shuffle) permreps(200) rseed(1) controls(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(v03_boot)   base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) cluster(region) bootstrap(standard) reps(200) rseed(1) controls(xdollar))
* --- Ibragimov-Muller: the degeneracy guard is a THRESHOLD on the cross-cluster spread ---
wdcell , dta("`DD'/dCc.dta") tag(v04_imp)    base(`BASE') out("`OO'") opts(att esprelength(2) espostlength(2) impvalues cluster(region) controls(xdollar))
wdcell , dta("`DD'/dCdet.dta") tag(v05_impdet) base(wooldid ydet cohort t ttre) out("`OO'") opts(att esprelength(2) espostlength(2) impvalues cluster(region) controls(xdollar))
wdcell , dta("`DD'/dCdet.dta") tag(v06_det)  base(wooldid ydet cohort t ttre) out("`OO'") opts(att att_i esprelength(3) espostlength(3) jointtests controls(xdollar))
dis "P11_DONE"
