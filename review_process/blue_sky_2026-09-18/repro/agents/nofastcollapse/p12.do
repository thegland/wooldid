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
local ES att esprelength(3) espostlength(3) jointtests

* ===== THE FRAMING CELL =====
* Is the V movement specific to the COLLAPSE, or is margins' delta-method V simply
* row-order sensitive at 1e-6 once a continuous control is present?
* nofastcollapse on the ORIGINAL row order vs a shuffled copy of the same data.
use "`DD'/dCc.dta", clear
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "`DD'/dCc_shuf.dta", replace

qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/y01_orig_nofc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/y01_shuf_nofc.txt")
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/y01_orig_fc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/y01_shuf_fc.txt")

* ===== the two control sets that were dropped on my syntax error =====
wdcell , dta("`DD'/dCc.dta") tag(y02_ccc)     base(`BASE') out("`OO'") opts(`ES' coarsecohortcontrols(c.xdollar))
wdcell , dta("`DD'/dCc.dta") tag(y03_cccabs)  base(`BASE') out("`OO'") opts(`ES' ccc_absorb(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(y04_ccci)    base(`BASE') out("`OO'") opts(`ES' coarsecohortcontrols(i.region))
wdcell , dta("`DD'/dCc.dta") tag(y05_inter)   base(`BASE') out("`OO'") opts(`ES' interactivecontrols(xdollar))
dis "P12_DONE"
