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

* ---- dP: the COMMONEST wooldid shape. balanced panel, i = unit, ONE row per (i,t) cell.
* the collapse then aggregates nothing: it is a pure re-sort of the estimation sample.
local gseed 90210
local gN 40
local gT 10
local gonsets 2004 2006 2008
local gnever 0.35
local gysd 0.30
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
save "`DD'/dP.dta", replace
* a row-shuffled copy of the SAME data
set seed 31337
gen double _s = runiform()
sort _s
drop _s
save "`DD'/dP_shuf.dta", replace

local OPT att att_i att_it att_itime esprelength(3) espostlength(3) jointtests

* (a) collapse vs nofastcollapse on the same data
wdcell , dta("`DD'/dP.dta") tag(z01_onerow) base(wooldid y unit t ttre) out("`OO'") opts(`OPT')

* (b) the pure row-permutation reference: nofastcollapse on original vs shuffled
qui use "`DD'/dP.dta", clear
qui wooldid y unit t ttre, `OPT' nofastcollapse
wddump , file("`OO'/z02_orig_nofc.txt")
qui use "`DD'/dP_shuf.dta", clear
qui wooldid y unit t ttre, `OPT' nofastcollapse
wddump , file("`OO'/z02_shuf_nofc.txt")
* (c) and the collapse arm on the shuffled data -- if the collapse's sort is what pins
*     row order, the two COLLAPSED runs should agree far more tightly than the two
*     uncollapsed ones do.
qui use "`DD'/dP.dta", clear
qui wooldid y unit t ttre, `OPT'
wddump , file("`OO'/z03_orig_fc.txt")
qui use "`DD'/dP_shuf.dta", clear
qui wooldid y unit t ttre, `OPT'
wddump , file("`OO'/z03_shuf_fc.txt")

* (d) how many cells did the collapse make?  verbose says so.
qui use "`DD'/dP.dta", clear
noisily wooldid y unit t ttre, att esprelength(2) espostlength(2) verbose
dis "Z_fcncell=" e(wdid_fcncell) " N=" e(N)
dis "P06_DONE"
