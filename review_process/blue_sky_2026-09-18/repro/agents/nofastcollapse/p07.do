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

* --- build the controls test bed on top of dC ---
use "`DD'/dC.dta", clear
* xcell is ALREADY constant within a (unit,t) cell but varies within a (cohort,t) cell.
* build a control that is constant within the (cohort,t) cell -- the collapse's own cell.
bysort cohort t (unit): gen double xcohortc = x1[1]
* a control on an administrative scale: household income in dollars
gen double xdollar = 52000 + 9000*x1 + 300*x2
* a second one, and a percentage-scale one
gen double xdollar2 = 31000 - 4000*x2
gen double xpct = 50 + 8*x1
label var xdollar "control on a dollar scale"
save "`DD'/dCc.dta", replace

local BASE wooldid y cohort t ttre
local ES att att_i att_it esprelength(3) espostlength(3) jointtests

wdcell , dta("`DD'/dCc.dta") tag(w01_c1)      base(`BASE') out("`OO'") opts(`ES' controls(x1))
wdcell , dta("`DD'/dCc.dta") tag(w02_cell)    base(`BASE') out("`OO'") opts(`ES' controls(xcohortc))
wdcell , dta("`DD'/dCc.dta") tag(w03_dollar)  base(`BASE') out("`OO'") opts(`ES' controls(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(w04_dollar2) base(`BASE') out("`OO'") opts(`ES' controls(xdollar xdollar2))
wdcell , dta("`DD'/dCc.dta") tag(w05_pct)     base(`BASE') out("`OO'") opts(`ES' controls(xpct))
wdcell , dta("`DD'/dCc.dta") tag(w06_fe)      base(`BASE') out("`OO'") opts(`ES' fe(region))
wdcell , dta("`DD'/dCc.dta") tag(w07_fvctrl)  base(`BASE') out("`OO'") opts(`ES' controls(i.region))
wdcell , dta("`DD'/dCc.dta") tag(w08_ccc)     base(`BASE') out("`OO'") opts(`ES' coarsecohortcontrols(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(w09_ct)      base(`BASE') out("`OO'") opts(att att_i esprelength(3) espostlength(3) jointtests contreat(dose) contreatcontrols(basic) controls(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(w10_rob)     base(`BASE') out("`OO'") opts(`ES' controls(xdollar) robust)
wdcell , dta("`DD'/dCc.dta") tag(w11_clu)     base(`BASE') out("`OO'") opts(`ES' controls(xdollar) cluster(region))

* --- does any DISPLAYED digit move?  print both tables, marked. ---
dis "=========== DISPLAY_A (collapse on) ==========="
use "`DD'/dCc.dta", clear
noisily wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xdollar xdollar2)
dis "=========== DISPLAY_B (nofastcollapse) ==========="
use "`DD'/dCc.dta", clear
noisily wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xdollar xdollar2) nofastcollapse
dis "=========== DISPLAY_END ==========="
dis "P07_DONE"
