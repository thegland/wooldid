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

use "`DD'/dCc.dta", clear
gen double xs0 = x1
gen double xs2 = 100*x1 + 500
gen double xs4 = 10000*x1 + 50000
gen double xs6 = 1000000*x1 + 5000000
label var xs6 "county population scale"
save "`DD'/dCs.dta", replace

local BASE wooldid y cohort t ttre
local ES att esprelength(3) espostlength(3) jointtests

* --- scale ladder: how does the disagreement scale with the control's magnitude? ---
wdcell , dta("`DD'/dCs.dta") tag(x01_s0) base(`BASE') out("`OO'") opts(`ES' controls(xs0))
wdcell , dta("`DD'/dCs.dta") tag(x02_s2) base(`BASE') out("`OO'") opts(`ES' controls(xs2))
wdcell , dta("`DD'/dCs.dta") tag(x03_s4) base(`BASE') out("`OO'") opts(`ES' controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x04_s6) base(`BASE') out("`OO'") opts(`ES' controls(xs6))

* --- can the residue change WHAT is reported?  near-singular clustered ES covariance ---
wdcell , dta("`DD'/dCs.dta") tag(x05_advc) base(`BASE') out("`OO'") opts(att esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x06_band) base(`BASE') out("`OO'") opts(att att_i esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests simultaneousbands pretrendpower(0.5) controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x07_sgc)  base(`BASE') out("`OO'") opts(att esprelength(3) espostlength(3) jointtests subgroups(sgc) controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x08_imp)  base(`BASE') out("`OO'") opts(att esprelength(3) espostlength(3) impvalues cluster(region) controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x09_pois) base(wooldid ycount cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests poisson controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x10_tt)   base(`BASE') out("`OO'") opts(`ES' timetrends controls(xs4))
dis "P08_DONE"
