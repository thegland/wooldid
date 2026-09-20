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
local ES esprelength(3) espostlength(3)

* ---- control-set variants: things the help says are eligible ----
wdcell , dta("`DD'/dC.dta") tag(m01_timetrend) base(`BASE') out("`OO'") opts(att att_i `ES' jointtests timetrends)
wdcell , dta("`DD'/dC.dta") tag(m02_ccc)       base(`BASE') out("`OO'") opts(att `ES' jointtests coarsecohortcontrols(x1))
wdcell , dta("`DD'/dC.dta") tag(m03_feunit)    base(`BASE') out("`OO'") opts(att att_i `ES' jointtests fe(unit))
wdcell , dta("`DD'/dC.dta") tag(m04_esfbp)     base(`BASE') out("`OO'") opts(att att_i `ES' jointtests esfixedbaseperiod adversarial)
wdcell , dta("`DD'/dC.dta") tag(m05_antic)     base(`BASE') out("`OO'") opts(att `ES' jointtests anticipation(1) oostreatedcontrols(2))
wdcell , dta("`DD'/dC.dta") tag(m06_fastnose)  base(`BASE') out("`OO'") opts(att att_i att_it att_itime `ES' fastnose)
wdcell , dta("`DD'/dC.dta") tag(m07_hist)      base(`BASE') out("`OO'") opts(att att_i att_it att_itime `ES' histogramcohorteffects summarystats gteffects)
wdcell , dta("`DD'/dC.dta") tag(m08_ifres)     base(wooldid y cohort t ttre if x1 > -1) out("`OO'") opts(att att_i `ES' jointtests)
wdcell , dta("`DD'/dC.dta") tag(m09_esrel)     base(`BASE') out("`OO'") opts(att `ES' jointtests esrelativeto(-2))

* ---- messy panel: entry/exit/gaps/missing y and x ----
wdcell , dta("`DD'/dCm.dta") tag(m10_messy)    base(`BASE') out("`OO'") opts(att att_i att_it att_itime `ES' jointtests controls(x1 x2))
wdcell , dta("`DD'/dCm.dta") tag(m11_messy_ct) base(`BASE') out("`OO'") opts(att att_i `ES' jointtests contreat(dose) contreatcontrols(basic) controls(x1))
wdcell , dta("`DD'/dCm.dta") tag(m12_messy_cwm) base(`BASE') out("`OO'") opts(att att_i `ES' jointtests customweightmultiplier(cwm0))

* ---- repeated cross-sections, i = unit ----
wdcell , dta("`DD'/dU.dta") tag(m13_rcs)       base(wooldid y unit t ttre) out("`OO'") opts(att att_i att_it att_itime `ES' jointtests)
wdcell , dta("`DD'/dU.dta") tag(m14_rcs_ct)    base(wooldid y unit t ttre) out("`OO'") opts(att att_i `ES' contreat(dose) contreatcontrols(basic))
wdcell , dta("`DD'/dU.dta") tag(m15_rcs_sg)    base(wooldid y unit t ttre) out("`OO'") opts(att `ES' jointtests subgroups(sgc) cluster(region))
dis "P03_DONE"
