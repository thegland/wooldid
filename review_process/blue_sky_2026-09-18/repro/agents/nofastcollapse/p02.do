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

* ---- continuous treatment, every contreatcontrols mode ----
wdcell , dta("`DD'/dC.dta") tag(k01_ct_basic) base(`BASE') out("`OO'") opts(att att_i att_it att_itime `ES' jointtests contreat(dose) contreatcontrols(basic))
wdcell , dta("`DD'/dC.dta") tag(k02_ct_byT)   base(`BASE') out("`OO'") opts(att att_i `ES' jointtests contreat(dose) contreatcontrols(basic byTime))
wdcell , dta("`DD'/dC.dta") tag(k03_ct_byC)   base(`BASE') out("`OO'") opts(att att_i `ES' jointtests contreat(dose) contreatcontrols(basic byCohort))
wdcell , dta("`DD'/dC.dta") tag(k04_ct_byTt)  base(`BASE') out("`OO'") opts(att att_i `ES' jointtests contreat(dose) contreatcontrols(basic byTtre))
wdcell , dta("`DD'/dC.dta") tag(k05_ct_ccno)  base(`BASE') out("`OO'") opts(att `ES' contreat(dose) contreatcontrols(basic) ccnoabsorb)
wdcell , dta("`DD'/dC.dta") tag(k06_ct_dyex)  base(`BASE') out("`OO'") opts(att `ES' contreat(dose) contreatcontrols(basic) contreatelasticitytype(dyex))
wdcell , dta("`DD'/dC.dta") tag(k07_ct_aw)    base(wooldid y cohort t ttre [aw=w]) out("`OO'") opts(att att_it `ES' contreat(dose) contreatcontrols(basic))

* ---- poisson ----
wdcell , dta("`DD'/dC.dta") tag(k08_pois)     base(wooldid ycount cohort t ttre) out("`OO'") opts(att att_i att_it att_itime `ES' jointtests poisson)
wdcell , dta("`DD'/dC.dta") tag(k09_pois_ct)  base(wooldid ycount cohort t ttre) out("`OO'") opts(att `ES' poisson contreat(dose) contreatcontrols(basic))
wdcell , dta("`DD'/dC.dta") tag(k10_pois_dyex) base(wooldid ycount cohort t ttre) out("`OO'") opts(att `ES' poisson contreat(dose) contreatcontrols(basic) contreatelasticitytype(dyex))
wdcell , dta("`DD'/dC.dta") tag(k11_pois_semi) base(wooldid ycount cohort t ttre) out("`OO'") opts(att `ES' poisson semielasticity)

* ---- customweightmultiplier with exact zeros ----
wdcell , dta("`DD'/dC.dta") tag(k12_cwm0)     base(`BASE') out("`OO'") opts(att att_i att_it att_itime `ES' jointtests customweightmultiplier(cwm0))
wdcell , dta("`DD'/dC.dta") tag(k13_cwm0_ct)  base(`BASE') out("`OO'") opts(att att_i `ES' contreat(dose) contreatcontrols(basic) customweightmultiplier(cwm0))

* ---- subgroups ----
wdcell , dta("`DD'/dC.dta") tag(k14_sg_ct)    base(`BASE') out("`OO'") opts(att `ES' jointtests subgroups(sgc) contreat(dose) contreatcontrols(basic))
wdcell , dta("`DD'/dC.dta") tag(k15_sg_clu)   base(`BASE') out("`OO'") opts(att `ES' jointtests subgroups(sgc) cluster(cohort))

dis "P02_DONE"
