clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* ---- a realistic-ish staggered panel: 150 units, 12 periods, 3 treated cohorts + never ----
set obs 150
gen id = _n
gen byte grp = mod(_n,5)
gen ttre = .
replace ttre = 5 if grp==1
replace ttre = 7 if grp==2
replace ttre = 9 if grp==3
gen byte sg = mod(_n,2)
gen double ui = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort id: gen t = _n
gen double y = 2 + ui + 0.1*t + rnormal() + 2*(t>=ttre & ttre<.)*(1+0.3*(t-ttre))
gen double w = 1 + runiform()
compress
save panel.dta, replace

di "WDCFG_BEGIN A"
cap noisily wooldid y id t ttre, att
di "WDRC=" _rc
di "WDCFG_END A"

di "WDCFG_BEGIN B"
cap noisily wooldid y id t ttre, att att_i att_it att_itime esprelength(3) espostlength(3) cluster(id) jointtests
di "WDRC=" _rc
di "WDCFG_END B"

di "WDCFG_BEGIN C"
cap noisily wooldid y id t ttre, att att_i esprelength(3) espostlength(3) cluster(id) jointtests subgroups(sg) adversarial esfixedbaseperiod pretrendpower(0.5) gteffects summarystats
di "WDRC=" _rc
di "WDCFG_END C"

di "WDCFG_BEGIN D"
cap noisily wooldid y id t ttre, att esprelength(2) espostlength(2) cluster(id) permute(shuffle) permreps(25)
di "WDRC=" _rc
di "WDCFG_END D"

di "WDCFG_BEGIN E"
cap noisily wooldid y id t ttre, att esprelength(2) espostlength(2) cluster(id) bootstrap(standard) reps(25) boottstat dots(0)
di "WDRC=" _rc
di "WDCFG_END E"

di "WDCFG_BEGIN F"
cap noisily wooldid y id t ttre [aw=w], att att_i esprelength(2) espostlength(2) cluster(id) contreat(dose) contreatcontrols(byTime) lattice(quantiles 3) makecfxplots
di "WDRC=" _rc
di "WDCFG_END F"

di "SENTINEL_DONE"
