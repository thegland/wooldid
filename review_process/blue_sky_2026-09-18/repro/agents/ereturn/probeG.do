clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
include dump.doi
local ES esprelength(3) espostlength(3)

* permute-legal panel: cohorts are contiguous id blocks; clusters nest inside one cohort
clear
set obs 80
gen long id = _n
gen byte grp = floor((id-1)/20)
gen double ttre = .
replace ttre = 4 if grp==1
replace ttre = 6 if grp==2
replace ttre = 8 if grp==3
gen long clu = grp*10 + ceil(mod(id-1,20)/4)
gen byte sg = mod(id,2)
gen double ufe = rnormal()
expand 10
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + 0.3*t + 1.5*post + rnormal()
save wd_perm.dta, replace
qui tab clu grp
dis "NCLUS=" r(r)

* P01 permute shuffle
use wd_perm.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) permute(shuffle) permreps(40) rseed(777)
dis "P01_RC=" _rc
edump P01_perm `=_rc'

* P02 permute + permtstat
use wd_perm.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) permute(shuffle) permreps(40) permtstat rseed(777)
dis "P02_RC=" _rc
edump P02_perm_tstat `=_rc'

* P03 permute + permanalyticse
use wd_perm.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod `ES' cluster(clu) permute(shuffle) permreps(40) permanalyticse rseed(777)
dis "P03_RC=" _rc
edump P03_perm_ase `=_rc'

* P04 permute + subgroups
use wd_perm.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) subgroups(sg) permute(shuffle) permreps(40) rseed(777)
dis "P04_RC=" _rc
edump P04_perm_sg `=_rc'

dis "PROBEG_DONE"
