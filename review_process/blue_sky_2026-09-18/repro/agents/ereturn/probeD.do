clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
include dump.doi
local ES esprelength(3) espostlength(3)

* D01 bootstrap standard
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) bootstrap(standard) reps(40) rseed(777)
edump D01_boot `=_rc'

* D02 bootstrap boottstat
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) bootstrap(standard) reps(40) boottstat rseed(777)
edump D02_boottstat `=_rc'

* D03 bootstrap bayes
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) bootstrap(bayes) reps(40) rseed(777)
edump D03_bootbayes `=_rc'

* D04 bootstrap + subgroups
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) subgroups(sg) bootstrap(standard) reps(40) rseed(777)
edump D04_boot_sg `=_rc'

* D05 bootstrap + simultaneousbands + gteffects
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i esfixedbaseperiod `ES' cluster(id) simultaneousbands gteffects bootstrap(standard) reps(40) rseed(777)
edump D05_boot_band_gt `=_rc'

* D06 permute shuffle
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) permute(shuffle) permreps(40) rseed(777)
edump D06_perm `=_rc'

* D07 permute + permtstat
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) permute(shuffle) permreps(40) permtstat rseed(777)
edump D07_perm_tstat `=_rc'

* D08 permute + subgroups
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id) subgroups(sg) permute(shuffle) permreps(40) rseed(777)
edump D08_perm_sg `=_rc'

* D09 permute + permanalyticse
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i esfixedbaseperiod `ES' cluster(id) permute(shuffle) permreps(40) permanalyticse rseed(777)
edump D09_perm_ase `=_rc'

* D10 permute + gteffects
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i esfixedbaseperiod `ES' cluster(id) gteffects permute(shuffle) permreps(40) rseed(777)
edump D10_perm_gt `=_rc'

dis "PROBED_DONE"
