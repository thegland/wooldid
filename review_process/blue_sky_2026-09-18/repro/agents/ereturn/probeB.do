clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
include dump.doi
local ES esprelength(3) espostlength(3)

* B01 event study, all4, fixed base
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES'
edump B01_es `=_rc'

* B02 same + fastnose  <-- template bug check
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' fastnose
edump B02_es_fastnose `=_rc'

* B03 same + subgroups
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' subgroups(sg)
edump B03_es_sg `=_rc'

* B04 same + contreat untreatedZero
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(untreatedZero)
edump B04_es_ct0 `=_rc'

* B05 contreat basic
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(basic)
edump B05_es_ctbasic `=_rc'

* B06 gteffects
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod `ES' gteffects
edump B06_es_gt `=_rc'

* B07 simultaneousbands
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' simultaneousbands
edump B07_es_simband `=_rc'

* B08 pretrendpower
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod `ES' pretrendpower(0.5)
edump B08_es_ptp `=_rc'

* B09 timetrends
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod `ES' timetrends
edump B09_es_tt `=_rc'

* B10 impvalues (i=coh, cluster finer)
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i esfixedbaseperiod `ES' impvalues cluster(id)
edump B10_es_imp `=_rc'

* B11 jointtests + adversarial
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' jointtests adversarial
edump B11_es_jtadv `=_rc'

* B12 poisson
use wd_base.dta, clear
cap noi wooldid ypos id t ttre , att att_i esfixedbaseperiod `ES' poisson
edump B12_es_pois `=_rc'

* B13 cluster finer than i
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(id)
edump B13_es_finerclus `=_rc'

* B14 non-fixed base (pooled reference)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime `ES'
edump B14_es_pooled `=_rc'

* B15 sg + contreat + ES (widest surface)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' subgroups(sg) contreat(dose) contreatcontrols(untreatedZero) jointtests
edump B15_wide `=_rc'

* B16 ES + robust
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' robust
edump B16_es_robust `=_rc'

* B17 ES + sg + simultaneousbands
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' subgroups(sg) simultaneousbands
edump B17_es_sg_simband `=_rc'

* B18 ES + anticipation
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod `ES' anticipation(1)
edump B18_es_antic `=_rc'

dis "PROBEB_DONE"
