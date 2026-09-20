clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
cap erase edump.txt
include dump.doi

* ---------- sanity: does Stata accept >32-char matrix colnames? ----------
tempname T
matrix `T' = (1,2)
cap matrix colnames `T' = abcdefghij0123456789012345678901 xx
dis "NAME32_RC=" _rc
matrix `T' = (1,2)
cap matrix colnames `T' = abcdefghij01234567890123456789012 xx
dis "NAME33_RC=" _rc
matrix `T' = (1,2)
cap matrix colnames `T' = abcdefghij01234567890123456789012 xx
if _rc==0 {
  local cn : colnames `T'
  dis "NAME33_STORED=`cn'"
}
tempname S
matrix `S' = 1
cap ereturn scalar thisisaverylongname012345678901234 = 1
dis "ERETSCA33_RC=" _rc
ereturn clear

include mkdata.doi
save wd_base.dta, replace

* ================= grid =================
local k = 0

* A01 default (att only)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att
edump A01_default `=_rc'

* A02 all four estimands
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime
edump A02_all4 `=_rc'

* A03 all4 + esfixedbaseperiod
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod
edump A03_esfixed `=_rc'

* A04 fastnose
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime fastnose
edump A04_fastnose `=_rc'

* A05 subgroups
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime subgroups(sg)
edump A05_sg `=_rc'

* A06 contreat
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime contreat(dose)
edump A06_contreat `=_rc'

* A07 gteffects
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i gteffects
edump A07_gteffects `=_rc'

* A08 simultaneousbands
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i simultaneousbands
edump A08_simband `=_rc'

* A09 jointtests + adversarial
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i jointtests adversarial
edump A09_jtadv `=_rc'

* A10 cluster on i (=id)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i cluster(id)
edump A10_clusid `=_rc'

* A11 cohort-level i, cluster finer (id)
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i cluster(id)
edump A11_finerclus `=_rc'

* A12 robust
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i robust
edump A12_robust `=_rc'

* A13 unconditionalse
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i unconditionalse
edump A13_uncond `=_rc'

* A14 customvce
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i customvce(cluster id)
edump A14_customvce `=_rc'

* A15 impvalues
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i impvalues cluster(id)
edump A15_impvalues `=_rc'

* A16 poisson
use wd_base.dta, clear
cap noi wooldid ypos id t ttre , att att_i poisson
edump A16_poisson `=_rc'

* A17 timetrends
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i timetrends
edump A17_timetrends `=_rc'

* A18 anticipation
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i anticipation(1)
edump A18_antic `=_rc'

* A19 summarystats + histogram
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i summarystats histogramcohorteffects
edump A19_sumhist `=_rc'

* A20 pretrendpower
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod pretrendpower(0.5)
edump A20_ptp `=_rc'

* A21 nofastcollapse
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i nofastcollapse
edump A21_nofc `=_rc'

* A22 fw weights
use wd_base.dta, clear
cap noi wooldid y id t ttre [fw=fwt], att att_i
edump A22_fw `=_rc'

* A23 aw weights
use wd_base.dta, clear
cap noi wooldid y id t ttre [aw=awt], att att_i
edump A23_aw `=_rc'

* A24 subgroups + contreat + all4 (name stress)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime subgroups(sg) contreat(dose)
edump A24_sgcont `=_rc'

* A25 level(90)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i level(90)
edump A25_level90 `=_rc'

* A26 level(97.5)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i level(97.5)
edump A26_level975 `=_rc'

* A27 fastnose + subgroups
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i subgroups(sg) fastnose
edump A27_fastnose_sg `=_rc'

* A28 fastnose + contreat
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i contreat(dose) fastnose
edump A28_fastnose_ct `=_rc'

dis "PROBEA_DONE"
