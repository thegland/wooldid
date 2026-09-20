clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
include dump.doi
local ES esprelength(3) espostlength(3)

* ---------------- contreat with coarse i ----------------
use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(untreatedZero)
edump C01_ct0 `=_rc'

use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(basic)
edump C02_ctbasic `=_rc'

use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(untreatedZero) subgroups(sg)
edump C03_ct_sg `=_rc'

use wd_base.dta, clear
cap noi wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(untreatedZero) fastnose
edump C04_ct_fastnose `=_rc'

* ---------------- naming stress: 7-digit and 8-digit subgroup values ------------
use wd_base.dta, clear
gen long sg7 = 1000000 + mod(id,2)
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' subgroups(sg7)
dis "C05_RC=" _rc
edump C05_sg7 `=_rc'

use wd_base.dta, clear
gen long sg8 = 10000000 + mod(id,2)
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' subgroups(sg8)
dis "C06_RC=" _rc
edump C06_sg8 `=_rc'

* 8-digit but WITHOUT att_itime (shorter suffix) - does it survive?
use wd_base.dta, clear
gen long sg8 = 10000000 + mod(id,2)
cap noi wooldid y id t ttre , att esfixedbaseperiod `ES' subgroups(sg8)
dis "C07_RC=" _rc
edump C07_sg8_attonly `=_rc'

* 9-digit
use wd_base.dta, clear
gen long sg9 = 100000000 + mod(id,2)
cap noi wooldid y id t ttre , att esfixedbaseperiod `ES' subgroups(sg9)
dis "C08_RC=" _rc
edump C08_sg9_attonly `=_rc'

* ---------------- staleness ----------------
use wd_base.dta, clear
qui wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' subgroups(sg) gteffects simultaneousbands jointtests
edump C10_wide `=_rc'
local cl1 `"`e(cmdline)'"'

* now a NARROWER successful run
cap noi wooldid y id t ttre , att
edump C11_narrow `=_rc'

* now a FAILING run on top of a successful one
use wd_base.dta, clear
qui wooldid y id t ttre , att att_i esfixedbaseperiod `ES' gteffects
local cl_good `"`e(cmdline)'"'
cap noi wooldid y id t ttre , att att_i att_it subgroups(nonexistentvar)
dis "C12_FAILRC=" _rc
dis `"C12_CMDLINE_AFTER_FAIL=`e(cmdline)'"'
dis "C12_CMD=`e(cmd)'"
edump C12_afterfail `=_rc'

* replay refusal after failure
cap noi wooldid
dis "C13_REPLAY_AFTER_FAIL_RC=" _rc

* clear the session marker the way a user might, then replay
macro drop WOOLDID_LASTRC
cap noi wooldid
dis "C14_REPLAY_AFTER_MACRODROP_RC=" _rc
dis `"C14_CMDLINE=`e(cmdline)'"'

* macro drop _all variant
use wd_base.dta, clear
qui wooldid y id t ttre , att att_i esfixedbaseperiod `ES'
cap noi wooldid y id t ttre , att att_i subgroups(nonexistentvar)
macro drop _all
cap noi wooldid
dis "C15_REPLAY_AFTER_MACRODROPALL_RC=" _rc

* ---------------- sentinel semantics ----------------
* what the -1e100 sentinel does to test / lincom / ereturn display
clear
set obs 100
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()
qui reg yy x1 x2
tempname bq vq
matrix `bq' = (1.0, 2.0)
matrix `vq' = (0.25, 0 \ 0, -1e100)
matrix colnames `bq' = good bad
matrix colnames `vq' = good bad
matrix rownames `vq' = good bad
ereturn post `bq' `vq', obs(100) dof(97)
dis "SENT_POST_RC=" _rc
dis "--- ereturn display ---"
cap noi ereturn display
dis "--- test bad ---"
cap noi test bad
dis "SENT_TEST_BAD_RC=" _rc
dis "SENT_TEST_BAD_p=" r(p) " chi2/F=" r(F) r(chi2)
dis "--- test good bad (joint) ---"
cap noi test good bad
dis "SENT_TESTJOINT_RC=" _rc
dis "SENT_TESTJOINT_p=" r(p)
dis "--- test good = bad ---"
cap noi test good = bad
dis "SENT_TESTEQ_RC=" _rc
dis "SENT_TESTEQ_p=" r(p)
dis "--- lincom good + bad ---"
cap noi lincom good + bad
dis "SENT_LINCOM_RC=" _rc  " se=" r(se)
dis "--- test good alone (should be unaffected) ---"
cap noi test good
dis "SENT_TESTGOOD_p=" r(p)
ereturn clear

dis "PROBEC_DONE"
