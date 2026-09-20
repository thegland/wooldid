clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
include dump.doi
local ES esprelength(3) espostlength(3)

* F01 permute (i=id, cluster(clu) which nests id)
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) permute(shuffle) permreps(40) rseed(777)
edump F01_perm `=_rc'

* F02 permute + permtstat
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) permute(shuffle) permreps(40) permtstat rseed(777)
edump F02_perm_tstat `=_rc'

* F03 permute + subgroups + jointtests
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) subgroups(sg) jointtests permute(shuffle) permreps(40) rseed(777)
edump F03_perm_sg_jt `=_rc'

* F04 permute + permanalyticse
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i esfixedbaseperiod `ES' cluster(clu) permute(shuffle) permreps(40) permanalyticse rseed(777)
edump F04_perm_ase `=_rc'

* ============ BLANKING-GUARD AUDIT ============
* identical spec, plain vs bootstrap vs permute: do the analytic joint/adversarial scalars survive?
local SPEC att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) subgroups(sg) jointtests adversarial

use wd_base.dta, clear
cap noi wooldid y id t ttre , `SPEC'
dis "G_PLAIN_RC=" _rc
foreach s in joint_att_allsgssame joint_att_i_allsgssame joint_pre_att_allsgssame joint_horiz joint_horiz_i joint_horiz_s0 joint_horiz_s1 adversarial_pval_att adversarial_pval_att_it {
  dis "GPLAIN|`s'|" %21x (e(`s'))
}
edump G01_plain_jtadv `=_rc'

use wd_base.dta, clear
cap noi wooldid y id t ttre , `SPEC' bootstrap(standard) reps(40) rseed(777)
dis "G_BOOT_RC=" _rc
foreach s in joint_att_allsgssame joint_att_i_allsgssame joint_pre_att_allsgssame joint_horiz joint_horiz_i joint_horiz_s0 joint_horiz_s1 adversarial_pval_att adversarial_pval_att_it {
  dis "GBOOT|`s'|" %21x (e(`s'))
}
edump G02_boot_jtadv `=_rc'

use wd_base.dta, clear
cap noi wooldid y id t ttre , `SPEC' permute(shuffle) permreps(40) rseed(777)
dis "G_PERM_RC=" _rc
foreach s in joint_att_allsgssame joint_att_i_allsgssame joint_pre_att_allsgssame joint_horiz joint_horiz_i joint_horiz_s0 joint_horiz_s1 adversarial_pval_att adversarial_pval_att_it {
  dis "GPERM|`s'|" %21x (e(`s'))
}
edump G03_perm_jtadv `=_rc'

* ============ bootstrap with reps < ncoef : rank of e(V) ============
use wd_base.dta, clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) bootstrap(standard) reps(8) rseed(777)
dis "H_RC=" _rc
if _rc == 0 {
  dis "H_K=" colsof(e(b)) " boot_completed=" e(boot_completed) " boot_N_V=" e(boot_N_V) " boot_v_posted=" e(boot_v_posted)
  mata: V=st_matrix("e(V)"); ev=Re(eigenvalues(V)); printf("H_rank=%g H_mineig=%g H_maxeig=%g\n", rank(V), min(ev), max(ev))
  cap noi test ES_att_post0 ES_att_post1 ES_att_post2
  dis "H_TEST_p=" r(p) " df=" r(df)
  cap noi test main_att main_att_i ES_att_pre4 ES_att_pre3 ES_att_pre2 ES_att_post0 ES_att_post1 ES_att_post2
  dis "H_TESTBIG_p=" r(p) " df=" r(df)
}
edump H01_boot_reps8 `=_rc'

dis "PROBEF_DONE"
