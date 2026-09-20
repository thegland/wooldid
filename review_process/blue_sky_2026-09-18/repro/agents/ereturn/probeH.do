clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
include dump.doi
local ES esprelength(3) espostlength(3)
local SPEC att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) subgroups(sg) jointtests adversarial
local SCAL joint_att_allsgssame joint_att_i_allsgssame joint_pre_att_allsgssame joint_horiz joint_horiz_i joint_horiz_s0 joint_horiz_s1 adversarial_pval_att adversarial_pval_att_it

use wd_perm.dta, clear
ereturn clear
cap noi wooldid y id t ttre , `SPEC'
dis "Q_PLAIN_RC=" _rc
foreach s of local SCAL {
  dis "QPLAIN|`s'|" %21x (e(`s'))
}
edump Q01_plain `=_rc'

use wd_perm.dta, clear
ereturn clear
cap noi wooldid y id t ttre , `SPEC' permute(shuffle) permreps(40) rseed(777)
dis "Q_PERM_RC=" _rc
foreach s of local SCAL {
  dis "QPERM|`s'|" %21x (e(`s'))
}
dis "QPERM_cmd=`e(cmd)' perm_method=`e(perm_method)'"
edump Q02_perm `=_rc'

* reps(8) vs reps(40): rank of e(V) and what the low-rep caution says
use wd_perm.dta, clear
ereturn clear
cap noi wooldid y id t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' cluster(clu) bootstrap(standard) reps(8) rseed(777)
dis "R8_RC=" _rc
if _rc==0 {
  mata: V=st_matrix("e(V)"); printf("R8 K=%g rank=%g mineig=%g\n", cols(V), rank(V), min(Re(eigenvalues(V))))
  dis "R8 boot_N_V=" e(boot_N_V) " boot_completed=" e(boot_completed) " vcvstate=`e(wdid_vcvstate)' properties=`e(properties)'"
  cap noi test main_att main_att_i
  dis "R8_TEST2_p=" r(p) " df=" r(df)
}
dis "PROBEH_DONE"
