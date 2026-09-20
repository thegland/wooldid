* ============ ATOMICITY: what does e() look like when an OPTIONAL block dies? ============
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
include "estate.doi"
global CKCOUT "vioT.txt"
cap erase "vioT.txt"
include "ckc.doi"
include "ckc2.doi"

cap prog drop STATE
prog def STATE
  args tag rc
  local vk = .
  cap confirm matrix e(V)
  if _rc == 0 local vk = colsof(e(V))
  local rk = .
  cap mata: st_local("rk", strofreal(rank(st_matrix("e(V)"))))
  local bk = .
  cap confirm matrix e(b)
  if _rc == 0 local bk = colsof(e(b))
  dis "ATOM|`tag'|rc=`rc'|cmd=`e(cmd)'|K=`bk'|VK=`vk'|rank=`rk'|vcvstate=`e(wdid_vcvstate)'|props=`e(properties)'|vce=`e(vce)'|bootmeth=`e(boot_method)'|permmeth=`e(perm_method)'|bootcomp=`e(boot_completed)'|bootfail=`e(boot_failed)'|permcomp=`e(perm_completed)'|lastrc=${WOOLDID_LASTRC}"
end

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"

* ---------------- 1. bootsave() to an existing file with no replace ----------------
dis "########## 1. bootsave onto an existing file ##########"
use D1.dta, clear
cap erase bsv.dta
qui save bsv.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(20) rseed(7) bootsave(bsv)
local rc1 = _rc
dis "BOOTSAVE_EXISTING_RC=`rc1'"
STATE BOOTSAVE_EXIST `rc1'
cap noi ckc BOOTSAVE_EXIST
cap noi wooldid
dis "REPLAY_AFTER_BOOTSAVEFAIL_RC=" _rc

* ---------------- 2. permsave() to an existing file with no replace ----------------
dis "########## 2. permsave onto an existing file ##########"
use D1.dta, clear
cap erase psv.dta
qui save psv.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) permute(shuffle) permreps(20) rseed(7) permsave(psv)
local rc2 = _rc
dis "PERMSAVE_EXISTING_RC=`rc2'"
STATE PERMSAVE_EXIST `rc2'
cap noi wooldid
dis "REPLAY_AFTER_PERMSAVEFAIL_RC=" _rc

* ---------------- 3. a bootstrap where every replication fails ----------------
dis "########## 3. all replications fail ##########"
use D1.dta, clear
gen byte clu1 = 1
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu1) bootstrap(standard) reps(6) rseed(7)
local rc3 = _rc
dis "ALLREPFAIL_RC=`rc3'"
STATE ALLREPFAIL `rc3'
cap noi ckc ALLREPFAIL

* ---------------- 4. low reps: the rank-deficient replication VCE (F7) ----------------
dis "########## 4. rank-deficient bootstrap VCE ##########"
use D1.dta, clear
foreach R in 5 8 12 20 40 {
  cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(`R') rseed(777)
  dis "BOOTR`R'_RC=" _rc
  STATE BOOTR`R' `=_rc'
  cap noi test ES_att_post0 ES_att_post1 ES_att_post2
  dis "BOOTR`R'_TEST3 p=" r(p) " df=" r(df) " F=" r(F) " chi2=" r(chi2) " rc=" _rc
  cap noi test main_att main_att_i ES_att_pre3 ES_att_pre2 ES_att_post0 ES_att_post1 ES_att_post2
  dis "BOOTR`R'_TEST7 p=" r(p) " df=" r(df) " rc=" _rc
  cap noi ckc BOOTR`R'
}
* the same with subgroups (K=60 at reps 40 -> rank 39)
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) subgroups(sg) bootstrap(standard) reps(40) rseed(777)
dis "BOOTSG40_RC=" _rc
STATE BOOTSG40 `=_rc'
cap noi test main_SG0_att main_SG1_att main_SG0_att_i main_SG1_att_i
dis "BOOTSG40_TEST p=" r(p) " df=" r(df) " rc=" _rc
cap noi ckc BOOTSG40

* ---------------- 5. impvalues: does e(df_r) match the table's own df? ----------------
dis "########## 5. impvalues df ##########"
use D1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) impvalues
dis "IMP_RC=" _rc
dis "IMP df_r=" e(df_r) " imdf_min=" e(imdf_min) " imdf_max=" e(imdf_max) " df_r_analytic=" e(df_r_analytic)
cap noi ckc IMPV
cap noi ckc2 IMPV
* the same quantity two ways: the printed table vs lincom
tempname MM
matrix `MM' = e(wdidmainresults)
dis "IMP main_att table: est=" `MM'[1,1] " se=" `MM'[1,2] " lb=" `MM'[1,5] " ub=" `MM'[1,6]
cap noi lincom main_att
dis "IMP lincom main_att: est=" r(estimate) " se=" r(se) " lb=" r(lb) " ub=" r(ub) " rc=" _rc

* the same check on a plain clustered run
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu)
dis "CLU_RC=" _rc " df_r=" e(df_r)
matrix `MM' = e(wdidmainresults)
dis "CLU main_att table: est=" `MM'[1,1] " se=" `MM'[1,2] " lb=" `MM'[1,5] " ub=" `MM'[1,6]
cap noi lincom main_att
dis "CLU lincom main_att: est=" r(estimate) " se=" r(se) " lb=" r(lb) " ub=" r(ub) " rc=" _rc
cap noi ckc CLUDF

* and under poisson, where e(df_r) is reportedly absent
cap noi wooldid ypos id t ttre, att att_i att_it att_itime `ES'
dis "POIS_RC=" _rc " df_r=" e(df_r) " [missing means absent]"
matrix `MM' = e(wdidmainresults)
dis "POIS main_att table: est=" `MM'[1,1] " se=" `MM'[1,2] " lb=" `MM'[1,5] " ub=" `MM'[1,6]
cap noi lincom main_att
dis "POIS lincom main_att: est=" r(estimate) " se=" r(se) " lb=" r(lb) " ub=" r(ub) " rc=" _rc
cap noi ckc POISDF

dis "ATOMIC_DONE"
