clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
do "`W'/pchk2.do"

* ---------- A. continuous treatment x permute ----------
cap prog drop ctchk
prog def ctchk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(unit) permute(shuffle) permreps(30) rseed(777) `opts'
  dis "CT|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "CT|`tag'|cfx_blank=" "`e(perm_cfx_blank)'" "|gt_blank=" "`e(perm_gt_blank)'"
  dis "CT|`tag'|bnames=" "`: colfullnames e(b)'"
  cap confirm matrix e(perm_comparator_b)
  if _rc == 0 mata: wdpcmp2("`tag'")
end
ctchk ct_plain  , opts(contreat(dose))
ctchk ct_basic  , opts(contreat(dose) contreatcontrols(basic))
ctchk ct_es     , opts(contreat(dose) esprelength(2) espostlength(2))
ctchk ct_lat    , opts(contreat(dose) lattice(all 3))
ctchk ct_cfx    , opts(contreat(dose) makecfxplots)
ctchk ct_sg     , opts(contreat(dose) subgroups(sg2))

* ---------- B. permonsetcohorts (i = cohort, cluster finer) ----------
cap prog drop oncchk
prog def oncchk
  syntax anything(name=tag) [, OPTS(string asis) IV(string) CL(string)]
  if "`iv'" == "" local iv cohort
  if "`cl'" == "" local cl unit
  cap noi wooldid y `iv' t ttre , att cluster(`cl') permute(shuffle) permreps(30) rseed(777) `opts'
  dis "ONC|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "ONC|`tag'|irebuilt=" e(perm_irebuilt) "|npossible=" e(perm_npossible) "|done=" e(perm_completed) "|fail=" e(perm_failed)
  cap confirm matrix e(perm_comparator_b)
  if _rc == 0 mata: wdpcmp2("`tag'")
end
oncchk onc_noflag                                         // cluster finer than i, NO assertion -> should refuse
oncchk onc_flag     , opts(permonsetcohorts)
oncchk onc_es       , opts(permonsetcohorts esprelength(2) espostlength(2))
oncchk onc_sg       , opts(permonsetcohorts subgroups(sg2))
oncchk onc_ctrl     , opts(permonsetcohorts controls(zbig))
oncchk onc_ateq     , opts(permonsetcohorts) cl(cohort)   // cluster == i: rebuild must be a no-op
oncchk onc_gt       , opts(permonsetcohorts gteffects)

* ---------- C. permstrata / permsave / weights+if together ----------
cap prog drop mchk
prog def mchk
  syntax anything(name=tag) [, OPTS(string asis) WT(string asis) IFQ(string asis)]
  cap noi wooldid y unit t ttre `wt' `ifq' , att cluster(unit) permute(shuffle) permreps(30) rseed(777) `opts'
  dis "MX|`tag'|RC=" _rc
  if _rc != 0 exit
  cap confirm matrix e(perm_comparator_b)
  if _rc == 0 mata: wdpcmp2("`tag'")
end
mchk strata     , opts(permstrata(region) esprelength(2) espostlength(2))
mchk wt_and_if  , wt([aw=w]) ifq(if keepreg==1) opts(esprelength(2) espostlength(2))
mchk wt_if_sg   , wt([aw=w]) ifq(if keepreg==1) opts(esprelength(2) espostlength(2) subgroups(sg2))
mchk pw_cwm     , wt([pw=w]) opts(customweightmultiplier(fwt) esprelength(2) espostlength(2))
mchk save1      , opts(permsave(p07perm) replace esprelength(2) espostlength(2))
cap confirm file "p07perm.dta"
dis "MX|save1|FILEEXISTS=" (_rc==0)
dis "P07_DONE"
