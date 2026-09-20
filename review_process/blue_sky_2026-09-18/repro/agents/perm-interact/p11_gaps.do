clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear

* ---- 1. customvce + permute with NO cluster(): is ado:128 reachable at all? ----
cap noi wooldid y unit t ttre , att customvce(cluster unit) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "GAP|cvce_truly_nocluster|RC=" _rc
cap noi wooldid y unit t ttre , att customvce(robust) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "GAP|cvce_robust_nocluster|RC=" _rc

* ---- 2. the help says fastnose+adversarial/jointtests is accepted ONLY under the
*         permute() fast-mode default (no permtstat, NO permanalyticse). Check. ----
cap prog drop g2
prog def g2
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0) `opts'
  dis "GAP|`tag'|RC=" _rc "|perm_analyticse=" e(perm_analyticse)
end
g2 fn_jt                , opts(fastnose jointtests)
g2 fn_adv               , opts(fastnose adversarial)
g2 fn_ase_jt            , opts(fastnose permanalyticse jointtests)
g2 fn_ase_adv           , opts(fastnose permanalyticse adversarial)
g2 ase_jt               , opts(permanalyticse jointtests)
g2 fn_ase_sd            , opts(fastnose permanalyticse permstepdown)
g2 fn_ase_gt            , opts(fastnose permanalyticse gteffects)
g2 fn_ase_simband       , opts(fastnose permanalyticse simultaneousbands)
g2 fn_ase_cfx           , opts(fastnose permanalyticse makecfxplots)

* ---- 3. does the same hole exist on the bootstrap twin? ----
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) bootstrap(standard) reps(30) rseed(777) dots(0) fastnose
dis "GAP|boot_fn|RC=" _rc
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) bootstrap(standard) reps(30) rseed(777) dots(0) fastnose boottstat
dis "GAP|boot_fn_ts|RC=" _rc
dis "P11_DONE"
