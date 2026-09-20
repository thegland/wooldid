clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear

cap prog drop ochk
prog def ochk
  syntax anything(name=tag) [, OPTS(string asis) WT(string asis) YV(string)]
  if "`yv'" == "" local yv y
  cap noi wooldid `yv' unit t ttre `wt' , att cluster(unit) permute(shuffle) permreps(30) rseed(777) dots(0) `opts'
  dis "ODD|`tag'|RC=" _rc
end
ochk plots_sg    , opts(makeplots saveplots(p09plt) replace subgroups(sg2) esprelength(2) espostlength(2))
ochk pois_aw_es  , yv(ycount) wt([aw=w]) opts(poisson esprelength(2) espostlength(2))
ochk tt_jt_adv   , opts(timetrends jointtests adversarial esprelength(2) espostlength(2))
ochk trip999     , opts(subgroups(sg2) esprelength(2) espostlength(-999))
ochk exact_big   , opts(permexact)
ochk convmix     , opts(permabsp permrawp)
ochk rawp        , opts(permrawp)
ochk sideR       , opts(permside(right))
ochk lvl975      , opts(level(97.5) gteffects)
ochk probs_noP   , opts(permute(probabilities))

* ---- permute(probabilities) with a real table ----
qui levelsof unit, local(uu)
local pp
foreach u of local uu {
  local pp `pp' `u' : 0.3 0.3 0.4 ;
}
cap noi wooldid y unit t ttre , att cluster(unit) permute(probabilities) permprobs(`pp') permreps(30) rseed(777) dots(0)
dis "ODD|probs|RC=" _rc "|npossible=" e(perm_npossible)

* ---- REPLAY after a permutation run ----
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "REPLAY_EST_RC=" _rc "|perm_method=" "`e(perm_method)'"
matrix BPRE = e(b)
local ppre = e(perm_p_att)
cap noi wooldid
dis "REPLAY_RC=" _rc
dis "REPLAY_SAME_P=" (e(perm_p_att) == `ppre')

* ---- honestdid downstream, fast mode vs permanalyticse ----
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0)
cap noi honestdid, pre(1/2) post(3/4) parallel(0)
dis "HONEST_FAST_RC=" _rc
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0) permanalyticse
cap noi honestdid, pre(1/2) post(3/4) parallel(0)
dis "HONEST_ASE_RC=" _rc

* ---- does a permutation run leak state into the next plain run? ----
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2)
dis "AFTER_RC=" _rc "|perm_method=" "`e(perm_method)'" "|vcvstate=" "`e(wdid_vcvstate)'"
matrix BPOST = e(b)
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2)
matrix BPOST2 = e(b)
mata: printf("LEAKCHK maxdiff=%21x\n", max(abs(st_matrix("BPOST") - st_matrix("BPOST2"))))
dis "P09_DONE"
