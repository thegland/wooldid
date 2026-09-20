* probe11: settle the permute()+permtstat 4.672e-2 gap from probe5.  Under permute() the
* sample, N and K are constant across replications, so an aweight rescales every t by the
* same factor and the RANKS of |t| are invariant -> the permutation p must be IDENTICAL
* between [aw=k] and physical duplication even though the analytic SE is not.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

use "p2_A.dta", clear
qui wooldid y i t ttre [aw=k], att att_i cluster(i) permute(shuffle) permreps(60) rseed(909) permtstat
dis "PERMTS_A:"
matlist e(wdidmainresults), format(%16.10g)
dis "A_perm_p_att=" %21x e(perm_p_att)
dis "A_perm_p_att_i=" %21x e(perm_p_att_i)

use "p2_B.dta", clear
qui wooldid y i t ttre, att att_i cluster(i) permute(shuffle) permreps(60) rseed(909) permtstat
dis "PERMTS_B:"
matlist e(wdidmainresults), format(%16.10g)
dis "B_perm_p_att=" %21x e(perm_p_att)
dis "B_perm_p_att_i=" %21x e(perm_p_att_i)
dis "PROBE11_DONE"
