* p18: replay after a permute() run, compared in SCALARS (p09 used a local, whose
* decimal round-trip loses precision - the same trap ado:7645 warns about).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) subgroups(sg2) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "RP|est|RC=" _rc "|perm_method=`e(perm_method)'"
scalar P0 = e(perm_p_att)
scalar S0 = e(perm_sd_att)
matrix B0  = e(b)
matrix PP0 = e(perm_p)
matrix CI0 = e(perm_ci)
matrix KK0 = e(perm_kernel)
cap noi wooldid
dis "RP|replay|RC=" _rc
scalar P1 = e(perm_p_att)
matrix B1  = e(b)
matrix PP1 = e(perm_p)
matrix CI1 = e(perm_ci)
matrix KK1 = e(perm_kernel)
dis "RP|p_att_pre=" %21x P0 "|p_att_post=" %21x P1 "|EQ=" (P0==P1)
mata:
 printf("RP|maxdiff b=%21x p=%21x ci=%21x kernel=%21x\n",
   max(abs(editmissing(st_matrix("B0"),0)  - editmissing(st_matrix("B1"),0))),
   max(abs(editmissing(st_matrix("PP0"),0) - editmissing(st_matrix("PP1"),0))),
   max(abs(editmissing(st_matrix("CI0"),0) - editmissing(st_matrix("CI1"),0))),
   max(abs(editmissing(st_matrix("KK0"),0) - editmissing(st_matrix("KK1"),0))))
end
* and the p09 poisson+aw refusal message
cap noi wooldid ycount unit t ttre [aw=w] , att cluster(unit) poisson esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "RP|pois_aw|RC=" _rc
dis "P18_DONE"
