clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
which wooldid
include "dgp.do"
desc, short
tab ttre, m
timer clear 1
timer on 1
cap noi wooldid y unit t ttre, cluster(unit) att esprelength(3) espostlength(3) jointtests permute(shuffle) permreps(99) rseed(7) permsave(p00) replace
dis "RUN_RC=" _rc
timer off 1
timer list 1
ereturn list
