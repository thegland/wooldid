*==========================================================================
* probeI.do -- permute(probabilities): does each cluster actually receive
*              the onset distribution the permprobs() table assigns to IT?
*   The engine pairs the probability matrix with the onset vector
*   POSITIONALLY (mata wooldid_permxp), relying on both being in ascending
*   cluster order.  Probe it three ways: table in ascending order, the same
*   table written in DESCENDING order, and the cluster ids relabelled so
*   that ascending order reverses.
*   Eight clusters are given a DETERMINISTIC target onset that is not the
*   one they actually have, so a misalignment is visible in permsave().
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

mkpanel
save panelI.dta, replace

* target onsets, columns ordered (2004 2006 2008 .)
*  s1->2008  s2->.  s3->2006  s4->2004  s5->.  s6->2008  s7->2006  s8->2004
*  s9,s10 -> half 2004, half never
local PASC (ttre 2004 2006 2008 .) (1 0 0 1 0) (2 0 0 0 1) (3 0 1 0 0) (4 1 0 0 0) (5 0 0 0 1) (6 0 0 1 0) (7 0 1 0 0) (8 1 0 0 0) (9 .5 0 0 .5) (10 .5 0 0 .5)
local PDESC (ttre 2004 2006 2008 .) (10 .5 0 0 .5) (9 .5 0 0 .5) (8 1 0 0 0) (7 0 1 0 0) (6 0 0 1 0) (5 0 0 0 1) (4 1 0 0 0) (3 0 1 0 0) (2 0 0 0 1) (1 0 0 1 0)

*--- I0 ascending
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(probabilities) permprobs(`PASC') permreps(24) rseed(20260917) dots(0) replace permsave(ps_I0)
dis "RC_I0=" _rc
cap dmp I0
cap asgdmp I0

*--- I1 the SAME table written in descending cluster order
use panelI.dta, clear
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(probabilities) permprobs(`PDESC') permreps(24) rseed(20260917) dots(0) replace permsave(ps_I1desc)
dis "RC_I1desc=" _rc
cap dmp I1desc
cap asgdmp I1desc

*--- I2 cluster ids relabelled 100-state, table rewritten by value
use panelI.dta, clear
qui replace state = 100 - state
local PREL (ttre 2004 2006 2008 .) (99 0 0 1 0) (98 0 0 0 1) (97 0 1 0 0) (96 1 0 0 0) (95 0 0 0 1) (94 0 0 1 0) (93 0 1 0 0) (92 1 0 0 0) (91 .5 0 0 .5) (90 .5 0 0 .5)
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(probabilities) permprobs(`PREL') permreps(24) rseed(20260917) dots(0) replace permsave(ps_I2rel)
dis "RC_I2rel=" _rc
cap dmp I2rel
cap asgdmp I2rel

dis "PROBEI_DONE"
