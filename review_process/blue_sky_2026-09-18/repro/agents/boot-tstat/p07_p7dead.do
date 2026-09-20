clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid
local gG 20
local gUnits 3
include "`WD'/gen_small.do"
gen double ysup = y
replace ysup = . if inlist(state, 4, 9, 14)

*==================================================================
* T1: boottstat (re-identification leg + the P7 assertion) with three
*     clusters that contribute NO rows to e(sample).
*==================================================================
*==================================================================
* T2: joint studentized test p -- recompute from the dump
*==================================================================
dis "=== T2 joint studentized p ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(199) boottstat rseed(88) att esprelength(3) espostlength(3) jointtests bootsave(d7) replace dots(0)
dis "T2_RC=" _rc
ereturn list
local JP = e(joint_horiz)
dis "T2_JOINTP=" %21x `JP'
preserve
use d7, clear
ds jointW_*
qui sum jointW_att if _rep == 0
local WOBS = r(mean)
qui count if _rep > 0 & jointW_att < .
local NJ = r(N)
qui count if _rep > 0 & jointW_att < . & jointW_att >= `WOBS'
local GE = r(N)
local PMINE = (1 + `GE')/(`NJ' + 1)
dis "T2_WOBS=" %21x `WOBS' "  NJ=`NJ'  GE=`GE'"
dis "T2_PMINE=" %21x `PMINE'
dis "T2_PDIFF=" %21x (`PMINE' - `JP')
restore
dis "P07_DONE"
