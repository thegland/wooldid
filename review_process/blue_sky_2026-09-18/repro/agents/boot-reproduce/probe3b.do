* probe3.do -- full reproducibility: same command twice in a session; rngstate replay;
*              option no-ops that must not move numbers.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 333
local gN 50
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
save base3b.dta, replace

cap prog drop dumpe
prog def dumpe
  syntax , tag(string)
  tempname b v
  mat `b' = e(b)
  mat `v' = e(V)
  local K = colsof(`b')
  dis "`tag': K=" `K'
  forvalues k = 1/`K' {
    local nm : word `k' of `: colfullnames `b''
    dis "`tag'|" "`nm'" "|b=" %21x `b'[1,`k'] "|v=" %21x `v'[`k',`k']
  }
  dis "`tag'|rngstate=" e(boot_rngstate)
  dis "`tag'|seed=" e(boot_seed)
end

dis "=============== P3a: same command twice, same session, rseed() ==============="
use base3b.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(0)
dis "R1_RC=" _rc
dumpe , tag(RUN1)
local rng1 = e(boot_rngstate)

use base3b.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(0)
dis "R2_RC=" _rc
dumpe , tag(RUN2)

dis "=============== P3b: count leg twice ==============="
use base3b.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(0)
dumpe , tag(RUN3)
use base3b.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(0)
dumpe , tag(RUN4)

dis "=============== P3c: e(boot_rngstate) replay (documented at sthlp:542) ==============="
* run WITHOUT rseed after a plain set seed, record the state, then replay it
use base3b.dta, clear
set seed 31415
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) dots(0)
dumpe , tag(NOSEED1)
local st1 = e(boot_rngstate)

use base3b.dta, clear
set rngstate `st1'
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) dots(0)
dumpe , tag(REPLAY1)
local st2 = e(boot_rngstate)
dis "STATES_EQUAL=" ("`st1'" == "`st2'")

dis "=============== P3d: option no-ops must not move numbers ==============="
use base3b.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(5)
dumpe , tag(DOTS5)
use base3b.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(0) bootdisplay(wald)
dumpe , tag(WALD)
use base3b.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) rseed(7777) dots(0) bootdf(normal)
dumpe , tag(DFNORM)

dis "=============== P3e: rseed given vs set seed beforehand ==============="
use base3b.dta, clear
set seed 7777
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(25) dots(0)
dumpe , tag(SETSEED7777)

dis "PROBE3_DONE"
