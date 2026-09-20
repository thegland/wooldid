* probe10.do -- SELF-FOUND reproducibility probes that rseed() alone does not cover:
*  (a) does the bootstrap depend on c(sortseed) at a fixed rseed()?
*  (b) what does the command leave the session RNG state at, and is that deterministic?
*  (c) rseed() vs no-rseed repetition semantics.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 1010
local gN 50
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
save base10.dta, replace

cap prog drop dumpse
prog def dumpse
  syntax , tag(string)
  tempname b v
  mat `b' = e(b)
  mat `v' = e(V)
  local K = colsof(`b')
  forvalues k = 1/`K' {
    local nm : word `k' of `: colfullnames `b''
    dis "`tag'|" "`nm'" "|b=" %21x `b'[1,`k'] "|se=" %21x sqrt(`v'[`k',`k'])
  }
  dis "`tag'|rng_after=" substr(c(rngstate),1,40)
end

local SD 4711
local R 30

dis "=== (a) sortseed dependence, COUNT leg ==="
use base10.dta, clear
set sortseed 1
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "A1_RC=" _rc
dumpse , tag(A1_SORTSEED1_COUNT)
use base10.dta, clear
set sortseed 987654321
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "A2_RC=" _rc
dumpse , tag(A2_SORTSEED987_COUNT)

dis "=== (a2) sortseed dependence, REID leg ==="
use base10.dta, clear
set sortseed 1
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "A3_RC=" _rc
dumpse , tag(A3_SORTSEED1_REID)
use base10.dta, clear
set sortseed 987654321
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "A4_RC=" _rc
dumpse , tag(A4_SORTSEED987_REID)

dis "=== (a3) sortseed dependence, UNCLUSTERED REID leg (plain bsample) ==="
use base10.dta, clear
set sortseed 1
cap noi wooldid ync unit t ttre, att att_i bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dumpse , tag(A5_SORTSEED1_UREID)
use base10.dta, clear
set sortseed 987654321
cap noi wooldid ync unit t ttre, att att_i bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dumpse , tag(A6_SORTSEED987_UREID)

dis "=== (c) no-rseed repetition: second run must differ; with rseed it must not ==="
set sortseed 20260916
use base10.dta, clear
set seed 55555
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') dots(0)
dumpse , tag(C1_NOSEED_FIRST)
use base10.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') dots(0)
dumpse , tag(C2_NOSEED_SECOND)
use base10.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dumpse , tag(C3_RSEED_FIRST)
use base10.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dumpse , tag(C4_RSEED_SECOND)
dis "PROBE10_DONE"
