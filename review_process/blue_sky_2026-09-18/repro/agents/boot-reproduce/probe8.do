* probe8.do -- UNCLUSTERED: do the two legs converge to the same att standard error?
* Count leg draws with `bsample, cluster(rowid) idcluster()` on the map frame; the
* re-identification leg draws with a plain `bsample` on the full frame.  Same seed,
* different draws.  Is the difference Monte Carlo, or a different estimator?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 888
local gN 40
local gT 7
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
save base8.dta, replace

local R 300
foreach f in u1 u2 u3 u4 {
  cap erase `f'.dta
}

cap prog drop seatt
prog def seatt
  syntax , tag(string)
  tempname b v
  mat `b' = e(b)
  mat `v' = e(V)
  local k = colnumb(`b',"main_att")
  dis "`tag'|b_att=" %21x `b'[1,`k'] "|se_att=" %21x sqrt(`v'[`k',`k']) "|se_dec=" %12.9f sqrt(`v'[`k',`k']) "|completed=" e(boot_completed) "|failed=" e(boot_failed)
end

dis "=== U1: att alone, NO cluster, seed A (COUNT leg) ==="
use base8.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) reps(`R') rseed(1001) dots(0) bootsave(u1) replace
dis "U1_RC=" _rc
seatt , tag(U1_COUNT_seedA)

dis "=== U2: att + att_i, NO cluster, seed A (REID leg) ==="
use base8.dta, clear
cap noi wooldid ync unit t ttre, att att_i bootstrap(standard) reps(`R') rseed(1001) dots(0) bootsave(u2) replace
dis "U2_RC=" _rc
seatt , tag(U2_REID_seedA)

dis "=== U3: att alone, NO cluster, seed B (COUNT leg) ==="
use base8.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) reps(`R') rseed(2002) dots(0) bootsave(u3) replace
dis "U3_RC=" _rc
seatt , tag(U3_COUNT_seedB)

dis "=== U4: att + att_i, NO cluster, seed B (REID leg) ==="
use base8.dta, clear
cap noi wooldid ync unit t ttre, att att_i bootstrap(standard) reps(`R') rseed(2002) dots(0) bootsave(u4) replace
dis "U4_RC=" _rc
seatt , tag(U4_REID_seedB)

dis "=== SDs computed straight from the saved draw files ==="
foreach f in u1 u2 u3 u4 {
  qui use "`f'.dta", clear
  qui sum main_att if _rep > 0
  dis "`f': n=" r(N) " mean=" %12.9f r(mean) " sd=" %12.9f r(sd)
  qui sum main_att if _rep == 0
  dis "`f': observed(rep0) = " %12.9f r(mean)
}
dis "PROBE8_DONE"
