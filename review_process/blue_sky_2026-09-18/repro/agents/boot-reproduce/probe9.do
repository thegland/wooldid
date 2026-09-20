* probe9.do -- further differential tests of "two code paths, one object".
*  (i)   att_i ALONE vs att att_i : att_i's draws must not move.
*  (ii)  all four estimands vs att alone (clustered) : att's draws must not move.
*  (iii) e(bootdist) vs bootsave() agreement.
*  (iv)  an `if' restriction: count leg vs reid leg draws.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 999
local gN 50
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui gen byte keepme = (mod(unit,10) != 0)
save base9.dta, replace

local R 30
local SD 5150
foreach f in q1 q2 q3 q4 q5 q6 {
  cap erase `f'.dta
}

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
end

dis "=== (i) att_i alone vs att att_i, cluster(unit) ==="
use base9.dta, clear
cap noi wooldid ync unit t ttre, att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(q1) replace
dis "Q1_RC=" _rc
dumpse , tag(Q1_ATTI_ALONE)
use base9.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(q2) replace
dis "Q2_RC=" _rc
dumpse , tag(Q2_ATT_ATTI)

dis "=== (ii) all four estimands vs att alone, cluster(unit) ==="
use base9.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(q3) replace
dis "Q3_RC=" _rc
dumpse , tag(Q3_ATT_ALONE)
use base9.dta, clear
cap noi wooldid ync unit t ttre, att att_it att_itime att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(q4) replace
dis "Q4_RC=" _rc
dumpse , tag(Q4_ALLFOUR)

dis "=== (iii) e(bootdist) vs bootsave() ==="
use base9.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
cap mat bd = e(bootdist)
dis "BOOTDIST_ROWS=" rowsof(bd) " COLS=" colsof(bd)
dis "BOOTDIST_colnames: " "`: colfullnames bd'"
forvalues j = 1/5 {
  dis "bootdist row `j' = " %21x bd[`j',1]
}

dis "=== (iv) `if' restriction: count leg vs reid leg ==="
use base9.dta, clear
cap noi wooldid ync unit t ttre if keepme, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(q5) replace
dis "Q5_RC=" _rc
dumpse , tag(Q5_IF_COUNT)
use base9.dta, clear
cap noi wooldid ync unit t ttre if keepme, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(q6) replace
dis "Q6_RC=" _rc
dumpse , tag(Q6_IF_REID)

do cmplib.do
cap noi cmpf , f1(q1) f2(q2) tag("(i) att_i alone vs att+att_i")
dis "CI_RC=" _rc
cap noi cmpf , f1(q3) f2(q4) tag("(ii) att alone vs all four")
dis "CII_RC=" _rc
cap noi cmpf , f1(q5) f2(q6) tag("(iv) if-restricted: count vs reid")
dis "CIV_RC=" _rc
dis "PROBE9_DONE"
