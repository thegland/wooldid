* probe6.do -- `ties': establish what it changes and that it does so deterministically.
* Few clusters => a draw that reproduces the original sample exactly occurs with
* non-trivial probability, so replicate == observed ties really happen.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 666
local gN 48
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui gen long cl4 = mod(unit,4) + 1
save base6.dta, replace

local R  60
local SD 8642
cap erase bT0.dta
cap erase bT1.dta
cap erase bT1b.dta

cap prog drop dumpz
prog def dumpz
  syntax , tag(string)
  dis "`tag'|boot_ties=" e(boot_ties)
  dis "`tag'|boot_completed=" e(boot_completed)
  cap mat z0 = e(boot_z0)
  if _rc == 0 {
    local K = colsof(z0)
    forvalues k = 1/`K' {
      local nm : word `k' of `: colfullnames z0'
      dis "`tag'|z0|" "`nm'" "|" %21x z0[1,`k']
    }
  }
  cap mat ci = e(boot_ci_bc)
  if _rc == 0 mat list ci, format(%21x)
  tempname b v
  mat `b' = e(b)
  mat `v' = e(V)
  local K = colsof(`b')
  forvalues k = 1/`K' {
    local nm : word `k' of `: colfullnames `b''
    dis "`tag'|" "`nm'" "|b=" %21x `b'[1,`k'] "|se=" %21x sqrt(`v'[`k',`k'])
  }
end

dis "=== T0: no ties option ==="
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(cl4) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootdisplay(bc) bootsave(bT0) replace
dis "T0_RC=" _rc
dumpz , tag(T0)

dis "=== T1: ties ==="
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(cl4) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootdisplay(bc) ties bootsave(bT1) replace
dis "T1_RC=" _rc
dumpz , tag(T1)

dis "=== T1b: ties, repeated (determinism) ==="
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(cl4) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootdisplay(bc) ties bootsave(bT1b) replace
dis "T1b_RC=" _rc
dumpz , tag(T1b)

dis "=== hand oracle for z0 from the saved draws ==="
cap noi {
  use bT0.dta, clear
  qui sum main_att if _rep == 0, meanonly
  local bobs = r(mean)
  qui count if _rep > 0 & !missing(main_att)
  local RR = r(N)
  qui count if _rep > 0 & main_att < `bobs'
  local nlt = r(N)
  qui count if _rep > 0 & main_att == `bobs'
  local neq = r(N)
  dis "HAND: R=" `RR' " n_lt=" `nlt' " n_eq(exact ties with observed)=" `neq'
  dis "HAND: z0_noties = " %21x invnormal(`nlt'/`RR')
  dis "HAND: z0_ties   = " %21x invnormal((`nlt' + 0.5*`neq')/`RR')
}
dis "HAND_RC=" _rc

dis "=== T2: ties under boottstat (help says it touches nothing displayed) ==="
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(cl4) bootstrap(standard) boottstat reps(`R') rseed(`SD') dots(0)
dis "T2_RC=" _rc
dumpz , tag(T2NOTIES)
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(cl4) bootstrap(standard) boottstat reps(`R') rseed(`SD') dots(0) ties
dis "T2b_RC=" _rc
dumpz , tag(T2TIES)

dis "=== T3: ties on a continuous statistic must be a no-op (help) ==="
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootdisplay(bc)
dumpz , tag(T3NOTIES)
use base6.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootdisplay(bc) ties
dumpz , tag(T3TIES)

dis "PROBE6_DONE"
