* probe5.do -- SELF-FOUND differential: a CONSTANT aweight of 1 is a mathematical no-op
* but it switches wooldid's bootstrap engine from the fweight COUNT leg to the
* PHYSICAL-EXPANSION leg (ado:392).  Both legs must give the same answer.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 555
local gN 50
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui gen double one = 1
save base5.dta, replace

local R  30
local SD 1357
cap erase bX.dta
cap erase bY.dta
cap erase bZ.dta
cap erase bW.dta

cap prog drop dumpe
prog def dumpe
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

dis "=== X: boottstat, NO weight, no cluster  -> fweight COUNT leg ==="
use base5.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) boottstat reps(`R') rseed(`SD') dots(0) bootsave(bX) replace
dis "X_RC=" _rc
dumpe , tag(X)

dis "=== Y: boottstat, [aw=one] (one==1), no cluster -> PHYSICAL leg ==="
use base5.dta, clear
cap noi wooldid ync unit t ttre [aw=one], att bootstrap(standard) boottstat reps(`R') rseed(`SD') dots(0) bootsave(bY) replace
dis "Y_RC=" _rc
dumpe , tag(Y)

dis "=== W: no boottstat, NO weight, no cluster -> aweight COUNT leg ==="
use base5.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bW) replace
dis "W_RC=" _rc
dumpe , tag(W)

dis "=== Z: no boottstat, [aw=one], no cluster -> aweight COUNT leg ==="
use base5.dta, clear
cap noi wooldid ync unit t ttre [aw=one], att bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bZ) replace
dis "Z_RC=" _rc
dumpe , tag(Z)

cap prog drop cmp2
prog def cmp2
  syntax , f1(string) f2(string) tag(string)
  qui use "`f1'.dta", clear
  qui ds
  qui rename * X1_*
  qui rename X1__rep _rep
  tempfile t1
  qui save `t1'
  qui use "`f2'.dta", clear
  qui ds
  local v2 `r(varlist)'
  qui merge 1:1 _rep using `t1', nogen
  local worst = 0
  local worstv
  foreach v of local v2 {
    if "`v'" == "_rep" continue
    cap confirm variable X1_`v'
    if _rc continue
    qui gen double __d = reldif(`v', X1_`v')
    qui replace __d = 0 if missing(`v') & missing(X1_`v')
    qui sum __d, meanonly
    if r(max) > `worst' {
      local worst = r(max)
      local worstv `v'
    }
    qui drop __d
  }
  dis "`tag': MAXRELDIF = " %12.4e `worst' " on " "`worstv'"
  cap confirm variable main_att
  if _rc == 0 {
    qui gen double __d = reldif(main_att, X1_main_att)
    qui sum __d, meanonly
    dis "`tag': main_att draws MAXRELDIF = " %12.4e r(max)
    list _rep main_att X1_main_att in 1/8, noobs
  }
end
cap noi cmp2 , f1(bX) f2(bY) tag("X(fweight count leg) vs Y(physical leg, aw=1)")
dis "CMPXY_RC=" _rc
cap noi cmp2 , f1(bW) f2(bZ) tag("W(count,no wt) vs Z(count,aw=1)")
dis "CMPWZ_RC=" _rc
cap noi cmp2 , f1(bX) f2(bW) tag("X(boottstat count) vs W(plain count)")
dis "CMPXW_RC=" _rc
dis "PROBE5_DONE"
