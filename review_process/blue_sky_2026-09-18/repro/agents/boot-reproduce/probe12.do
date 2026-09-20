* probe12.do -- more "two code paths, one object" differentials around bootstrap():
*  (a) the cell-collapse path vs nofastcollapse UNDER bootstrap(): same draws, same answer?
*  (b) point estimates must be unchanged by bootstrap() (sthlp:532).
*  (c) e(V_analytic) under bootstrap() must equal e(V) from the same run without bootstrap().
*  (d) bootsave() dataset characteristics vs the help (sthlp:542/557).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 1313
local gN 50
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
local gpercell 3
local gsub 2
include "`GEN'/gen_core.do"
save base12.dta, replace
qui count
dis "NOBS=" r(N)

cap prog drop dumpall
prog def dumpall
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
cap prog drop dumpanal
prog def dumpanal
  syntax , tag(string)
  tempname b v
  mat `b' = e(b)
  cap mat `v' = e(V_analytic)
  if _rc {
    dis "`tag'|NO_V_ANALYTIC"
    exit
  }
  local K = colsof(`b')
  forvalues k = 1/`K' {
    local nm : word `k' of `: colfullnames `b''
    dis "`tag'|" "`nm'" "|b=" %21x `b'[1,`k'] "|seANA=" %21x sqrt(`v'[`k',`k'])
  }
end

local R 30
local SD 8080

dis "=== (b,c) plain analytic run ==="
use base12.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) esprelength(2) espostlength(2)
dis "P_RC=" _rc
dumpall , tag(PLAIN)

dis "=== bootstrap run, default collapse ==="
use base12.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(c1) replace
dis "B1_RC=" _rc
dumpall , tag(BOOT_COLLAPSE)
dumpanal , tag(BOOT_COLLAPSE)

dis "=== bootstrap run, nofastcollapse ==="
use base12.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) bootstrap(standard) reps(`R') rseed(`SD') dots(0) nofastcollapse bootsave(c2) replace
dis "B2_RC=" _rc
dumpall , tag(BOOT_NOCOLLAPSE)
dumpanal , tag(BOOT_NOCOLLAPSE)

dis "=== plain analytic run, nofastcollapse ==="
use base12.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) nofastcollapse
dis "P2_RC=" _rc
dumpall , tag(PLAIN_NOCOLLAPSE)

dis "=== (d) bootsave characteristics ==="
use c1.dta, clear
char list _dta[]
qui ds
dis "c1 vars: " "`r(varlist)'"
dis "c1 N=" _N

do cmplib.do
cap noi cmpf , f1(c1) f2(c2) tag("collapse vs nofastcollapse draws")
dis "CC_RC=" _rc
dis "PROBE12_DONE"
