* probe7.do -- SELF-FOUND probes:
*  (a) ROW-ORDER invariance of the bootstrap draws.  The count leg draws on a map frame
*      it SORTS by cluster (ado:8950); the re-identification leg draws on a copy of the
*      user's data IN ITS OWN ROW ORDER (ado:9341).  Does the same command + same rseed()
*      give the same answer when the user's rows are in a different order?
*  (b) NESTED CONSISTENCY of reps(): reps(20) vs the first 20 of reps(40) at one seed.
*  (c) the KEEP-i regime of the re-identification leg (cluster finer than i) -- T2/T3 code
*      shipped 2026-09-16 -- count vs reid agreement on att.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 777
local gN 60
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
save base7.dta, replace
* a shuffled copy: identical DATA, different ROW ORDER
use base7.dta, clear
set seed 99
qui gen double shufkey = runiform()
sort shufkey
drop shufkey
save base7shuf.dta, replace

local R  30
local SD 3141
foreach f in s1 s2 s3 s4 n20 n40 k1 k2 {
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

dis "=============== (a) ROW-ORDER invariance ==============="
use base7.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(s1) replace
dis "S1_RC=" _rc
dumpse , tag(S1_COUNT_SORTED)
use base7shuf.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(s2) replace
dis "S2_RC=" _rc
dumpse , tag(S2_COUNT_SHUFFLED)
use base7.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(s3) replace
dis "S3_RC=" _rc
dumpse , tag(S3_REID_SORTED)
use base7shuf.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(s4) replace
dis "S4_RC=" _rc
dumpse , tag(S4_REID_SHUFFLED)

dis "=============== (b) NESTED CONSISTENCY of reps() ==============="
use base7.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(`SD') dots(0) bootsave(n20) replace
dis "N20_RC=" _rc
use base7.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(40) rseed(`SD') dots(0) bootsave(n40) replace
dis "N40_RC=" _rc

dis "=============== (c) KEEP-i regime: i = cohort, cluster = unit (finer) ==============="
use base7.dta, clear
cap noi wooldid ync cohort t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(k1) replace
dis "K1_RC=" _rc
dumpse , tag(K1_COUNT_FINERCLUS)
use base7.dta, clear
cap noi wooldid ync cohort t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(k2) replace
dis "K2_RC=" _rc
dumpse , tag(K2_REID_FINERCLUS)

dis "=============== comparisons ==============="
cap prog drop cmp2
prog def cmp2
  syntax , f1(string) f2(string) tag(string) [maxrep(integer 0)]
  qui use "`f1'.dta", clear
  if `maxrep' > 0 qui keep if _rep <= `maxrep'
  qui ds
  qui rename * X1_*
  qui rename X1__rep _rep
  tempfile t1
  qui save `t1'
  qui use "`f2'.dta", clear
  if `maxrep' > 0 qui keep if _rep <= `maxrep'
  qui ds
  local v2 `r(varlist)'
  qui merge 1:1 _rep using `t1', nogen
  local worst = 0
  local worstv
  local nw = 0
  local wworst = 0
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
    if substr("`v'",1,2) == "w_" {
      local ++nw
      qui gen double __e = abs(`v' - X1_`v')
      qui sum __e, meanonly
      if r(max) > `wworst' local wworst = r(max)
      qui drop __e
    }
    qui drop __d
  }
  dis "`tag': MAXRELDIF = " %12.4e `worst' " on " "`worstv'"
  dis "`tag': W-DRAW cols = " `nw' " MAXABSDIFF = " %12.4e `wworst'
end

cap noi cmp2 , f1(s1) f2(s2) tag("(a1) COUNT leg: sorted vs shuffled rows")
dis "C1_RC=" _rc
cap noi cmp2 , f1(s3) f2(s4) tag("(a2) REID  leg: sorted vs shuffled rows")
dis "C2_RC=" _rc
cap noi cmp2 , f1(s1) f2(s3) tag("(a3) sorted:   COUNT vs REID")
dis "C3_RC=" _rc
cap noi cmp2 , f1(s2) f2(s4) tag("(a4) shuffled: COUNT vs REID")
dis "C4_RC=" _rc
cap noi cmp2 , f1(n20) f2(n40) tag("(b) reps(20) vs first 20 of reps(40)") maxrep(20)
dis "C5_RC=" _rc
cap noi cmp2 , f1(k1) f2(k2) tag("(c) keep-i regime: COUNT vs REID")
dis "C6_RC=" _rc
dis "PROBE7_DONE"
