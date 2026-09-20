*=====================================================================
* pA: skeptic repro of FUZZ-01 (elasticity header in binary mode),
*     FUZZ-03 (level 99.9/99.99), FUZZ-04 (estimates replay),
*     FUZZ-06 (dots(-5) on the two engines).
* SELF-CONTAINED: generates its own data.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz-refute/gen_mypanel.do"
local mN 120
local mseed 4242
include "`GEN'"
tempfile base
qui save "`base'", replace
qui count
dis "NOBS=" r(N)
tab ttre, m

*---------------------------------------------------------------------
* FUZZ-01 : binary mode + contreatelasticitytype()
*---------------------------------------------------------------------
dis "=================== F01 BINARY MODE ==================="
cap noi wooldid y i t ttre , att robust
dis "F01_BASE_RC=" _rc
scalar b0 = e(b)[1,1]
matrix Vb = e(V)
scalar s0 = sqrt(Vb[1,1])
dis "F01_BASE_ATT=" %21x b0 "  SE=" %21x s0
dis "F01_BASE_cetype=[" e(wdid_cetype) "] cetyperaw=[" e(wdid_cetyperaw) "]"

foreach ct in dydx dyex eydx eyex {
  dis "---- contreatelasticitytype(`ct') ----"
  cap noi wooldid y i t ttre , att robust contreatelasticitytype(`ct')
  dis "F01_`ct'_RC=" _rc
  if _rc == 0 {
    scalar bx = e(b)[1,1]
    matrix Vx = e(V)
    scalar sx = sqrt(Vx[1,1])
    dis "F01_`ct'_ATT=" %21x bx "  SE=" %21x sx
    dis "F01_`ct'_SAME=" (bx==b0) (sx==s0)
    dis "F01_`ct'_ecet=[" e(wdid_cetype) "] eraw=[" e(wdid_cetyperaw) "]"
  }
}
* domain screen without a parent?
cap noi wooldid y i t ttre , att robust contreatelasticitytype(garbage)
dis "F01_GARBAGE_RC=" _rc
* does the replay carry the same header?
cap noi wooldid y i t ttre , att robust contreatelasticitytype(eyex)
dis "F01_REPLAYSETUP_RC=" _rc
cap noi wooldid
dis "F01_REPLAY_RC=" _rc
* comparison: what do the OTHER contreat companions do in the same slot?
foreach o in "lattice(all 3)" "contreatcontrols(basic)" "contreatwithin" "makecfxplots" "cfxlattice(dose)" "latticeignoreweights" "ccnoabsorb" "contreatpoly(2)" "cfxplottypes(0)" "makecfxplotsbysg" "cfxplotnose" "semielasticity" {
  cap noi wooldid y i t ttre , att robust `o'
  dis "F01_COMPANION_RC [`o'] = " _rc
}

*---------------------------------------------------------------------
* FUZZ-03 : level()
*---------------------------------------------------------------------
dis "=================== F03 LEVEL ==================="
use "`base'", clear
foreach L in 95 90 99 99.5 97.5 90.25 99.95 68.27 99.9 99.99 99.8 99.7 99.75 84.13 95.4 92.5 {
  cap noi wooldid y i t ttre , att robust level(`L')
  dis "F03_LEVEL_`L'_RC=" _rc "  e(level)=" e(level) " c(level)=" c(level)
}
* workaround check: set level externally and omit level()
set level 99.9
cap noi wooldid y i t ttre , att robust
dis "F03_EXTERNAL999_RC=" _rc "  e(level)=" e(level)
mat lr = r(table)
set level 95
* level under the engines
cap noi wooldid y i t ttre , att robust bootstrap(standard) reps(12) rseed(5) level(99.9) dots(0)
dis "F03_BOOT999_RC=" _rc

*---------------------------------------------------------------------
* FUZZ-04 : estimates replay
*---------------------------------------------------------------------
dis "=================== F04 REPLAY ==================="
use "`base'", clear
cap noi wooldid y i t ttre , att robust
dis "F04_EST_RC=" _rc
cap noi estimates store m1
dis "F04_STORE_RC=" _rc
cap noi estimates replay m1
dis "F04_REPLAY_RC=" _rc
cap noi estimates table m1
dis "F04_TABLE_RC=" _rc
cap noi estimates restore m1
dis "F04_RESTORE_RC=" _rc
cap noi wooldid
dis "F04_BARE_RC=" _rc
cap noi wooldid ,
dis "F04_COMMA_RC=" _rc
cap noi wooldid , level(95)
dis "F04_COMMALEVEL_RC=" _rc
* does a stored-but-not-current replay work at all for other e-class cmds?
qui reg y x1
cap noi estimates store r1
cap noi estimates replay r1
dis "F04_REGREPLAY_RC=" _rc

*---------------------------------------------------------------------
* FUZZ-06 : dots(-5)
*---------------------------------------------------------------------
dis "=================== F06 DOTS ==================="
use "`base'", clear
cap noi wooldid y i t ttre , att robust permute(shuffle) permreps(20) rseed(9) dots(-5)
dis "F06_PERM_NEG_RC=" _rc
cap noi wooldid y i t ttre , att robust permute(shuffle) permreps(20) rseed(9) dots(0)
dis "F06_PERM_ZERO_RC=" _rc
cap noi wooldid y i t ttre , att robust bootstrap(standard) reps(20) rseed(9) dots(-5)
dis "F06_BOOT_NEG_RC=" _rc
* rseed is claimed to be validated on both; check the analogous claim
cap noi wooldid y i t ttre , att robust permute(shuffle) permreps(20) rseed(notanumber) dots(0)
dis "F06_PERM_BADSEED_RC=" _rc
dis "=================== DONE pA ==================="
