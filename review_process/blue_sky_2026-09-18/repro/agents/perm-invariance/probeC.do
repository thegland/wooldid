*==========================================================================
* probeC.do -- permexact: the enumeration is exhaustive, so the reported
*              p-values MUST be invariant to cluster relabeling and row
*              order.  Plus RNG hygiene checks.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

*--- 12 counties in 6 states: 2 adopt 2004, 2 adopt 2006, 2 never.
*    6!/(2!2!2!) = 90 distinct assignments -> permexact enumerates them all.
prog def mksmall
  clear
  set seed 4242
  set sortseed 4242
  qui set obs 12
  qui gen long unit  = _n
  qui gen long state = ceil(unit/2)
  qui gen double ai  = rnormal()
  qui gen double ttre = .
  qui replace ttre = 2004 if state <= 2
  qui replace ttre = 2006 if state == 3 | state == 4
  qui expand 10
  qui bysort unit: gen int t = 1999 + _n
  qui gen double eps  = rnormal()*0.5
  qui gen byte   post = (t >= ttre) & !missing(ttre)
  qui gen double y    = ai + 0.2*(t-2000) + 0.8*post + eps
  qui drop eps
  qui order unit state t ttre y
end

local OPTX att att_i cluster(state) permute(shuffle) permexact rseed(20260917) dots(0) replace

mksmall
save panelC.dta, replace

cap noi wooldid y unit t ttre , `OPTX' permsave(ps_C0)
dis "RC_C0=" _rc
dis "REPS_C0=" e(perm_reps)
dmp C0
asgdmp C0

*--- row shuffle
use panelC.dta, clear
qui gen double _rr = mod(_n*2654435761, 1000003)
sort _rr
qui drop _rr
cap noi wooldid y unit t ttre , `OPTX' permsave(ps_C1rowsh)
dis "RC_C1rowsh=" _rc
cap dmp C1rowsh
cap asgdmp C1rowsh

*--- cluster relabel, order preserving
use panelC.dta, clear
qui replace state = state*13 + 5
cap noi wooldid y unit t ttre , `OPTX' permsave(ps_C2crelU)
dis "RC_C2crelU=" _rc
cap dmp C2crelU
cap asgdmp C2crelU

*--- cluster relabel, order REVERSING
use panelC.dta, clear
qui replace state = 100 - state
cap noi wooldid y unit t ttre , `OPTX' permsave(ps_C3crelD)
dis "RC_C3crelD=" _rc
cap dmp C3crelD
cap asgdmp C3crelD

*--- i relabel, order REVERSING
use panelC.dta, clear
qui replace unit = 500 - unit
cap noi wooldid y unit t ttre , `OPTX' permsave(ps_C4irelD)
dis "RC_C4irelD=" _rc
cap dmp C4irelD
cap asgdmp C4irelD

*======================= RNG HYGIENE =======================
use panelC.dta, clear
local rngA = c(rngstate)
cap noi wooldid y unit t ttre , `OPTX'
dis "RC_H1=" _rc
local rngB = c(rngstate)
dis "HYG_EXACT_STATE_UNCHANGED=" ("`rngA'" == "`rngB'")

*--- sampled draws with an rseed: does the state move?  (documented: yes)
use panelC.dta, clear
local rngC = c(rngstate)
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "RC_H2=" _rc
local rngD = c(rngstate)
dis "HYG_SAMPLED_STATE_UNCHANGED=" ("`rngC'" == "`rngD'")
dis "HYG_PERMSEED=" e(perm_seed)

*--- the documented replay claim (sthlp:542): a run WITHOUT rseed() can be
*    replayed exactly by `set rngstate e(perm_rngstate)' before repeating it
use panelC.dta, clear
set seed 31415
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(30) dots(0)
dis "RC_H3=" _rc
local rng0 = e(perm_rngstate)
dmp H3a
set rngstate `rng0'
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(30) dots(0)
dis "RC_H3b=" _rc
dmp H3b

*--- set seed honoured: same user seed before two identical no-rseed runs
use panelC.dta, clear
set seed 99999
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(30) dots(0)
dmp H4a
set seed 99999
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(30) dots(0)
dmp H4b

*--- does a permute() run with rseed() leave a LATER unrelated draw predictable?
use panelC.dta, clear
set seed 12345
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permreps(25) rseed(555) dots(0)
qui gen double zz1 = runiform() in 1
local z1 = zz1[1]
qui drop zz1
use panelC.dta, clear
set seed 67890
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permreps(25) rseed(555) dots(0)
qui gen double zz2 = runiform() in 1
local z2 = zz2[1]
dis "HYG_POSTDRAW_EQUAL=" (`z1' == `z2') " z1=" %21x `z1' " z2=" %21x `z2'

dis "PROBEC_DONE"
