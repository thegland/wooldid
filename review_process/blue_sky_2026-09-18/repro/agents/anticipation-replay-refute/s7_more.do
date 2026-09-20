clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 100
which wooldid

*=========== CASE P: permute(), 25 reps, cluster(uid) ======================
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
log using "P_live.txt", text replace
cap noi wooldid y uid t ttre, att cluster(uid) permute(shuffle) permreps(25) rseed(777) dots(0)
log close
dis "@@P rc=" _rc
log using "P_replay.txt", text replace
cap noi wooldid
log close
dis "@@Pr rc=" _rc

*=========== CASE C: continuous treatment, cohort-constant dose ============
clear all
use "dat.dta", clear
gen double dose = 0
replace dose = 1.5 if grp==1
replace dose = 3.0 if grp==2
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
log using "C_live.txt", text replace
cap noi wooldid y uid t ttre, att contreat(dose) contreatcontrols(basic) cluster(uid)
log close
dis "@@C rc=" _rc
log using "C_replay.txt", text replace
cap noi wooldid
log close
dis "@@Cr rc=" _rc
local emac : e(macros)
local esca : e(scalars)
dis "@@Cmac: `emac'"
dis "@@Csca: `esca'"

*=========== NOVEL: anticipation(k) vs an externally shifted ttre ==========
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui wooldid y uid t ttre, att att_i att_it att_itime robust anticipation(2)
mat A = e(b)
scalar aN = e(N)
qui wooldid y uid t ttre, att att_i att_it att_itime robust anticipation(0)
mat A0 = e(b)
clear all
use "dat.dta", clear
gen ttre2 = ttre - 2
qui wooldid y uid t ttre2, att att_i att_it att_itime robust
mat S = e(b)
scalar sN = e(N)
dis "@@AN antic_N=" aN " shift_N=" sN
local cn : colfullnames A
local sn : colfullnames S
dis "@@ANAMES antic=[`cn']"
dis "@@SNAMES shift=[`sn']"
local nc = colsof(A)
local worst = 0
forvalues j = 1/`nc' {
  local a = A[1,`j']
  local s = S[1,`j']
  local d = abs(`a'-`s')/max(abs(`s'),1e-12)
  if `d' > `worst' local worst = `d'
  dis "@@ACMP j=`j' a=" %21x `a' " s=" %21x `s' " rel=" %12.4e `d'
}
dis "@@AWORST reltol=" %12.4e `worst'
* sanity: anticipation(2) must differ from anticipation(0)
local d0 = abs(A[1,1]-A0[1,1])
dis "@@ASANITY |antic2 - antic0| = " %12.4e `d0'

*=========== NOVEL: anticipation with an event study ======================
clear all
use "dat.dta", clear
qui wooldid y uid t ttre, att robust anticipation(2) esprelength(2) espostlength(3)
mat EA = e(b)
local ean : colfullnames EA
clear all
use "dat.dta", clear
gen ttre2 = ttre - 2
qui wooldid y uid t ttre2, att robust esprelength(2) espostlength(3)
mat ES = e(b)
local esn : colfullnames ES
dis "@@ESN antic=[`ean']"
dis "@@ESN shift=[`esn']"
local nc = colsof(EA)
local w2 = 0
if "`ean'" == "`esn'" {
  forvalues j = 1/`nc' {
    local d = abs(EA[1,`j']-ES[1,`j'])/max(abs(ES[1,`j']),1e-12)
    if `d' > `w2' local w2 = `d'
  }
  dis "@@ESWORST reltol=" %12.4e `w2'
}
else dis "@@ESWORST NAME_MISMATCH"

dis "S7_DONE"
