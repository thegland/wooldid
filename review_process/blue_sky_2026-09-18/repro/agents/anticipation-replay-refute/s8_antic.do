clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 100
which wooldid

*=========== NOVEL A: anticipation(k) vs an externally shifted ttre ========
use "dat.dta", clear
qui wooldid y uid t ttre, att att_i att_it att_itime robust anticipation(2)
mat A = e(b)
global aN = e(N)
qui wooldid y uid t ttre, att att_i att_it att_itime robust anticipation(0)
mat A0 = e(b)
use "dat.dta", clear
gen ttre2 = ttre - 2
qui wooldid y uid t ttre2, att att_i att_it att_itime robust
mat S = e(b)
global sN = e(N)
dis "@@AN antic_N=$aN shift_N=$sN"
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
local d0 = abs(A[1,1]-A0[1,1])
dis "@@ASANITY |antic2-antic0| = " %12.4e `d0'

*=========== NOVEL B: same with an event study =============================
use "dat.dta", clear
qui wooldid y uid t ttre, att robust anticipation(2) esprelength(2) espostlength(3)
mat EA = e(b)
local ean : colfullnames EA
use "dat.dta", clear
gen ttre2 = ttre - 2
qui wooldid y uid t ttre2, att robust esprelength(2) espostlength(3)
mat ES = e(b)
local esn : colfullnames ES
dis "@@ESN antic=[`ean']"
dis "@@ESN shift=[`esn']"
local w2 = 0
if "`ean'" == "`esn'" {
  local nc = colsof(EA)
  forvalues j = 1/`nc' {
    local d = abs(EA[1,`j']-ES[1,`j'])/max(abs(ES[1,`j']),1e-12)
    if `d' > `w2' local w2 = `d'
  }
  dis "@@ESWORST reltol=" %12.4e `w2'
}
else dis "@@ESWORST NAME_MISMATCH"

*=========== NOVEL C: NEGATIVE anticipation (delayed onset) ================
use "dat.dta", clear
cap noi wooldid y uid t ttre, att robust anticipation(-2)
dis "@@NEG rc=" _rc
if _rc == 0 {
  mat NA = e(b)
  global nN = e(N)
  use "dat.dta", clear
  gen ttre3 = ttre + 2
  qui wooldid y uid t ttre3, att robust
  mat NS = e(b)
  dis "@@NEGN antic_N=$nN shift_N=" e(N)
  local d = abs(NA[1,1]-NS[1,1])/max(abs(NS[1,1]),1e-12)
  dis "@@NEGCMP a=" %21x NA[1,1] " s=" %21x NS[1,1] " rel=" %12.4e `d'
}

*=========== NOVEL D: does level() actually move the gteffects CIs? ========
use "dat.dta", clear
qui wooldid y uid t ttre, att robust gteffects level(95)
mat G95 = e(gteffects)
scalar t95 = invttail(e(df_r), .025)
scalar df95 = e(df_r)
qui wooldid y uid t ttre, att robust gteffects level(80)
mat G80 = e(gteffects)
scalar t80 = invttail(e(df_r), .10)
local gn : colfullnames G95
dis "@@GN [`gn']"
dis "@@GDF df=" df95 " t95=" t95 " t80=" t80 " ratio=" t80/t95
local nr = rowsof(G95)
local ncg = colsof(G95)
dis "@@GDIM r=`nr' c=`ncg'"
* print the first 3 rows of both
forvalues r = 1/3 {
  local line95
  local line80
  forvalues c = 1/`ncg' {
    local line95 "`line95' " + string(G95[`r',`c'],"%12.6g")
    local line80 "`line80' " + string(G80[`r',`c'],"%12.6g")
  }
}
mat list G95, format(%12.6g) nohalf
mat list G80, format(%12.6g) nohalf

*=========== NOVEL E: does level() move the ES interval / e(b_es)? =========
use "dat.dta", clear
qui wooldid y uid t ttre, att robust esprelength(2) espostlength(3) level(95)
mat V95 = e(V_es_att)
mat M95 = e(wdidmainresults)
qui wooldid y uid t ttre, att robust esprelength(2) espostlength(3) level(80)
mat M80 = e(wdidmainresults)
local mn : colfullnames M95
dis "@@MN [`mn']"
dis "@@M95:"
mat list M95, format(%12.6g) nohalf
dis "@@M80:"
mat list M80, format(%12.6g) nohalf

dis "S8_DONE"
