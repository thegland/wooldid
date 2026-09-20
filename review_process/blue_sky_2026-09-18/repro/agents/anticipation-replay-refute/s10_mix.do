clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 120
which wooldid

*=========== A: anticipation(k) + esfixedbaseperiod vs shifted ttre ========
use "dat.dta", clear
cap noi qui wooldid y uid t ttre, att robust anticipation(2) esfixedbaseperiod esrelativeto(-1) esprelength(2) espostlength(3)
dis "@@FA rc=" _rc
mat FA = e(b)
global faN = e(N)
local fan : colfullnames FA
use "dat.dta", clear
gen ttre2 = ttre - 2
cap noi qui wooldid y uid t ttre2, att robust esfixedbaseperiod esrelativeto(-1) esprelength(2) espostlength(3)
dis "@@FS rc=" _rc
mat FS = e(b)
local fsn : colfullnames FS
dis "@@FN antic=[`fan']"
dis "@@FN shift=[`fsn']"
dis "@@FNOBS antic_N=$faN shift_N=" e(N)
if "`fan'" == "`fsn'" {
  local w = 0
  local nc = colsof(FA)
  forvalues j = 1/`nc' {
    local d = abs(FA[1,`j']-FS[1,`j'])/max(abs(FS[1,`j']),1e-12)
    if `d' > `w' local w = `d'
  }
  dis "@@FWORST reltol=" %12.4e `w'
}
else dis "@@FWORST NAME_MISMATCH"

*=========== B: anticipation() with cohort-level i (coarse) ================
use "dat.dta", clear
egen long ci = group(ttre), missing
cap noi qui wooldid y ci t ttre, att robust anticipation(2) fe(uid)
dis "@@CA rc=" _rc
mat CA = e(b)
use "dat.dta", clear
gen ttre2 = ttre - 2
egen long ci = group(ttre2), missing
cap noi qui wooldid y ci t ttre2, att robust fe(uid)
dis "@@CS rc=" _rc
mat CS = e(b)
local d = abs(CA[1,1]-CS[1,1])/max(abs(CS[1,1]),1e-12)
dis "@@CWORST a=" %21x CA[1,1] " s=" %21x CS[1,1] " rel=" %12.4e `d'

*=========== C: CS-2 check - unspecified level() under the engines =========
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 90
qui wooldid y uid t ttre, att robust bootstrap(standard) reps(25) rseed(101) dots(0)
mat MB = e(wdidmainresults)
local mbn : colfullnames MB
dis "@@CS2B [`mbn'] elevel=" e(level) " clevel_now=" c(level)
mat list MB, format(%12.6g)
set level 95
qui wooldid y uid t ttre, att robust bootstrap(standard) reps(25) rseed(101) dots(0)
mat MB95 = e(wdidmainresults)
local mb95n : colfullnames MB95
dis "@@CS2B95 [`mb95n'] elevel=" e(level)
mat list MB95, format(%12.6g)

*=========== D: replay with NO dataset in memory (sthlp promise) ===========
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui wooldid y uid t ttre, att robust esprelength(2) espostlength(3)
drop _all
cap noi wooldid
dis "@@NOD rc=" _rc " eN=" e(N)

dis "S10_DONE"
