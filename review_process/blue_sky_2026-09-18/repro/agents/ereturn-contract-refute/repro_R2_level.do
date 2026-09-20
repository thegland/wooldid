*==============================================================================
* REPRO (skeptic re-derivation of EC-2): wooldid refuses 40 of the 90 one-decimal
* confidence levels between 90.0 and 99.9, all of which `regress' and `set level'
* accept; and a session-wide `set level 99.9' makes every wooldid call fail.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_R2_level.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103
which wooldid

* --- an ordinary staggered panel ---
clear
set obs 80
gen long id = _n
gen double ui = rnormal()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2006 if mod(id,4)==1
replace ttre = 2008 if mod(id,4)==2
expand 10
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 2 + 0.3*(t-ttre) if t>=ttre & ttre<.

* ---------- 1. EVERY one-decimal level in [90,100): regress vs wooldid ----------
local nbad 0
local badl ""
forval k = 900/999 {
  local L = `k'/10
  cap reg y i.t, level(`L')
  local rr = _rc
  cap set level `L'
  local rs = _rc
  set level 95
  cap wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) level(`L')
  local rw = _rc
  if `rw' != 0 {
    local ++nbad
    local badl "`badl' `L'"
    if `rr' == 0 & `rs' == 0 dis "MISMATCH level(`L'): regress rc=`rr', set level rc=`rs', wooldid rc=`rw'"
  }
}
dis "WOOLDID REFUSED `nbad' of 100 one-decimal levels in [90,100)"
dis "REFUSED:`badl'"

* ---------- 2. the mechanism ----------
foreach L in 95 99.9 94.9 92.5 {
  local x = `L'
  dis "level(`L') macro-expands to |`x'|"
}

* ---------- 3. a session-wide `set level' bricks wooldid entirely ----------
set level 95
cap wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) level(92.5)
dis "prior good run rc=" _rc " e(level)=" e(level)
mat A = e(wdidmainresults)
dis "  interval columns: `: colnames A'"
set level 99.9
dis "set level 99.9 rc=" _rc " -> c(level)=" c(level)
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2)
dis "wooldid with NO level() option rc=" _rc
dis "  e(level) still = " e(level) "  (the PREVIOUS run's level)"
mat B = e(wdidmainresults)
dis "  interval columns still: `: colnames B'"
cap reg y i.t
dis "regress in the same session rc=" _rc
set level 95
dis "REPRO_R2_DONE"
