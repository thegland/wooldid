clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103

which wooldid

* ---- my own panel, nothing borrowed ----
clear
set obs 80
gen long id = _n
gen byte clu = ceil(id/4)
gen double ui = rnormal()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2006 if mod(id,4)==1
replace ttre = 2008 if mod(id,4)==2
expand 10
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 2 + 0.3*(t-ttre) if t>=ttre & ttre<.

* ---- A. macro expansion of a real ----
foreach L in 90 95 99 97.5 99.5 99.9 90.1 99.99 68.3 92.5 91.1 93.7 {
  local x = `L'
  dis "EXPAND `L' -> |`x'|"
}

* ---- B. does regress / set level accept them? ----
foreach L in 90 90.1 91.1 92.5 93.7 95 97.5 99 99.5 99.9 99.99 68.3 {
  cap reg y i.t if id<=40, level(`L')
  local r1 = _rc
  cap set level `L'
  local r2 = _rc
  dis "LITERAL `L' regrc=`r1' setlevrc=`r2' clevel=" c(level)
  set level 95
}

* ---- B2. the SAME levels routed through a `real' local, as syntax does ----
foreach L in 90 90.1 91.1 92.5 93.7 95 97.5 99 99.5 99.9 99.99 68.3 {
  local x = `L'
  cap set level `x'
  dis "VIAREAL `L' -> |`x'| setlevrc=" _rc
  set level 95
}

* ---- C. exhaustive one-decimal scan of the MECHANISM (no wooldid; fast) ----
local badlist ""
forval k = 100/9999 {
  local L = `k'/100
  if `L' < 10 continue
  local x = `L'
  cap set level `x'
  if _rc != 0 local badlist "`badlist' `L'"
  set level 95
}
dis "MECH_BAD_COUNT=" wordcount("`badlist'")
dis "MECH_BAD: `badlist'"

* ---- D. confirm with wooldid itself on a handful ----
foreach L in 90 90.1 92.5 97.5 99 99.5 99.9 99.99 68.3 91.1 93.7 94.9 {
  cap wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) level(`L')
  dis "WDLEV `L' rc=" _rc " elevel=" e(level)
}
dis "GRID_DONE"

* ---- E. the no-option route after set level 99.9 ----
set level 95
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) level(92.5)
dis "PRE_RC=" _rc
dis "PRE_LEVEL=" e(level)
mat lb = e(wdidmainresults)
dis "PRE_COLS=" "`: colnames lb'"
cap noisily set level 99.9
dis "SETL999_RC=" _rc " clevel=" c(level)
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2)
dis "NOOPT_RC=" _rc
dis "POST_LEVEL=" e(level)
dis "POST_CMD=" e(cmd)
mat lb2 = e(wdidmainresults)
dis "POST_COLS=" "`: colnames lb2'"
cap noisily reg y i.t if id<=40
dis "REG_NOOPT_RC=" _rc " reglevel=" e(level)
set level 95
dis "ALL_DONE"
