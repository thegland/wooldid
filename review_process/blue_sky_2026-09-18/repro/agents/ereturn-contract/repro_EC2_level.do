* EC-2: level(99.9) is refused r(198) with a message that contradicts what was typed,
* and `set level 99.9' makes wooldid unusable with no level() option at all.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
replace ttre = 11 if grp==4
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + 0.25*t + 1.2*post + rnormal()

foreach L in 90 95 99 97.5 99.5 68.5 99.9 90.1 99.99 99.95 68.3 92.5 {
  cap noi regress y t, level(`L')
  local a = _rc
  cap set level `L'
  local b = _rc
  set level 95
  cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) level(`L')
  local c = _rc
  dis "level(`L'): regress rc=`a'  set level rc=`b'  WOOLDID rc=`c'"
}
set level 95
dis "--- and with no level() option at all, after set level 99.9 ---"
set level 99.9
dis "set level 99.9 rc=" _rc "  c(level)=" c(level)
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3)
dis "WOOLDID with NO level() option: rc=" _rc
set level 95

dis ""
dis "======== the mechanism ========"
foreach L in 90 97.5 99.5 68.3 99.95 92.5 99.9 90.1 99.99 {
  local x = `L'
  dis "  level(`L') : a Stata `real' expands as `x'"
}
dis "  -> the three that fail are exactly the three whose %18.0g macro expansion"
dis "     shows more than two decimals, which is what `set level' refuses"
