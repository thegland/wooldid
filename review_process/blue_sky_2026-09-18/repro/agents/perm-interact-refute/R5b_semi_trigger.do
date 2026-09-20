* R5b_semi_trigger.do -- SELF-CONTAINED.
* (a) FRAME/DATA STATE after a wooldid run that dies inside margins (r459).
* (b) LEVEL LADDER: shift the same outcome up by c = 0..10; where does semielasticity
*     start to work, relative to min(y)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

set obs 60
gen int firm = _n
gen double fe_i = rnormal(0, 1.4)
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
gen double ycross = fe_i + fe_t + effect + eps
drop rel effect
gen double USERMARKER = 12345
forvalues c = 0/10 {
  gen double yc`c' = ycross + `c'
}
summ ycross
dis "LADDER_min_raw=" r(min) " max=" r(max)

*========== (a) frame / variable state after a failing run
qui frame
dis "FRAME_before=" r(currentframe)
qui ds
dis "VARS_before=" r(varlist)
cap noi wooldid ycross firm year onset, att cluster(firm) semielasticity
dis "FAIL_RC=" _rc
qui frame
dis "FRAME_after=" r(currentframe)
cap qui ds
dis "DS_after_rc=" _rc
cap dis "VARS_after=" r(varlist)
cap confirm variable ycross
dis "HAS_ycross_after=" (_rc==0)
cap confirm variable USERMARKER
dis "HAS_USERMARKER_after=" (_rc==0)
cap count
dis "N_after=" r(N)
frame dir

*========== same, after the fastnose variant
cap noi wooldid ycross firm year onset, att cluster(firm) semielasticity fastnose
dis "FAIL2_RC=" _rc
qui frame
dis "FRAME2_after=" r(currentframe)
cap confirm variable ycross
dis "HAS_ycross_after2=" (_rc==0)
frame dir

*========== (b) ladder
forvalues c = 0/10 {
  qui summ yc`c'
  local mn = r(min)
  ereturn clear
  cap noi wooldid yc`c' firm year onset, att cluster(firm) semielasticity
  local rc = _rc
  local a = .
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    local a = WR[1,1]
  }
  dis "LADDER c=`c' miny=`mn' rc=`rc' att=`a'"
  qui frame
  local fr = r(currentframe)
  ereturn clear
  cap noi wooldid yc`c' firm year onset, att cluster(firm) semielasticity fastnose
  local rcf = _rc
  cap confirm matrix e(b)
  local hb = (_rc==0)
  local ec "`e(cmd)'"
  dis "LADDERFN c=`c' miny=`mn' rc=`rcf' hasb=`hb' ecmd=`ec' frame=`fr'"
}
dis "R5 done"
