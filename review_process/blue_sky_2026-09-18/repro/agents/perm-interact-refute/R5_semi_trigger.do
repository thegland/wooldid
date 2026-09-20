* R5_semi_trigger.do -- SELF-CONTAINED.
* (a) LEVEL LADDER: shift the same outcome up by c = 0,1,...,10 and record where
*     semielasticity starts to work, against min(y) at that c.  If it still fails at a c
*     where min(y) > 0, "the outcome takes non-positive values" is the wrong trigger.
* (b) rc-0 HUNT: a 40-unit x 8-year, two-onset panel -- the shape on which the accused
*     agent saw `semielasticity fastnose' return rc 0 with a printed all-missing table.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

*---------------------------------------------------------------- panel 1 (60x12)
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
summ ycross
dis "LADDER_min_raw=" r(min) " max=" r(max)

forvalues c = 0/10 {
  cap drop yc
  qui gen double yc = ycross + `c'
  qui summ yc
  local mn = r(min)
  ereturn clear
  cap noi wooldid yc firm year onset, att cluster(firm) semielasticity
  local rc = _rc
  local a = .
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    local a = WR[1,1]
  }
  dis "LADDER c=`c' miny=`mn' rc=`rc' att=`a'"
  ereturn clear
  cap noi wooldid yc firm year onset, att cluster(firm) semielasticity fastnose
  local rcf = _rc
  cap confirm matrix e(b)
  local hb = (_rc==0)
  dis "LADDERFN c=`c' miny=`mn' rc=`rcf' hasb=`hb' ecmd=`e(cmd)'"
}

*---------------------------------------------------------------- panel 2 (40x8, 2 onsets)
clear
set seed 33117
set sortseed 33117
set obs 40
gen int u = _n
gen double a_i = rnormal(0, 1.1)
gen byte g2 = mod(_n, 3)
gen int onset = .
replace onset = 2003 if g2 == 1
replace onset = 2005 if g2 == 2
expand 8
bysort u: gen int year = 2000 + _n
gen double e2 = rnormal(0, 1.3)
gen int rel = year - onset
gen double eff = 0
replace eff = 0.6 + 0.4*rel if onset < . & rel >= 0
gen double y2 = a_i + .15*(year-2000) + eff + e2
drop rel eff
summ y2
dis "P2_min=" r(min) " max=" r(max)

ereturn clear
cap noi wooldid y2 u year onset, att cluster(u)
dis "P2_base_RC=" _rc
ereturn clear
cap noi wooldid y2 u year onset, att cluster(u) semielasticity fastnose
dis "P2_fn_RC=" _rc
cap confirm matrix e(b)
dis "P2_fn_hasb=" (_rc==0)
dis "P2_fn_ecmd=" "`e(cmd)'"
cap confirm matrix e(wdidmainresults)
if _rc == 0 {
  matrix WR2 = e(wdidmainresults)
  matrix list WR2
}
ereturn clear
cap noi wooldid y2 u year onset, att cluster(u) semielasticity
dis "P2_plain_RC=" _rc
ereturn clear
cap noi wooldid y2 u year onset, att cluster(u) semielasticity permute(shuffle) permreps(30) rseed(7)
dis "P2_perm_RC=" _rc
dis "---- R5 done ----"
