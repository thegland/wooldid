* R7_rc0_hunt.do -- SELF-CONTAINED.
* Hunt for the claimed "rc 0 with a printed all-missing table and NO e(b)" branch of
* `semielasticity fastnose' on a sign-crossing outcome.  Eight independently drawn
* 40-unit x 8-year two-onset panels; record rc, whether e(b) posted, and e(cmd).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"

forvalues s = 1/8 {
  clear
  set seed `=33000 + `s'*137'
  set sortseed `=33000 + `s'*137'
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
  qui summ y2
  local mn = r(min)
  local mx = r(max)

  ereturn clear
  cap noi wooldid y2 u year onset, att cluster(u) semielasticity fastnose
  local rcf = _rc
  cap confirm matrix e(b)
  local hb = (_rc==0)
  local ec "`e(cmd)'"
  dis "HUNT s=`s' miny=`mn' maxy=`mx' fastnose_rc=`rcf' hasb=`hb' ecmd=`ec'"
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WH = e(wdidmainresults)
    matrix list WH
  }

  ereturn clear
  cap noi wooldid y2 u year onset, att cluster(u) semielasticity
  dis "HUNT s=`s' plain_rc=" _rc

  ereturn clear
  cap noi wooldid y2 u year onset, att cluster(u) semielasticity permute(shuffle) permreps(30) rseed(7)
  dis "HUNT s=`s' perm_rc=" _rc
}
dis "R7 done"
