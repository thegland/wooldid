* p5_perm_refusals.do -- do rows wooldid itself drops cause HARD permute() refusals,
* and do they change the EXACT permutation p-value?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootsample"

cap prog drop mkpanel
prog def mkpanel
  syntax , nstates(integer) ntreat(integer)
  clear
  set obs `nstates'
  gen int state = _n
  gen double sfx = rnormal()*0.4
  gen int onset = .
  replace onset = 2003 + mod(_n-1,5) if _n <= `ntreat'
  expand 2
  bysort state: gen int unit = (state-1)*2 + _n
  gen double ufx = rnormal()*0.3
  expand 10
  bysort unit: gen int year = 1999 + _n
  gen byte post = (onset < . & year >= onset)
  gen double y = 1 + .10*(year-1999) + sfx + ufx + .5*post + rnormal()*0.6
end

*=========================================== A: onset heterogeneity among DEAD rows only
mkpanel , nstates(20) ntreat(8)
qui su year
dis "A_YEARMIN=" r(min) " YEARMAX=" r(max)
* one stale record for a unit in state 7: a different onset, and no outcome
gen byte stale = (state == 7 & unit == 13 & year == r(min))
qui count if stale
dis "A_STALEROWS=" r(N)
replace onset = 2009 if stale
replace y = .        if stale
tempfile a0
save `a0', replace
dis "=== A1: FULL DATA, permute() ==="
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(50) rseed(1) dots(0)
dis "A1_RC=" _rc
dis "=== A2: the ONE stale row dropped (oracle) ==="
use `a0', clear
drop if stale
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(50) rseed(1) dots(0)
dis "A2_RC=" _rc
dis "A2_N=" e(N)
dis "=== A3: the ordinary (non-permute) run on the FULL data ==="
use `a0', clear
cap noi wooldid y unit year onset, cluster(state)
dis "A3_RC=" _rc
dis "A3_N=" e(N)

*=========================================== B: extended missing on a dead row
use `a0', clear
replace onset = .  if stale
replace onset = .a if stale
dis "=== B1: FULL DATA, one dead row with .a in the onset ==="
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(50) rseed(1) dots(0)
dis "B1_RC=" _rc
dis "=== B2: that row dropped ==="
drop if stale
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(50) rseed(1) dots(0)
dis "B2_RC=" _rc

*=========================================== C: permonsetcohorts bijection (T1, 1 day old)
clear
set seed 77
set obs 5
gen int coh = _n
gen int onset = 2003 + _n - 1
replace onset = . if coh == 5
expand 6
bysort coh: gen int clus = (coh-1)*6 + _n
expand 10
bysort clus: gen int year = 1999 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .1*(year-1999) + rnormal()*0.5 + .5*post
tempfile c0
save `c0', replace
dis "=== C1: clean permonsetcohorts run ==="
cap noi wooldid y coh year onset, cluster(clus) permute(shuffle) permreps(50) permonsetcohorts rseed(2) dots(0)
dis "C1_RC=" _rc
dis "C1_N=" e(N)
use `c0', clear
gen byte stale = (clus == 1 & year == 2000)
qui count if stale
dis "C_STALEROWS=" r(N)
replace onset = 2008 if stale
replace y = .        if stale
tempfile c1
save `c1', replace
dis "=== C2: plus ONE dead row carrying a stray onset ==="
cap noi wooldid y coh year onset, cluster(clus) permute(shuffle) permreps(50) permonsetcohorts rseed(2) dots(0)
dis "C2_RC=" _rc
dis "=== C3: that dead row dropped ==="
use `c1', clear
drop if stale
cap noi wooldid y coh year onset, cluster(clus) permute(shuffle) permreps(50) permonsetcohorts rseed(2) dots(0)
dis "C3_RC=" _rc
dis "C3_N=" e(N)
dis "=== C4: the ordinary run on the FULL data ==="
use `c1', clear
cap noi wooldid y coh year onset, cluster(clus)
dis "C4_RC=" _rc
dis "C4_N=" e(N)

*=========================================== D: permexact with dead clusters
clear
set seed 313
set obs 7
gen int unit = _n
gen double ufx = rnormal()*0.3
gen int onset = .
replace onset = 2003 if unit == 1
replace onset = 2004 if unit == 2
expand 12
bysort unit: gen int year = 1997 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .08*(year-1997) + ufx + .6*post + rnormal()*0.5
tempfile d0
save `d0', replace
* units 1 and 7 report no outcome at all
replace y = . if inlist(unit,1,7)
dis "=== D1: FULL DATA, permexact ==="
cap noi wooldid y unit year onset, cluster(unit) permute(shuffle) permexact rseed(3) dots(0)
dis "D1_RC=" _rc
dis "D1_N=" e(N)
cap mat li e(perm_p)
dis "=== D2: dead rows dropped (oracle), permexact ==="
use `d0', clear
drop if inlist(unit,1,7)
cap noi wooldid y unit year onset, cluster(unit) permute(shuffle) permexact rseed(3) dots(0)
dis "D2_RC=" _rc
dis "D2_N=" e(N)
cap mat li e(perm_p)

dis "P5_ALLDONE"
