* NOTE (added after p5): SECTIONS B AND C OF THIS FILE ARE INVALID and are superseded by
* p5_perm_refusals.do. Two bugs in THIS do-file, not in wooldid: (i) mkpanel numbers years
* from 2000, so the `year == 1999` targeting in B and C matched no rows (the log shows
* "0 real changes made"); (ii) section B reuses the tempfile name `c0`, clobbering the
* section C dataset. Section C's design was also wrong (two never-treated cohort levels),
* which permonsetcohorts correctly refuses. Sections A and D of this file are VALID and are
* the evidence for HV-PERM-01 and HV-VAL-01.

* p3_perm_assert.do
*  A) permute(): is the randomization design drawn over clusters that contribute
*     ZERO rows to e(sample)?  permsave() names every cluster it assigns to.
*  B) permute() HARD REFUSALS driven by rows wooldid itself would drop:
*     onset heterogeneity inside a cluster that exists only among dead rows.
*  C) permonsetcohorts (T1, one day old) bijection check on the same unscreened sample.
*  D) wooldid_validate's parse-time asserts evaluate on the RAW user `if'.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootsample"

cap prog drop mkpanel
prog def mkpanel
  syntax , nstates(integer) ntreat(integer) [deadto(integer 0)]
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
  if `deadto' > 0 replace y = . if state <= `deadto'
end

*============================================ A: permute design universe
mkpanel , nstates(20) ntreat(8) deadto(4)
tempfile a0
save `a0', replace
dis "=== A1: FULL DATA (states 1-4 report no outcome) ==="
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(200) rseed(4321) dots(0) permsave(psaveA) replace
dis "A1_RC=" _rc
dis "A1_N=" e(N)
dis "A1_PERMREPS=" e(perm_reps)
cap dis "A1_PERMNPOSS=" e(perm_npossible)
mat pa = e(b)
dis "A1_B=" pa[1,1]
cap mat li e(perm_p)
ereturn list

dis "=== A2: DEAD ROWS DROPPED (oracle) ==="
use `a0', clear
drop if missing(y)
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(200) rseed(4321) dots(0) permsave(psaveB) replace
dis "A2_RC=" _rc
dis "A2_N=" e(N)
cap dis "A2_PERMNPOSS=" e(perm_npossible)
mat pb = e(b)
dis "A2_B=" pb[1,1]

dis "=== A3: which clusters did the design assign onsets to? ==="
use psaveA.dta, clear
qui ds onset_*
dis "A3_FULL_ASGCOLS=" `: word count `r(varlist)''
dis "A3_FULL_NAMES=" "`r(varlist)'"
use psaveB.dta, clear
qui ds onset_*
dis "A3_ORACLE_ASGCOLS=" `: word count `r(varlist)''
dis "A3_ORACLE_NAMES=" "`r(varlist)'"

*============================================ B: refusal from dead rows only
* state 7 is live and homogeneous IN SAMPLE, but carries one stale row with a
* different onset whose outcome is missing (a record with no usable data).
use `a0', clear
replace onset = 2009 if state == 7 & unit == 13 & year == 1999
replace y = . if state == 7 & unit == 13 & year == 1999
tempfile b0
save `b0', replace
dis "=== B1: FULL DATA, permute() ==="
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(50) rseed(1) dots(0)
dis "B1_RC=" _rc
dis "=== B2: the same run with that one dead row dropped ==="
use `b0', clear
drop if missing(y)
cap noi wooldid y unit year onset, cluster(state) permute(shuffle) permreps(50) rseed(1) dots(0)
dis "B2_RC=" _rc
dis "=== B3: and the ordinary (non-permute) run on the full data ==="
use `b0', clear
cap noi wooldid y unit year onset, cluster(state)
dis "B3_RC=" _rc
dis "B3_N=" e(N)

*============================================ C: permonsetcohorts bijection on dead rows
* cohort-level i (i == the onset grouping). One dead row gives a cohort a second onset.
clear
set seed 77
set obs 6
gen int coh = _n
gen int onset = 2003 + _n - 1
replace onset = . if coh > 4
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
replace onset = 2008 if clus == 1 & year == 1999
replace y = . if clus == 1 & year == 1999
dis "=== C2: same, plus ONE dead row carrying a stray onset ==="
cap noi wooldid y coh year onset, cluster(clus) permute(shuffle) permreps(50) permonsetcohorts rseed(2) dots(0)
dis "C2_RC=" _rc
dis "=== C3: and with that dead row dropped ==="
drop if missing(y)
cap noi wooldid y coh year onset, cluster(clus) permute(shuffle) permreps(50) permonsetcohorts rseed(2) dots(0)
dis "C3_RC=" _rc

*============================================ D: parse-time asserts on the raw `if'
* Realistic admin coding: unusable records carry the sentinel id -9 and no outcome.
mkpanel , nstates(12) ntreat(6)
tempfile d0
save `d0', replace
expand 2 if _n <= 20
sort state unit year
gen byte junk = 0
replace junk = 1 in 1/20
replace unit = -9 if junk == 1
replace y = .    if junk == 1
tempfile d1
save `d1', replace
dis "=== D1: 20 unusable rows coded unit == -9, y == . ==="
cap noi wooldid y unit year onset, cluster(state)
dis "D1_RC=" _rc
dis "=== D2: identical analysis sample, those rows dropped ==="
use `d1', clear
drop if junk == 1
cap noi wooldid y unit year onset, cluster(state)
dis "D2_RC=" _rc
dis "D2_N=" e(N)
dis "=== D3: full data but the user supplies the if wooldid would have applied ==="
use `d1', clear
cap noi wooldid y unit year onset if junk == 0, cluster(state)
dis "D3_RC=" _rc
dis "D3_N=" e(N)

dis "P3_ALLDONE"
