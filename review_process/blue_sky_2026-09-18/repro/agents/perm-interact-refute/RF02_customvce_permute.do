* RF02_customvce_permute.do
* SELF-CONTAINED repro for PX-02 (skeptic's independent version).
* CLAIM UNDER TEST: permute()+customvce() is unreachable, so the customvce() forwarding
* line in the permute replication string (ado:128) is dead code.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 771103
set sortseed 771103

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
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel if onset < . & rel >= 0
gen double y = 2.0 + fe_i + .12*(year-1997) + effect + eps
drop rel effect

foreach spec in 1 2 3 4 5 {
  ereturn clear
  if `spec' == 1 cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(30) customvce(cluster firm)
  if `spec' == 2 cap noi wooldid y firm year onset, att permute(shuffle) permreps(30) customvce(cluster firm)
  if `spec' == 3 cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(30) customvce(robust) safety(off)
  if `spec' == 4 cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(30) customvce()
  if `spec' == 5 cap noi wooldid y firm year onset, att cluster(firm year) permute(shuffle) permreps(30) customvce(cluster firm)
  dis "SPEC`spec'_RC=" _rc
}
dis "---- RF02 done ----"
