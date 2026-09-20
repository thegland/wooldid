*==============================================================================
* SELF-CONTAINED REPRO -- SAFETY-01 (my own confirmation, own data)
*
* CLAIM: wooldid's si8 maxvar preflight is NOT a lower bound on the default
* (reghdfe) engine's variable demand.  It refuses, with r(907), runs that
* complete and return the BIT-IDENTICAL answer the same default-safety call
* returns when maxvar is raised.  This contradicts ado:1341-1348 and the
* project record of 2026-09-16 ("zero false positives").
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_SAFETY01_si8_falsepositive.do
* Reads its log next to itself.  Needs nothing but the pinned build.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
dis "MAXVAR_AT_START=" c(maxvar)      // Stata's default is 5000

*------------------------------------------------------------------
* Generator: N units x 8 years, the FIRST 20 units ever-treated in
* three staggered cohorts.  The treated structure is IDENTICAL for
* every N; only the number of NEVER-TREATED CONTROL units changes.
*------------------------------------------------------------------
cap prog drop mkpanel
prog def mkpanel
  args nunits
  clear
  set obs `nunits'
  gen long unit = _n
  gen byte tg = (unit <= 20)
  expand 8
  bysort unit: gen int t = 1999 + _n
  gen double ttre = .
  replace ttre = 2002 + mod(unit,3) if tg == 1
  gen double ui = rnormal()
  bysort unit (t): replace ui = ui[1]
  gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
end

*==================================================================
* PART 1 -- the refusal threshold is driven by the NEVER-TREATED
* CONTROL COUNT, which contributes no regressor on the default engine
*==================================================================
mkpanel 4000
dis "=== P1a: 4000 units (3980 never-treated controls), default safety ==="
cap noi wooldid y unit t ttre, att cluster(unit)
dis "P1a_RC=" _rc

mkpanel 4900
dis "=== P1b: 4900 units -- SAME 20 treated units, 920 more pure controls ==="
cap noi wooldid y unit t ttre, att cluster(unit)
dis "P1b_RC=" _rc
dis "   ^ adding never-treated CONTROLS -- which add precision and add no"
dis "     regressor on the reghdfe path -- is what triggers the refusal."

*==================================================================
* PART 2 -- the refused run is not merely 'runnable'; it is BIT-IDENTICAL
* to the default-safety answer with maxvar headroom
*==================================================================
mkpanel 4800
save "_repro_s01.dta", replace
dis "=== P2a: 4800 units, DEFAULT safety at maxvar " c(maxvar) " ==="
cap noi wooldid y unit t ttre, att cluster(unit)
dis "P2a_RC=" _rc

dis "=== P2b: identical call + safety(off), SAME maxvar ==="
use "_repro_s01.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) safety(off)
dis "P2b_RC=" _rc
if _rc == 0 {
  mat bB = e(b)
  mat VB = e(V)
  scalar hB = bB[1,1]
  scalar hS = sqrt(VB[1,1])
  dis "P2b_N=" e(N)
  dis "P2b_ATT_HEX=" %21x hB
  dis "P2b_SE_HEX="  %21x hS
}

*==================================================================
* PART 3 -- the demand IS real on the reg/poisson path.  The fix has
* to be engine-conditional; simply dropping the level counts would
* reintroduce a late r(907) for customvce()/unconditionalse users.
*==================================================================
dis "MAXVAR_FOR_PART3=" c(maxvar)
mkpanel 5200
dis "=== P3a: default safety (expect the si8 refusal) ==="
cap noi wooldid y unit t ttre, att cluster(unit)
dis "P3a_RC=" _rc
dis "=== P3b: safety(off), DEFAULT engine (reghdfe) -- FALSE POSITIVE ==="
cap noi wooldid y unit t ttre, att cluster(unit) safety(off)
dis "P3b_RC=" _rc
if _rc == 0 dis "P3b_N=" e(N)
dis "=== P3c: safety(off), customvce -> reg path -- the demand is REAL ==="
cap noi wooldid y unit t ttre, att customvce(cluster unit) safety(off)
dis "P3c_RC=" _rc

* stash the safety(off) answer on disk, then reset the session so maxvar can be raised
clear
set obs 1
gen double attB = hB
gen double seB  = hS
save "_s01_b.dta", replace
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 20000
dis "=== P2c: identical call, DEFAULT safety, maxvar " c(maxvar) " ==="
use "_repro_s01.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit)
dis "P2c_RC=" _rc
if _rc == 0 {
  mat bC = e(b)
  mat VC = e(V)
  dis "P2c_N=" e(N)
  dis "P2c_ATT_HEX=" %21x bC[1,1]
  dis "P2c_SE_HEX="  %21x sqrt(VC[1,1])
  scalar cAtt = bC[1,1]
  scalar cSe  = sqrt(VC[1,1])
  preserve
  use "_s01_b.dta", clear
  dis "P2b_ATT_HEX_RECALL=" %21x attB[1]
  dis "P2b_SE_HEX_RECALL="  %21x seB[1]
  dis "ATT_BITIDENTICAL=" cond(attB[1] == cAtt, "YES", "NO")
  dis "SE_BITIDENTICAL="  cond(seB[1]  == cSe,  "YES", "NO")
  restore
}

dis "REPRO_SAFETY01_DONE"
