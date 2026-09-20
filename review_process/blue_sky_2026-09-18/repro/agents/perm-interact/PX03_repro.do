*==============================================================================
* PX-03 self-contained repro.
*
* semielasticity on an outcome that is NOT STRICTLY POSITIVE fails inside
* margins, and wooldid never diagnoses it.  Four observed failure modes, all on
* ordinary staggered panels, and a constant shift of the same outcome to be
* strictly positive makes every one of them run clean:
*
*   (1) rc 0, a full results table with an all-missing att row and R2 lines,
*       and NO e(b) posted at all            [panel A, + fastnose]
*   (2) r(3598), raw Mata traceback naming an internal tempvar
*                                            [panel A, no fastnose]
*   (3) r(459) "inconsistent estimation sample levels 0 and 1 of factor
*       __000007" - Stata's own margins error naming a wooldid tempvar
*                                            [panel B]
*   (4) r(459) "the observed comparison statistic for att could not be read from
*       the rep-0 run; permutation inference is impossible" - which blames the
*       comparator run when the OBSERVED run produced no coefficients either
*                                            [panel A, + permute()]
*
* permute() matters here because permute() SELECTS fastnose by default, so (1)
* is the branch a permutation user is on and (4) is what they are told.
*
* The only external dependency is the blue-sky review's shared data generator,
* at the absolute path in the include below.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260916
set sortseed 20260916

cap prog drop shw
prog def shw
  syntax anything(name=tag) [, OPTS(string asis) YV(string)]
  cap noi wooldid `yv' unit t ttre , att cluster(unit) `opts'
  local rc = _rc
  cap confirm matrix e(b)
  local hb = (_rc==0)
  dis "PX03|`tag'|RC=" `rc' "|HAS_e(b)=" `hb' "|e(cmd)=`e(cmd)'|bnames=`: colfullnames e(b)'"
end

*--------------------------------------------------------------------------
* PANEL A - 40 units, 8 years, onsets 2003/2005, heterogeneous effects
*--------------------------------------------------------------------------
local gseed   20260916
local gN      40
local gT      8
local gt0     2000
local gonsets 2003 2005
local gysd    0.30
local gtsd    0.40
local gaisd   1
local gtaumode hetero
local gdose   cellconst
local gwtype  unit
local gpois   1
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
gen double ypos = y + 10                     // identical shape, strictly positive
qui sum y
dis "PX03|panelA y range = " r(min) " .. " r(max)

dis _n "### A0. sanity: the plain run on this outcome is fine"
shw A_plain , yv(y)

dis _n "### A1. THE WORST CASE: semielasticity + fastnose -> rc 0, table printed, NO e(b)"
shw A_semi_fastnose , yv(y) opts(semielasticity fastnose)
dis "### ^^ scroll up: every column of the att row above is '.', an R2 line was"
dis "### printed, rc is 0, and e(b) does not exist.  Nothing was said."

dis _n "### A2. what a permute() user sees (permute selects fastnose by default)"
shw A_semi_permute , yv(y) opts(semielasticity permute(shuffle) permreps(30) rseed(777) dots(0))

dis _n "### A3. without fastnose: raw Mata traceback"
shw A_semi_full , yv(y) opts(semielasticity)

dis _n "### A4. the SAME outcome shifted strictly positive: all three are clean"
shw A_pos_fastnose , yv(ypos) opts(semielasticity fastnose)
shw A_pos_permute  , yv(ypos) opts(semielasticity permute(shuffle) permreps(30) rseed(777) dots(0))
shw A_pos_full     , yv(ypos) opts(semielasticity)

*--------------------------------------------------------------------------
* PANEL B - a plainer inline panel: the third failure mode
*--------------------------------------------------------------------------
clear
set seed 20260916
set sortseed 20260916
set obs 40
gen long unit = _n
gen double ttre = .
replace ttre = 2003 if unit <= 12
replace ttre = 2005 if unit > 12 & unit <= 24
expand 8
bysort unit: gen int t = 1999 + _n
gen double a_i = rnormal()
gen double g_t = 0.4*rnormal()
bysort t (unit): replace g_t = g_t[1]
gen double tau = 2.0 + rnormal()*0.5
bysort unit (t): replace tau = tau[1]
gen double y = a_i + g_t + tau*(ttre<. & t>=ttre) + rnormal()*0.6
gen double ypos = y + 10
qui sum y
dis "PX03|panelB y range = " r(min) " .. " r(max)
dis _n "### B1. semielasticity + fastnose -> r(459) 'inconsistent estimation sample"
dis "###     levels 0 and 1 of factor __000007' - Stata's margins error, naming a"
dis "###     wooldid tempvar, with no wooldid sentence"
shw B_semi_fastnose , yv(y) opts(semielasticity fastnose)
dis _n "### B2. same outcome, shifted positive: clean"
shw B_pos_fastnose , yv(ypos) opts(semielasticity fastnose)
dis "PX03_DONE"
