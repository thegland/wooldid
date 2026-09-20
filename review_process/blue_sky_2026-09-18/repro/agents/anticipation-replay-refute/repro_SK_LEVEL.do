*==============================================================================
* SELF-CONTAINED repro: level() / WOOLDID_LEVEL0 restore
* Run:  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_SK_LEVEL.do
* Reads repro_SK_LEVEL.log afterwards; every result line is tagged "##".
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*---- data: ordinary staggered panel, 60 units x 12 periods, 2 cohorts + never
set obs 60
gen long uid = _n
gen byte grp = mod(_n,3)
gen ttre = .
replace ttre = 5 if grp==1
replace ttre = 8 if grp==2
gen double ufe = rnormal()
expand 12
bysort uid: gen int t = _n
gen double y = ufe + 0.10*t + rnormal()
replace y = y + 0.8*(t - ttre + 1) if ttre < . & t >= ttre
tempfile D
save `D'

*==============================================================================
* 0. Stata itself accepts every one of these session levels.
*==============================================================================
foreach L in 99.9 99.99 99.91 32.22 {
  cap set level `L'
  dis "## SETLEVEL `L' rc=" _rc " c(level)=" c(level)
}
set level 95

*==============================================================================
* 1. level(99.9) and level(99.99) (the help's own upper bound) are REFUSED.
*    level(99.91) is accepted.  (option path)
*==============================================================================
foreach L in 99.9 99.99 99.91 {
  use `D', clear
  set level 95
  macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
  cap noi wooldid y uid t ttre, att robust level(`L')
  dis "## OPT level(`L') rc=" _rc " e(N)=" e(N) " e(level)=" e(level)
}

*==============================================================================
* 2. With NO level() option at all, a legal non-integer SESSION level kills it.
*    32.22 and 99.91 are in this set although 99.91 passes as an option in (1).
*==============================================================================
foreach L in 99.9 99.99 99.91 32.22 95 {
  use `D', clear
  macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
  set level `L'
  cap noi wooldid y uid t ttre, att robust
  dis "## SESS `L' (no level() typed) rc=" _rc " e(N)=" e(N) " LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC]"
  set level 95
}

*==============================================================================
* 3. A run whose estimation SUCCEEDS returns rc 198, loses the caller's level,
*    and leaves the shell globals set.
*==============================================================================
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 99.9
dis "## C0 c(level) before = " c(level)
cap noi wooldid y uid t ttre, att robust level(95)
mat B = e(b)
dis "## C1 rc=" _rc " e(N)=" e(N) " e(cmd)=[`e(cmd)'] e(level)=" e(level) " b_att=" B[1,1]
dis "## C2 c(level) AFTER = " c(level) "   (caller asked for 99.9)"
dis "## C3 LEVEL0=[$WOOLDID_LEVEL0] CLEANUP=[$WOOLDID_CLEANUP] LASTRC=[$WOOLDID_LASTRC]"

*==============================================================================
* 4. The replay-staleness guard is NOT armed by that failure, so a bare
*    `wooldid` afterwards REPLAYS the pre-failure results.
*==============================================================================
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
qui wooldid y uid t ttre, att robust
dis "## D0 good run e(N)=" e(N) " LASTRC=[$WOOLDID_LASTRC]"
set level 99.9
cap noi wooldid y uid t ttre, att_i robust
dis "## D1 FAILED call rc=" _rc " LASTRC=[$WOOLDID_LASTRC]  <- guard NOT armed"
cap noi wooldid
dis "## D2 bare replay rc=" _rc " e(N)=" e(N) "  <- pre-failure results redisplayed"
*---- control: an ordinary failure DOES arm the guard
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
qui wooldid y uid t ttre, att robust
cap noi wooldid y uid t ttre, att robust bogusoption
dis "## D3 ordinary failure rc=" _rc " LASTRC=[$WOOLDID_LASTRC]  <- guard armed"
cap noi wooldid
dis "## D4 bare replay rc=" _rc "  <- correctly refused"

*==============================================================================
* 5. NEW: while WOOLDID_LEVEL0 is stale, every later FAILING call has its true
*    return code replaced by 198 and the shell's rc-gated cleanup epilogue is
*    skipped (so LASTRC is never updated either).
*==============================================================================
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
cap noi wooldid y uid t ttre, att robust cluster(nosuchvar)
dis "## E0 baseline bad-cluster rc=" _rc " (expect 111) LASTRC=[$WOOLDID_LASTRC]"
set level 99.9
cap noi wooldid y uid t ttre, att robust level(95)
dis "## E1 poisoning call rc=" _rc " LEVEL0=[$WOOLDID_LEVEL0]"
cap noi wooldid y uid t ttre, att robust cluster(nosuchvar)
dis "## E2 SAME bad-cluster rc=" _rc " (was 111) LASTRC=[$WOOLDID_LASTRC] LEVEL0=[$WOOLDID_LEVEL0]"

set level 95
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
dis "## REPRO_SK_LEVEL_DONE"
