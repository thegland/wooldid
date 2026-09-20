*==============================================================================
* repro_AR3.do -- SELF-CONTAINED.
* Claim: the replay-staleness guard is SELF-POISONING.  Typing a bare `wooldid'
* to redisplay results, at a moment when some other command owns e(), produces
* (a) an error about a missing VARLIST rather than "no wooldid results to
*     replay", and
* (b) records itself as a "failed wooldid call", after which even `estimates
*     restore' of a perfectly valid stored wooldid model is refused with
*     "the last wooldid call failed with r(198); the stored results predate that
*      failed call. Refit wooldid before replaying."
* Nothing is stale: the failed call WAS the user's own replay attempt.  The only
* documented escape is a full refit.
*
* Run: stata-mp -e do repro_AR3.do  and read the three marked blocks.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

quietly {
    set obs 120
    gen long i = _n
    gen byte cohort = mod(_n,4)
    gen double ttre = .
    replace ttre = 2004 if cohort==1
    replace ttre = 2006 if cohort==2
    replace ttre = 2008 if cohort==3
    gen double a_i = rnormal()
    gen double x1 = rnormal()
    expand 12
    bysort i: gen int t = 1999 + _n
    gen double y = a_i + 0.1*(t-2000) + 0.8*(ttre<. & t>=ttre) + 0.3*rnormal()
}

dis as result "==== 1. a perfectly good wooldid run ===="
wooldid y i t ttre, att robust espre(2) espost(2)
estimates store GOOD
dis "LASTRC after the good run = [$WOOLDID_LASTRC]   (empty = replay allowed)"

dis as result "==== 2. the user runs something else, then types `wooldid' to look again ===="
qui regress y x1
dis "e(cmd) is now [" e(cmd) "]"
cap noi wooldid
dis "rc = " _rc
dis "LASTRC is now = [$WOOLDID_LASTRC]   <-- the replay ATTEMPT poisoned the marker"

dis as result "==== 3. the user recovers the model explicitly -- still refused ===="
estimates restore GOOD
dis "e(cmd) after restore = [" e(cmd) "]"
cap noi wooldid
dis "rc = " _rc
dis "  ^ the message blames a failed estimation; no estimation failed."
dis "REPRO_AR3_DONE"
