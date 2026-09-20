*==============================================================================
* repro_AR1.do  -- SELF-CONTAINED.
* Claim: a bare `wooldid' replay reprints the tables byte-for-byte but DROPS the
* cautions that disown the numbers in them.  sthlp:1458 says the dropped cautions
* "describe the estimation run, not the results"; the two shown here describe the
* results.
*
* Run:  stata-mp -e do repro_AR1.do    then read repro_AR1_live.txt vs
*       repro_AR1_replay.txt (or just the diff printed at the end of the log).
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local HERE "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

* ---- a plain staggered panel: 120 units, 12 years, three cohorts ------------
quietly {
    set obs 120
    gen long i = _n
    gen byte cohort = mod(_n,4)
    gen double ttre = .
    replace ttre = 2004 if cohort==1
    replace ttre = 2006 if cohort==2
    replace ttre = 2008 if cohort==3
    gen double a_i = rnormal()
    expand 12
    bysort i: gen int t = 1999 + _n
    gen double y = a_i + 0.1*(t-2000) + 0.8*(ttre<. & t>=ttre) + 0.3*rnormal()
}

* commands issued from inside a program so the do-file echo does not pollute
cap prog drop LIVE
prog def LIVE
  wooldid y i t ttre, att att_it cluster(i) espre(3) espost(3) jointtests adversarial
end
cap prog drop REPLAY
prog def REPLAY
  wooldid
end

log using "`HERE'/repro_AR1_live.txt", replace text name(L1)
LIVE
log close L1

log using "`HERE'/repro_AR1_replay.txt", replace text name(L2)
REPLAY
log close L2

dis as result "==== lines present LIVE but absent on REPLAY ===="
* print them by hand so the finding is visible without leaving Stata
type "`HERE'/repro_AR1_live.txt", lines(24)
dis as result "==== the replay starts here instead ===="
type "`HERE'/repro_AR1_replay.txt", lines(14)
dis "REPRO_AR1_DONE"
