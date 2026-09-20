* r00_setup.do -- assert the pinned build, build my own panel, time one bootstrap rep
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* ---- my own generator: 30 states x 4 counties, 18 periods, state-level staggered onset
clear
set obs 30
gen state = _n
gen onset = .
replace onset = 8  + mod(_n, 7) if mod(_n,3)==0        // ~10 treated states
expand 4
bysort state: gen cty = _n
gen county = state*10 + cty
expand 18
bysort county: gen t = _n
gen double ushock = 0
bysort county (t): replace ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double sshock = 0
bysort state t: gen byte first_st = (_n==1)
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
label var y "outcome"
compress
save "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-interact-refute/panelA.dta", replace
tab onset, m
dis "NSTATES_TREATED=" 
qui levelsof state if onset<., local(ts)
dis "TREATEDSTATES=" `: word count `ts''

timer clear
timer on 1
cap noi wooldid y county t onset, cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) jointtests pretrendpower(80)
dis "ANALYTIC_RC=" _rc
timer off 1
timer list 1

timer on 2
cap noi wooldid y county t onset, cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) jointtests pretrendpower(80) bootstrap(standard) reps(10) rseed(9001)
dis "BOOT10_RC=" _rc
timer off 2
timer list 2
