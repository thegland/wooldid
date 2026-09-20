* repro_cwm_missing_error.do -- SELF-CONTAINED.
* (1) customweightmultiplier() missing values mark the row OUT of the estimation screen
*     (ado:704).  With a multiplier that administrative data supplies only for program
*     participants, every never-treated row is marked out and wooldid stops with an error
*     that states something FALSE about the user's data -- "you have no never-treated
*     control group" -- while the dataset plainly has nine never-treated counties.
*     Nothing names customweightmultiplier(); both suggested remedies are wrong.
* (2) the negative-cwm refusal at ado:718-723 tells the user "To ignore this error,
*     specify safety(off)".  Following that advice dies r(402) "negative weights
*     encountered" with no wooldid context at all.
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_cwm_missing_error.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 20
gen int county = _n
gen int onset = .
replace onset = 2011 if inrange(county,1,4)
replace onset = 2013 if inrange(county,5,8)
replace onset = 2015 if inrange(county,9,11)
gen double cfe = rnormal(0,1)
* program enrolment size: a POSITIVE, time-constant administrative attribute that exists
* only for counties that ever ran the program.  Missing for the nine never-treated ones.
gen double progsize = .
replace progsize = 200 + 40*county if onset < .
expand 10
bysort county: gen int year = 2007 + _n
expand 3
bysort county year: gen int estab = _n
gen double sw = exp(rnormal(1.0,0.5))
gen byte post = (onset<. & year>=onset)
gen double y = 3 + cfe + 0.12*(year-2007) + 2.4*post + rnormal(0,0.8)
qui count if onset >= .
dis "NEVER_TREATED_ROWS=" r(N)
qui levelsof county if onset >= ., local(nt)
dis "NEVER_TREATED_COUNTIES=" wordcount("`nt'")
qui count if progsize < 0
dis "NEGATIVE_MULTIPLIER_ROWS=" r(N)
qui sum progsize
dis "MULTIPLIER_MIN=" r(min) "  MAX=" r(max)

*--- (a) baseline: runs fine with no multiplier
cap noi wooldid y county year onset [aw=sw], att
dis "A_RC=" _rc
dis "A_eN=" e(N)

*--- (b) THE FINDING: the same run with the participants-only multiplier
cap noi wooldid y county year onset [aw=sw], att customweightmultiplier(progsize)
dis "B_RC=" _rc

*--- (c) the first remedy the message offers
cap noi wooldid y county year onset [aw=sw], att customweightmultiplier(progsize) safety(off)
dis "C_RC=" _rc

*--- (d) the ACTUAL remedy, which no message names: fill the multiplier in on controls
gen double progsize0 = cond(missing(progsize), 0, progsize)
cap noi wooldid y county year onset [aw=sw], att customweightmultiplier(progsize0)
dis "D_RC=" _rc
dis "D_eN=" e(N)

*--- (e) SECOND FINDING: safety(off) is advertised as the bypass for a negative
*        multiplier, and does not bypass it.
gen double signedmult = cond(missing(progsize), 1, progsize/100)
replace signedmult = -0.5 if county == 3 & year <= 2009
cap noi wooldid y county year onset [aw=sw], att customweightmultiplier(signedmult)
dis "E_RC=" _rc
cap noi wooldid y county year onset [aw=sw], att customweightmultiplier(signedmult) safety(off)
dis "F_RC=" _rc

dis "MISSINGCWM_DONE"
