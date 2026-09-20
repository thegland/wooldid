*==============================================================================
* SELF-CONTAINED REPRO -- FP-NOPRE-LATEENTRY (skeptic's own version)
*
* ado:12053  nopre = max(ttre <= mint) by i, where mint = min(t) by i -- i.e. the
*            screen is PER UNIT against that unit's own first observed period.
*            A unit that enters the panel late and whose onset sits comfortably
*            inside the sample window is dropped as well.
* sthlp:73   describes the screen only as "cohorts treated before the sample
*            begins are dropped".  A help-wide grep for "first observed",
*            "no pre-treatment", "without a pre" and "first usable" finds only
*            sthlp:824, which is the anticipation() twin of the same screen --
*            and that twin DOES print a caution (ado:12036-12046).
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_nopre_own.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606
set sortseed 606
which wooldid

* ---- ARM 1: units 11-16 enter the panel at t=4, which is exactly their onset
clear
set obs 30
gen int i = _n
gen int ttre = .
replace ttre = 5 if i <= 10
replace ttre = 4 if i > 10 & i <= 16
gen double alpha = 50 + i
expand 8
bysort i: gen int t = _n
drop if i > 10 & i <= 16 & t < 4
gen byte post = (t >= ttre) & ttre < .
gen double eff = cond(i > 10 & i <= 16, 6, 2)
gen double y = alpha + 0.4*t + eff*post + rnormal()*0.5
qui count
dis "ARM1 rows = " r(N) "  (sample window is t = 1..8; every onset is inside it)"
tab ttre, m
dis _n "########## ARM 1: late entrants, onset == first observed period ##########"
cap noi wooldid y i t ttre, att
dis "RC_ARM1=" _rc
if _rc==0 dis "ARM1 att = " %13.7f e(b)[1,1] "   e(N) = " e(N) "   e(N_dropped) = " e(N_dropped)

* ---- ARM 2: the same units observed from t=3, so they have ONE pre period
clear
set obs 30
gen int i = _n
gen int ttre = .
replace ttre = 5 if i <= 10
replace ttre = 4 if i > 10 & i <= 16
gen double alpha = 50 + i
expand 8
bysort i: gen int t = _n
drop if i > 10 & i <= 16 & t < 3
gen byte post = (t >= ttre) & ttre < .
gen double eff = cond(i > 10 & i <= 16, 6, 2)
gen double y = alpha + 0.4*t + eff*post + rnormal()*0.5
dis _n "########## ARM 2: same units, one pre period each ##########"
cap noi wooldid y i t ttre, att
dis "RC_ARM2=" _rc
if _rc==0 dis "ARM2 att = " %13.7f e(b)[1,1] "   e(N) = " e(N) "   e(N_dropped) = " e(N_dropped)

dis _n "EXPECTED: ARM1's header reports the dropped count but names no cause, and no"
dis "caution fires -- while the anticipation() twin of the same screen prints one."
dis "The estimand also changes: the six dropped cohorts carry a true effect of 6"
dis "against the others' 2."
dis "REPRO_NOPRE_DONE"
