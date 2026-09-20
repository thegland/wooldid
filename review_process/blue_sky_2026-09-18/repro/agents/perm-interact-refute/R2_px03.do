* R2_px03.do -- SELF-CONTAINED. Skeptic's independent test of PX-03.
* Claim under test: semielasticity on a NON-POSITIVE outcome fails, and under fastnose
* the run returns rc 0 with a printed table and NO e(b).
* I also separate two candidate triggers:
*   (a) the raw outcome takes non-positive values      <- their stated mechanism
*   (b) the model's PREDICTION crosses/approaches zero <- margins' eydx actually divides by this
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260917
set sortseed 20260917

*------------------------------------------------------------------ my own DGP
set obs 60
gen int firm = _n
gen double fe_i = rnormal(0, 1.4)
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
gen double size = 3 + rnormal(0, .8)
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
* ycross : a z-score-like outcome that straddles zero (realistic: standardized test score)
gen double ycross = fe_i + fe_t + effect + eps
* ypos   : the SAME outcome shifted strictly positive
gen double ypos = ycross + 12
* ymix   : strictly positive PREDICTIONS (mean ~ 12) but a handful of NEGATIVE raw values,
*          produced by a fat-tailed idiosyncratic shock. Tests trigger (a) vs (b).
gen double shock = 0
replace shock = -30 if mod(_n, 197) == 0
gen double ymix = ypos + shock
drop rel effect

summ ycross ypos ymix
count if ycross <= 0
count if ymix <= 0
count if ypos <= 0

*--------------------------------------------------- reporting helper
cap prog drop rpt
prog def rpt
  args tag
  dis "`tag'_ecmd=" "`e(cmd)'"
  cap confirm matrix e(b)
  dis "`tag'_hasb=" (_rc==0)
  if _rc == 0 dis "`tag'_bnames=" "`: colnames e(b)'"
  cap confirm matrix e(wdidmainresults)
  dis "`tag'_hasmain=" (_rc==0)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    matrix list WR
  }
end

*============ 1. baseline, no semielasticity, sign-crossing outcome
cap noi wooldid ycross firm year onset, att cluster(firm)
dis "T1_RC=" _rc
rpt T1

*============ 2. semielasticity + fastnose on the sign-crossing outcome  <- the accused
cap noi wooldid ycross firm year onset, att cluster(firm) semielasticity fastnose
dis "T2_RC=" _rc
rpt T2

*============ 3. semielasticity alone on the sign-crossing outcome
cap noi wooldid ycross firm year onset, att cluster(firm) semielasticity
dis "T3_RC=" _rc
rpt T3

*============ 4. semielasticity under permute() (fast-mode default)
cap noi wooldid ycross firm year onset, att cluster(firm) semielasticity permute(shuffle) permreps(30) rseed(99)
dis "T4_RC=" _rc
rpt T4

*============ 5. same outcome shifted strictly positive: all three paths
cap noi wooldid ypos firm year onset, att cluster(firm) semielasticity fastnose
dis "T5_RC=" _rc
rpt T5
cap noi wooldid ypos firm year onset, att cluster(firm) semielasticity
dis "T6_RC=" _rc
rpt T6
cap noi wooldid ypos firm year onset, att cluster(firm) semielasticity permute(shuffle) permreps(30) rseed(99)
dis "T7_RC=" _rc
rpt T7

*============ 6. TRIGGER SEPARATION: negative raw values, positive predictions
cap noi wooldid ymix firm year onset, att cluster(firm) semielasticity fastnose
dis "T8_RC=" _rc
rpt T8
cap noi wooldid ymix firm year onset, att cluster(firm) semielasticity
dis "T9_RC=" _rc
rpt T9

*============ 7. does plain fastnose (no semielasticity) post e(b) on ycross?
cap noi wooldid ycross firm year onset, att cluster(firm) fastnose
dis "T10_RC=" _rc
rpt T10
dis "---- R2 done ----"
