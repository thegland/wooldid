*! SELF-CONTAINED REPRO -- state-leak agent, finding 1
*!
*! When the outcome is missing on every never-treated row -- an ordinary
*! administrative-data shape -- wooldid exits r(459) with a message that names the
*! DESIGN ("you have no never-treated control group ... modify the sample window
*! manually") and never says that rows were removed for missingness.  Both remedies
*! the message offers are wrong for this data, and the one it names first,
*! safety(off), proceeds on a silently different sample (half the rows, 9 of 12
*! clusters, identified off not-yet-treated comparisons only) and reports it with
*! no caution: the point estimate moves 28% here, 48% on a DGP with dynamics.
*!
*! Run:  cd <this directory>
*!       "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_missing_control_outcome.do
*!       then read repro_missing_control_outcome.log
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

* ------------------------------------------------------------------
* An ordinary staggered panel: 60 counties, 12 years, three treated cohorts
* (2005 / 2007 / 2009) and 15 NEVER-TREATED counties.
clear
set obs 60
gen int unit  = _n
gen int state = mod(_n-1, 12) + 1
gen int ttre = .
replace ttre = 2005 if mod(_n,4) == 1
replace ttre = 2007 if mod(_n,4) == 2
replace ttre = 2009 if mod(_n,4) == 3
gen double ufe = rnormal()
expand 12
bys unit: gen int year = 2000 + _n
gen byte post = (year >= ttre) & !missing(ttre)
gen double y = 1 + ufe + 0.03*(year-2000) + 0.5*post + rnormal()

qui count if missing(ttre)
dis as text "never-treated ROWS in memory: " r(N)

* ------------------------------------------------------------------
* LEG 1 (control): with a complete outcome the run succeeds.
capture noisily wooldid y unit year ttre, cluster(state) att
dis "LEG1_CONTROL_RC=" _rc
dis "LEG1_N=" e(N)
capture confirm matrix e(wdidmainresults)
if _rc == 0 {
  tempname M1
  matrix `M1' = e(wdidmainresults)
  dis "LEG1_ATT=" %18.12g `M1'[1,1] "  SE=" %18.12g `M1'[1,2]
}

* ------------------------------------------------------------------
* LEG 2: the outcome is observed only for participants.  The DESIGN has not
* changed at all -- the never-treated counties are still there, still never
* treated, still inside the sample window.  Only the outcome column is missing
* on their rows.
qui replace y = . if missing(ttre)
qui count if missing(y)
dis as text "rows whose outcome is now missing: " r(N)

capture noisily wooldid y unit year ttre, cluster(state) att
dis "LEG2_RC=" _rc

* ------------------------------------------------------------------
* LEG 3: the first remedy the message offers.  It returns a number, with no
* caution, on a sample that contains NO untreated unit at all.
capture noisily wooldid y unit year ttre, cluster(state) att safety(off)
dis "LEG3_SAFETYOFF_RC=" _rc
dis "LEG3_N=" e(N) "  N_dropped=" e(N_dropped)
capture confirm matrix e(wdidmainresults)
if _rc == 0 {
  tempname M3
  matrix `M3' = e(wdidmainresults)
  dis "LEG3_ATT=" %18.12g `M3'[1,1] "  SE=" %18.12g `M3'[1,2] "  P=" %18.12g `M3'[1,4]
}

* ------------------------------------------------------------------
* LEG 4: the honest thing a user would do if told the truth -- drop the rows with
* no outcome -- reproduces the SAME r459, which is the correct answer for a design
* that genuinely has no controls.
qui drop if missing(y)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "LEG4_RC=" _rc

* ==================================================================
* MEASURED (2026-09-17, pinned build c1d8979):
*   LEG1  rc 0    e(N)=720   att = .251427121182  se = .100824299164
*   LEG2  rc 459
*         Caution: estimability checks found time periods exist in sample without
*                  a control unit, causing them to be dropped.
*         Error:   you have no never-treated control group and no not-yet-treated
*                  control group whose treatment falls outside the sample window.
*         To ignore this error, specify safety(off) or modify the sample window
*         manually.
*   LEG3  rc 0    e(N)=360  e(N_dropped)=360  N clusters = 9
*         att = .182200591071  se = .0610866413393  p = .0175291582039
*         -- no caution of any kind about the missing control group.
*         The point estimate moves 28% (this DGP has a constant treatment effect).
*         On the richer DGP used by p09_misc.do / p15b_safetyoff.do, which adds a
*         dynamic term, the same three legs give att = .795871565626 intact,
*         att = .415533924404 under safety(off) -- a 48% move.
*   LEG4  rc 459
*
* WHY THE EXPECTATION HOLDS: the markout at ado:683
*   qui replace `ifvar' = 0 if `i' >= . | `t' >= . | `y' >= . | `w' >= . | `w' <= 0
* runs long BEFORE the control-group screen at ado:952
*   qui count if `ttre' == . & `ifstatement'  ->  r(N)==0  ->  error 459
* so by the time the screen looks, wooldid has already removed every control row
* itself.  wooldid knows the count -- it is the same quantity the help says
* surfaces "in the header's dropped-observation count" (sthlp:898, sthlp:927,
* sthlp:1411) -- but that count is only printed at ~ado:1497, after the regression,
* so on this path it is never shown.  The refusal is right; the diagnosis is not,
* and neither remedy it offers addresses the actual cause.
dis "REPRO_DONE"
