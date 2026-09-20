*==========================================================================
* repro_lateentry_silentdrop.do -- SELF-CONTAINED (no includes, no helpers).
*
* CLAIM: wooldid silently removes every treated `i' whose OWN first observed
* period is at or after its treatment onset.  In an unbalanced panel that is
* not "cohorts treated before the sample begins are dropped" (sthlp:73) - the
* sample begins in 2000 and these firms are treated in 2004 - it is a
* unit-level rule keyed on `gegen min(t), by(i)' (ado:12018-12020, 12055-12059).
* Nothing is printed: no caution, no named count, no row in e(gteffects).
* The identical drop performed twelve lines earlier for anticipation()
* (ado:12045) DOES print "Caution: anticipation(#) places the effective onset
* of N unit(s) at or before their first usable pre-period; those unit(s) and
* their M observations are dropped from the sample."
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_lateentry_silentdrop.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* An administrative register, 2000-2010.
*   firms  1-40 : never treated, present throughout
*   firms 41-60 : treated in 2004, present throughout
*   firms 61-70 : treated in 2004, but they only ENTER the register in 2004
*                 (they exist, they are treated, they have 7 post-treatment years)
* Their treatment effect is 2.5x the others', so dropping them moves the number.
capture program drop mkreg
program define mkreg
  syntax , ENTry(int)
  clear
  set obs 70
  gen long unit = _n
  gen byte late = (_n > 60)
  gen double ttre = cond(_n > 40, 2004, .)
  gen double ai = rnormal()
  expand 11
  bysort unit: gen int year = 1999 + _n
  drop if late == 1 & year < `entry'
  tempvar ts
  bysort year: gen double `ts' = rnormal() if _n == 1
  bysort year: replace `ts' = `ts'[1]
  gen double te = 0
  replace te = 1 + 0.25*(year - ttre) if year >= ttre & ttre < .
  replace te = te*2.5 if late == 1 & year >= ttre & ttre < .
  gen double y = 2 + ai + 0.5*`ts' + te + rnormal()*0.3
end

foreach e in 2000 2003 2004 2005 {
  dis as result _n "############ late entrants first observed in `e' ############"
  mkreg , entry(`e')
  quietly count
  local nobs = r(N)
  tempvar tg
  bysort unit: gen byte `tg' = (_n == 1)
  quietly count if `tg' == 1 & ttre < .
  dis as result "INPUT: `nobs' observations, " r(N) " treated firms"
  drop `tg'
  noisily wooldid y unit year ttre, att att_it att_i att_itime gteffects ///
       esprelength(2) espostlength(2)
  dis as result "  e(N) = " e(N) "   e(N_dropped) = " e(N_dropped)
  preserve
  clear
  quietly svmat double e(gteffects), names(col)
  bysort i: keep if _n == 1
  quietly count
  dis as result "  DISTINCT TREATED FIRMS WITH A ROW IN e(gteffects): " r(N)
  restore
  dis as result "  att = " %9.6f e(b)[1, colnumb(e(b), "main_att")]
}

dis as result _n "Read the four blocks together: with entry in 2000 or 2003 all 30"
dis as result "treated firms appear; with entry in 2004 or 2005 only 20 do, e(N_dropped)"
dis as result "jumps, att moves, and the ONLY trace in the printed output is the"
dis as result "generic '(N obs dropped from initial sample)' line."
