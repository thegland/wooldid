*==========================================================================
* repro_weightcols.do -- SELF-CONTAINED.
* CLAIM: the gteffects help says the weight columns are "missing for estimands
* not requested".  Requesting att_i (or att_itime) ALONE returns a fully
* populated w_att_it column although att_it was not requested, because ado:1643
* generates w_att_it whenever any of att_it / att_itime / att_i is asked for and
* the builder (ado:1758-1762) fills a column for every weight variable that
* exists.  Benign - the values are the right att_it weights - but the sentence
* is not accurate.
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_weightcols.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 60
gen long unit = _n
gen double ttre = .
replace ttre = 2004 if _n > 30 & _n <= 45
replace ttre = 2006 if _n > 45
gen double ai = rnormal()
expand 9
bysort unit: gen int year = 1999 + _n
tempvar ts
bysort year: gen double `ts' = rnormal() if _n == 1
bysort year: replace `ts' = `ts'[1]
gen double te = 0
replace te = 1 + 0.25*(year - ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`ts' + te + rnormal()*0.3
tempfile d
save `d'

foreach req in "att" "att_it" "att_i" "att_itime" {
  use `d', clear
  quietly wooldid y unit year ttre, `req' gteffects esprelength(2) espostlength(2)
  preserve
  clear
  quietly svmat double e(gteffects), names(col)
  local line ""
  foreach v in w_att w_att_it w_att_i w_att_itime {
    capture confirm variable `v'
    if _rc == 0 {
      quietly count if `v' < .
      local line "`line'  `v'=`r(N)'"
    }
    else local line "`line'  `v'=ABSENT"
  }
  dis as result "requested [`req'] :`line'"
  restore
}
dis as result _n "Expected from the help: every column but the requested one's is missing."
dis as result "Observed: [att_i] and [att_itime] both return a fully populated w_att_it."
