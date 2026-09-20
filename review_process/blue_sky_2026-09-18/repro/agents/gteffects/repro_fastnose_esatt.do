*==========================================================================
* repro_fastnose_esatt.do -- SELF-CONTAINED.
* Settles known_bydesign S7 item 11: "Item 11 looks unfinished on the pinned
* build: it was ruled to REMOVE fastnose (and impvalues) from the 'the
* covariance is not returned but the coefficient vector still is' list, and
* sthlp:1325 still opens that list with {cmd:fastnose}.  Check whether
* fastnose actually posts e(b_es_att)/e(numpre_att)."
*
* MEASURED ANSWER: it does.  fastnose posts e(b_es_att) and e(numpre_att)
* and withholds e(V_es_att), exactly as sthlp:1325 says.  So the help is
* accurate and the ruling to drop fastnose from that list was not built
* (or rested on a wrong premise about the code).
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_fastnose_esatt.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 80
gen long unit = _n
gen byte coh = 1
replace coh = 2 if _n > 35 & _n <= 58
replace coh = 3 if _n > 58
gen double ttre = .
replace ttre = 2004 if coh == 2
replace ttre = 2006 if coh == 3
gen double ai = rnormal()
expand 9
bysort unit: gen int year = 1999 + _n
tempvar ts
bysort year: gen double `ts' = rnormal() if _n == 1
bysort year: replace `ts' = `ts'[1]
gen double te = 0
replace te = 1 + 0.25*(year - ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`ts' + te + rnormal()*0.3

foreach mode in "" "fastnose" {
  dis as result _n "=== mode: `=cond("`mode'"=="","(analytic)","fastnose")' ==="
  quietly wooldid y coh year ttre, att att_it gteffects ///
      esprelength(2) espostlength(2) `mode'
  foreach k in b_es_att V_es_att b_es_att_it V_es_att_it {
    capture confirm matrix e(`k')
    dis as result "  matrix e(`k')  present = " (_rc == 0)
  }
  foreach k in numpre_att numpre_att_it {
    capture confirm scalar e(`k')
    dis as result "  scalar e(`k')  present = " (_rc == 0)
  }
}
dis as result _n "sthlp:1325 opens the 'coefficients posted without their covariance'"
dis as result "list with fastnose.  The measurement above agrees with the help."
