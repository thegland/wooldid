* ============================================================================
* SKEPTIC REPRO for hv-noparent-F1
* CLAIM UNDER TEST: nine options are accepted, do nothing, say nothing without
*   their parent, contradicting sthlp:1181 "No accepted option is a silent no-op".
* This probe tests (a) the BEHAVIOUR half by paired runs, (b) whether any output
*   is printed, (c) whether saveplots() writes files, and (d) the separate and
*   stronger cfxplottypes() VALIDATION-GAP claim.
* Self-contained. Generates its own binary staggered panel.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f1_noparent.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4471
set sortseed 4471

which wooldid

clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen sgv = mod(_n,2)
gen double ufe = rnormal()
gen double dose0 = runiform()*3
expand 10
bysort unit: gen int t = _n
gen double y = 2 + ufe + 0.1*t + rnormal()
replace y = y + 0.8*(t>=ttre) if ttre<.
gen double dose = dose0*(t>=ttre & ttre<.)
gen double w = 1 + runiform()
tempfile dat
save `dat', replace

* clear out any old artifacts so "file written?" is a real test
capture erase p1_ES.dta
capture erase p1_ES_att.gph
capture erase p1_cfx.dta
capture erase p2_ES.dta
capture erase p2_ES_att.gph

* ---------------- BASELINE ---------------------------------------------------
use `dat', clear
capture noisily wooldid y unit t ttre, att robust
dis "BASE_RC=" _rc
matrix MB = e(wdidmainresults)
matrix VB = e(V)
matrix bB = e(b)
dis "BASE_ATT=" %21x el(MB,1,1)

* ---------------- one leg per orphan option ----------------------------------
* leg tag, then the option text
local L1 "saveplots(p1)"
local L2 "makeplots saveplots(p2)"
local L3 "ccnoabsorb"
local L4 "contreatpoly(3)"
local L5 "suppresssgprimaryeffects"
local L6 "cfxplottypes(post 1)"
local L7 "cfxplotnose"
local L8 "makecfxplotsbysg"
local L9 "latticeignoreweights"

forval k = 1/9 {
  use `dat', clear
  dis "@@@@ LEG`k' OPT=[`L`k''] @@@@"
  capture noisily wooldid y unit t ttre, att robust `L`k''
  dis "LEG`k'_RC=" _rc
  matrix M`k' = e(wdidmainresults)
  dis "LEG`k'_MRELDIF_MAIN=" mreldif(M`k', MB)
  capture matrix V`k' = e(V)
  capture dis "LEG`k'_MRELDIF_V=" mreldif(V`k', VB)
  capture matrix b`k' = e(b)
  capture dis "LEG`k'_MRELDIF_B=" mreldif(b`k', bB)
  dis "@@@@ ENDLEG`k' @@@@"
}

* ---------------- did saveplots() write anything? ----------------------------
foreach f in p1_ES.dta p1_ES_att.gph p1_cfx.dta p2_ES.dta p2_ES_att.gph {
  capture confirm file "`f'"
  dis "FILE_`f'_EXISTS=" (_rc==0)
}

* ---------------- CONTROL: saveplots WITH a horizon request ------------------
use `dat', clear
capture noisily wooldid y unit t ttre, att robust makeplots esprelength(2) espostlength(2) saveplots(p3) replace
dis "CTRL_RC=" _rc
foreach f in p3_ES.dta p3_ES_att.gph {
  capture confirm file "`f'"
  dis "FILE_`f'_EXISTS=" (_rc==0)
}

* ============================================================================
* THE VALIDATION-GAP LEG: sthlp:974 says, with NO parent-option qualifier,
*   "a token that is neither post nor pre must be a number, and one that is
*    neither (postpost, 2post2, a variable name) stops the run with an error
*    naming the entry, as does a numeric token that is missing or not an integer."
* ============================================================================
use `dat', clear
capture noisily wooldid y unit t ttre, att robust cfxplottypes(psot)
dis "VG1_noparent_badtoken_RC=" _rc
use `dat', clear
capture noisily wooldid y unit t ttre, att robust cfxplottypes(1.5)
dis "VG2_noparent_nonint_RC=" _rc
use `dat', clear
capture noisily wooldid y unit t ttre, att robust cfxplottypes(nosuchvariable)
dis "VG3_noparent_varname_RC=" _rc
use `dat', clear
capture noisily wooldid y unit t ttre, att robust cfxplottypes(postpost 2post2)
dis "VG4_noparent_multi_RC=" _rc

* with contreat() but WITHOUT makecfxplots -- still no parent for the plots
use `dat', clear
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(basic) cfxplottypes(psot)
dis "VG5_contreat_nomakecfx_RC=" _rc

* the parent present: the documented error should fire
use `dat', clear
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(basic) makecfxplots cfxplottypes(psot)
dis "VG6_fullparent_RC=" _rc

* and the contreatpoly() contrast, which sthlp:917 says IS validated regardless
use `dat', clear
capture noisily wooldid y unit t ttre, att robust contreatpoly(0)
dis "VG7_contreatpoly0_noparent_RC=" _rc

dis "DONE_F1"
