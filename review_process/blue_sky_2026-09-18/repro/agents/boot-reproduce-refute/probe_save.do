* =================================================================
* probe_save.do -- a probe NOBODY suggested: the THREE records of one
* bootstrap run must agree.
*   (1) e(bootdist)            -- the replication coefficient draws in e()
*   (2) bootsave() <file>.dta  -- the same draws on disk
*   (3) the file's _rep == 0 row vs e(b) -- the observed comparator
* plus the documented dataset characteristics
* _dta[wooldid_engine] / [wooldid_seed] / [wooldid_rngstate] (sthlp, rseed()).
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_save.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 246810
set sortseed 246810

set obs 24
gen long unit = _n
gen int ttre = .
replace ttre = 2004 if unit <= 8
replace ttre = 2006 if unit > 8 & unit <= 16
gen double ufe = rnormal()
expand 8
bysort unit: gen int year = 2000 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.5, 0))
gen double y  = 1 + ufe + 0.1*(year-2000) + te + rnormal()*0.7
save "sk_sv.dta", replace

cap erase "sksv.dta"
cap noi wooldid y unit year ttre, att att_i cluster(unit) bootstrap(standard) reps(50) rseed(90210) bootsave(sksv) dots(0)
dis "SV_RC=" _rc
matrix D = e(bootdist)
local dn : colnames D
dis "SV_DISTNAMES=`dn'"
dis "SV_DISTROWS=" rowsof(D) "  COLS=" colsof(D)
matrix Bobs = e(b)
local bn : colnames Bobs
dis "SV_BNAMES=`bn'"
dis "SV_B_ATT="   %21x el(Bobs,1,1)
dis "SV_B_ATTI="  %21x el(Bobs,1,2)
scalar sv_state = 0
local esd = e(boot_rngstate)
local esed = e(boot_seed)
dis "SV_E_SEED=`esed'"

preserve
use "sksv.dta", clear
dis "SV_FILEOBS=" _N
qui ds
dis "SV_FILEVARS=" "`r(varlist)'"
qui sum _rep
dis "SV_REPMIN=" r(min) " REPMAX=" r(max)
* rep-0 row must carry the observed estimates
qui levelsof main_att if _rep == 0, local(f0att)
qui levelsof main_att_i if _rep == 0, local(f0atti)
dis "SV_FILE0_ATT="  %21x `f0att'
dis "SV_FILE0_ATTI=" %21x `f0atti'
* the draws themselves
gen long ord = _n
qui sum _rep
forvalues r = 1/5 {
  qui levelsof main_att if _rep == `r', local(v)
  dis "SV_FILEREP`r'_ATT=" %21x `v'
}
dis "SV_CHAR_ENGINE="  `"`: char _dta[wooldid_engine]'"'
dis "SV_CHAR_SEED="    `"`: char _dta[wooldid_seed]'"'
local cs : char _dta[wooldid_rngstate]
dis "SV_CHAR_STATE_PRESENT=" (`"`cs'"' != "")
restore

forvalues r = 1/5 {
  dis "SV_EDISTREP`r'_ATT=" %21x el(D,`r',1)
}

dis "PROBE_SAVE_DONE"
