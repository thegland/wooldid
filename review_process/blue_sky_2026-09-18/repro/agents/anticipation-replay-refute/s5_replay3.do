clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

use "dat.dta", clear

*================= AR-3: bare wooldid while another cmd owns e() =============
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
qui wooldid y uid t ttre, att robust
est store GOOD
dis "@@G0 good ecmd=[`e(cmd)'] eN=" e(N) " LASTRC=[$WOOLDID_LASTRC]"

qui regress y x1
dis "@@G1 after regress ecmd=[`e(cmd)'] LASTRC=[$WOOLDID_LASTRC]"

cap noi wooldid
dis "@@G2 bare-wooldid-foreign rc=" _rc " LASTRC=[$WOOLDID_LASTRC] ecmd=[`e(cmd)']"

est restore GOOD
dis "@@G3 restored ecmd=[`e(cmd)'] eN=" e(N) " LASTRC=[$WOOLDID_LASTRC]"

cap noi wooldid
dis "@@G4 bare-after-restore rc=" _rc

*----- control: is the same true with NO e() at all? -----
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
cap noi wooldid
dis "@@G5 bare-no-e rc=" _rc " LASTRC=[$WOOLDID_LASTRC]"

*----- control: does a foreign e() alone (never a wooldid run) behave same? ---
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui wooldid y uid t ttre, att robust
est store GOOD2
qui regress y x1
cap noi wooldid
dis "@@G6 bare-foreign rc=" _rc " LASTRC=[$WOOLDID_LASTRC]"
* recovery: does ANY successful wooldid clear it?
qui wooldid y uid t ttre, att robust
dis "@@G7 after-refit LASTRC=[$WOOLDID_LASTRC]"
est restore GOOD2
cap noi wooldid
dis "@@G8 bare-after-refit-and-restore rc=" _rc

*================= rc masking / cleanup skip under a stale LEVEL0 ===========
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
* baseline: what rc does a bad cluster() varname give on a clean session?
cap noi wooldid y uid t ttre, att robust cluster(nosuchvar)
dis "@@H0 baseline-badcluster rc=" _rc " LASTRC=[$WOOLDID_LASTRC]"
qui est dir
dis "@@H0b estdir=[`r(names)'] curframe:"
qui frame
dis "@@H0c frame=[`r(currentframe)']"

* poison
set level 99.9
cap noi wooldid y uid t ttre, att robust level(95)
dis "@@H1 poison rc=" _rc " LEVEL0=[$WOOLDID_LEVEL0] clevel=" c(level)
* now the same bad-cluster call
cap noi wooldid y uid t ttre, att robust cluster(nosuchvar)
dis "@@H2 poisoned-badcluster rc=" _rc " LASTRC=[$WOOLDID_LASTRC] LEVEL0=[$WOOLDID_LEVEL0]"
qui est dir
dis "@@H2b estdir=[`r(names)']"
qui frame
dis "@@H2c frame=[`r(currentframe)'] emptycells=[`c(emptycells)']"
qui frame dir
dis "@@H2d framedir=[`r(frames)']"

macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95
dis "S5_DONE"
