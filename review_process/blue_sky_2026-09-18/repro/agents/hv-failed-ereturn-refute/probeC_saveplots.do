*! probeC_saveplots.do -- SKEPTIC's own finding: saveplots() has no pre-flight
*! existence check, so RERUNNING A DO-FILE without `replace' destroys the run at
*! the very last step -- exactly the failure mode the author fixed for bootsave()
*! on 2026-09-05 (ado:6276-6278) and permsave() (ado:5866-5877).
*!   C1: bootstrap(standard) reps(12) -- the whole replication loop is lost.
*!   C2: analytic makecfxplots -- the r(602) fires at ado:2494, BEFORE wooldid's
*!       own ereturn post at ado:4051, so e() is left holding reghdfe's internal
*!       saturated model.  A file-already-exists rerun is thus a HV-EPOST-01 trigger.
*! Self-contained.  Run:
*!   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probeC_saveplots.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771113
set sortseed 771113
which wooldid

local fl : dir "." files "pcstem*"
foreach f of local fl {
  erase "`f'"
}
local fl : dir "." files "pdstem*"
foreach f of local fl {
  erase "`f'"
}

clear
set obs 40
gen int unit   = _n
gen byte region = ceil(unit/5)
gen int onset = .
replace onset = 2015 if region <= 2
replace onset = 2017 if region > 2 & region <= 4
gen double ai = rnormal()
expand 8
bysort unit: gen int year = 2011 + _n
gen double x1 = rnormal()
gen double dose = runiform()*2
gen byte post = (year >= onset) & !missing(onset)
gen double y  = 2 + ai + 0.12*(year-2011) + 1.5*post*(0.5+dose) + 0.4*x1 + rnormal()
xtset unit year

*=====================================================================
* C1 -- bootstrap: the whole replication loop is lost on a rerun
*=====================================================================
dis "=============== C1 RUN 1 (bootstrap 12 reps, writes pcstem_*) ==============="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(pcstem) bootstrap(standard) reps(12) rseed(4242)
dis "C1_RUN1_RC=" _rc
dis "C1_RUN1_ECMD=" e(cmd) "  BOOTM=" e(boot_method) "  KB=" colsof(e(b))
local fl : dir "." files "pcstem*"
dis `"C1_FILES=`fl'"'

dis "=============== C1 RUN 2 (identical line, files exist) ==============="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(pcstem) bootstrap(standard) reps(12) rseed(4242)
dis "C1_RUN2_RC=" _rc
dis "C1_RUN2_ECMD=" e(cmd)
dis "C1_RUN2_BOOTM=" e(boot_method)
dis "C1_RUN2_BOOTREPS=" e(boot_reps)
local kb = .
cap local kb = colsof(e(b))
dis "C1_RUN2_KB=`kb'"
dis "C1_RUN2_B11=" %21x e(b)[1,1]
dis "C1_RUN2_LASTRC=[${WOOLDID_LASTRC}]"
dis "=== C1 replay of that result ==="
cap noi wooldid
dis "C1_REPLAY_RC=" _rc

*=====================================================================
* C2 -- analytic makecfxplots: the r(602) fires BEFORE wooldid's own post
*=====================================================================
dis "=============== C2 RUN 1 (analytic cfx plots, writes pdstem_*) ==============="
cap noi wooldid y unit year onset, att contreat(dose) makecfxplots saveplots(pdstem)
dis "C2_RUN1_RC=" _rc
dis "C2_RUN1_ECMD=" e(cmd) "  KB=" colsof(e(b))
local c2good = e(b)[1,1]
dis "C2_RUN1_B11=" %21x `c2good'
est store c2good
local fl : dir "." files "pdstem*"
dis `"C2_FILES=`fl'"'

dis "=============== C2 RUN 2 (identical line, files exist) ==============="
cap noi wooldid y unit year onset, att contreat(dose) makecfxplots saveplots(pdstem)
dis "C2_RUN2_RC=" _rc
dis "C2_RUN2_ECMD=" e(cmd)
local kb = .
cap local kb = colsof(e(b))
dis "C2_RUN2_KB=`kb'"
dis "C2_RUN2_B11=" %21x e(b)[1,1]
dis "C2_RUN2_N=" e(N)
local cn
cap local cn : colfullnames e(b)
dis "C2_RUN2_COL1=" word("`cn'", 1)
dis "C2_RUN2_LASTRC=[${WOOLDID_LASTRC}]"
dis "C2_TRUTH_B11=" %21x `c2good'
dis "=== C2 replay ==="
cap noi wooldid
dis "C2_REPLAY_RC=" _rc

*=====================================================================
* C3 -- control: the same rerun WITH replace
*=====================================================================
dis "=============== C3 CONTROL (replace) ==============="
cap noi wooldid y unit year onset, att contreat(dose) makecfxplots saveplots(pdstem) replace
dis "C3_RC=" _rc
dis "C3_ECMD=" e(cmd)
dis "C3_LASTRC=[${WOOLDID_LASTRC}]"

dis "=============== PROBE C COMPLETE ==============="
