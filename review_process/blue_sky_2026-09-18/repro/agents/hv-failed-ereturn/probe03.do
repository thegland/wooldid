*=====================================================================
* probe03.do -- HV-FAILED-ERETURN: consequence surface
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 60
gen int unit = _n
gen byte grp  = mod(_n,5)
gen double onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i = rnormal()
gen byte state = 1 + mod(_n,10)
gen double dose0 = 0.5 + runiform()
expand 10
bys unit: gen int year = 1999 + _n
gen byte post = !missing(onset) & year >= onset
gen double x1 = rnormal()
gen double tau = 0
replace tau = 0.8*(year-onset+1) if post
gen double y = a_i + 0.10*(year-2000) + tau + 0.4*x1 + 0.30*rnormal()
gen double dosep = dose0*post

*---------------- 1. GOOD RUN, stored -------------------------------
dis "=== STEP1 good run ==="
wooldid y unit year onset, att
est store goodwd
dis "S1_att=" _b[main_att]
dis "S1_cmd=" "`e(cmd)'"

*---------------- 2. FAILING RUN ------------------------------------
dis "=== STEP2 failing run (recon repro) ==="
cap noisily wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
dis "S2_rc=" _rc
dis "S2_cmd=" "`e(cmd)'"
dis "S2_kb=" colsof(e(b))
dis "S2_N=" e(N)

*---------------- 3. what the user sees next ------------------------
dis "=== STEP3a replay: bare wooldid ==="
cap noisily wooldid
dis "S3a_rc=" _rc
dis "S3a_LASTRC=$WOOLDID_LASTRC"

dis "=== STEP3b predict ==="
cap noisily predict double xbhat
dis "S3b_rc=" _rc

dis "=== STEP3c test/lincom on an internal name ==="
local nm : colnames e(b)
local n1 : word 2 of `nm'
dis "S3c_firstname=`n1'"
cap noisily lincom _b[`n1']
dis "S3c_rc=" _rc

dis "=== STEP3d scripted harvest (the silent path) ==="
* a completely ordinary loop: run a spec, capture, harvest e()
matrix RES = J(1,3,.)
cap wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
matrix RES[1,1] = e(b)[1,1]
matrix RES[1,2] = e(N)
matrix RES[1,3] = _rc
dis "S3d_harvested_b11=" RES[1,1]
dis "S3d_harvested_N="   RES[1,2]
matrix list RES

dis "=== STEP3e estimates table / esttab ==="
cap est store failwd
dis "S3e_store_rc=" _rc
cap noisily estimates table goodwd failwd, keep(main_att)
dis "S3e_rc=" _rc

dis "=== STEP3f can the good run be recovered? ==="
cap noisily est restore goodwd
dis "S3f_rc=" _rc
dis "S3f_cmd=" "`e(cmd)'"
cap dis "S3f_att=" _b[main_att]

dis "=== STEP3g does wooldid replay work after est restore? ==="
cap noisily wooldid
dis "S3g_rc=" _rc

*---------------- 4. PEER COMPARISON --------------------------------
dis "=== STEP4 peers: does a failed peer clobber e()? ==="
* good baseline
regress y x1
dis "S4_base_cmdline=" "`e(cmdline)'"
* jwdid: a spec that fails after its internal regression
cap noisily jwdid y x1, ivar(unit) tvar(year) gvar(onset)
dis "S4_jwdid_ok_rc=" _rc
dis "S4_jwdid_cmd=" "`e(cmd)'"
regress y x1
cap noisily jwdid y, ivar(unit) tvar(year) gvar(onset) cluster(nosuchvar)
dis "S4_jwdidfail_rc=" _rc
dis "S4_jwdidfail_cmd=" "`e(cmd)'"
dis "S4_jwdidfail_cmdline=" substr(`"`e(cmdline)'"',1,70)

regress y x1
cap noisily csdid y, ivar(unit) time(year) gvar(onset) nosuchoption
dis "S4_csdidfail_rc=" _rc
dis "S4_csdidfail_cmd=" "`e(cmd)'"

dis "PROBE03 DONE"
