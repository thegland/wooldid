*=====================================================================
* probe04.do -- HV-FAILED-ERETURN: residue, replay-pair, peers
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
gen double ypos = exp(y/4)
qui ds
local VARS0 "`r(varlist)'"
dis "VARS0=`VARS0'"

*--- 1. the discriminating replay pair --------------------------------
dis "=== P4_1 replay after PARSE failure ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, attt
dis "P4_parse_rc=" _rc
dis "P4_parse_cmd=" "`e(cmd)'"
cap noisily wooldid
dis "P4_replay_after_parse_rc=" _rc

dis "=== P4_1b replay after POST-REGRESSION failure ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
dis "P4_post_rc=" _rc
dis "P4_post_cmd=" "`e(cmd)'"
cap noisily wooldid
dis "P4_replay_after_post_rc=" _rc

*--- 2. residue: estimates dir, variables, emptycells, frame ----------
dis "=== P4_2 residue after each failure class ==="
foreach spec in "contreat(dosep) contreatcontrols(untreatedZero)" "ccc_absorb(x1)" {
  qui wooldid y unit year onset, att
  qui est store base_ok
  cap noisily wooldid y unit year onset, att `spec'
  dis "P4_2|spec=`spec'|rc=" _rc "|cmd=`e(cmd)'|emptycells=`c(emptycells)'"
  qui frame
  dis "P4_2|frame=`r(currentframe)'"
  qui ds
  local V1 "`r(varlist)'"
  local extra : list V1 - VARS0
  dis "P4_2|extravars=`extra'"
  qui est dir
  dis "P4_2|estdir=`r(names)'"
  qui est drop base_ok
}

dis "=== P4_2c poisson residue ==="
qui wooldid y unit year onset, att
cap noisily wooldid ypos unit year onset, att poisson ccc_absorb(x1)
dis "P4_2c|rc=" _rc "|cmd=`e(cmd)'"
qui ds
local V1 "`r(varlist)'"
local extra : list V1 - VARS0
dis "P4_2c|extravars=`extra'"
qui est dir
dis "P4_2c|estdir=`r(names)'"

*--- 3. is the NEXT good run harmed? ----------------------------------
dis "=== P4_3 clean-vs-after-failure reproducibility ==="
clear all
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
clear
set obs 60
gen int unit = _n
gen byte grp  = mod(_n,5)
gen double onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i = rnormal()
gen double dose0 = 0.5 + runiform()
expand 10
bys unit: gen int year = 1999 + _n
gen byte post = !missing(onset) & year >= onset
gen double x1 = rnormal()
gen double tau = 0
replace tau = 0.8*(year-onset+1) if post
gen double y = a_i + 0.10*(year-2000) + tau + 0.4*x1 + 0.30*rnormal()
gen double dosep = dose0*post
qui wooldid y unit year onset, att
dis "P4_3_clean=" %21x _b[main_att]
cap wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
qui wooldid y unit year onset, att
dis "P4_3_after=" %21x _b[main_att]

*--- 4. the r461 harvest (small kb, most deceptive) -------------------
dis "=== P4_4 r461 harvest ==="
qui wooldid y unit year onset, att
dis "P4_4_good_N=" e(N) "  good_b11=" e(b)[1,1]
cap wooldid y unit year onset, att ccc_absorb(x1)
dis "P4_4_rc=" _rc
dis "P4_4_fail_N=" e(N) "  fail_b11=" e(b)[1,1] "  kb=" colsof(e(b))
cap count if e(sample)
dis "P4_4_nsamp=" r(N)

*--- 5. peers that fail after internal work ---------------------------
dis "=== P4_5 peers ==="
regress y x1
dis "P4_5_base=" "`e(cmdline)'"
cap noisily did_imputation y unit year onset, horizons(3) pretrends(200)
dis "P4_5_didimp_rc=" _rc
dis "P4_5_didimp_cmd=" "`e(cmd)'"
dis "P4_5_didimp_cmdline=" substr(`"`e(cmdline)'"',1,70)

regress y x1
cap noisily jwdid y, ivar(unit) tvar(year) gvar(onset) group(nosuchvar)
dis "P4_5_jwdid_rc=" _rc
dis "P4_5_jwdid_cmd=" "`e(cmd)'"
dis "P4_5_jwdid_cmdline=" substr(`"`e(cmdline)'"',1,70)

dis "PROBE04 DONE"
