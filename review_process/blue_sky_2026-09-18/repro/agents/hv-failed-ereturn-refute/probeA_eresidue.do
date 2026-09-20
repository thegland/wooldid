*! probeA_eresidue.do -- SKEPTIC's own repro of HV-EPOST-01 / HV-EPOST-02
*! Self-contained: generates its own panel, never touches the project test suite.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probeA_eresidue.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771113
set sortseed 771113

dis "=============== BUILD ASSERTION ==============="
which wooldid

*---------------------------------------------------------------------
* 0.  My own panel.  40 counties in 8 regions x 8 years, 2012-2019.
*     Policy adopted AT THE REGION LEVEL (realistic admin shape):
*     regions 1-2 in 2015, regions 3-4 in 2017, regions 5-8 never.
*---------------------------------------------------------------------
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
gen long regyear = region*10000 + year
gen double x1 = rnormal()
gen double dose = runiform()*2
gen byte post = (year >= onset) & !missing(onset)
gen double y  = 2 + ai + 0.12*(year-2011) + 1.5*post + 0.4*x1 + rnormal()
* a count-ish outcome that ALSO contains negatives: a realistic poisson mis-spec
gen double ynet = y - 3
xtset unit year
compress
save "panelA.dta", replace

*---------------------------------------------------------------------
* 1.  A GOOD run.  Record its e() fingerprint.
*---------------------------------------------------------------------
dis "=============== SEC1 GOOD RUN ==============="
wooldid y unit year onset, att
dis "SEC1_CMD=" e(cmd)
dis "SEC1_CMDLINE=" e(cmdline)
dis "SEC1_KB=" colsof(e(b))
dis "SEC1_N=" e(N)
local goodb = e(b)[1,1]
dis "SEC1_B11=" %21x `goodb'
est store good
dis "SEC1_LASTRC=[${WOOLDID_LASTRC}]"

*---------------------------------------------------------------------
* 2.  Trigger battery.  Each trigger runs from the SAME starting e()
*     (the good run), restored beforehand.  Report what e() holds after.
*---------------------------------------------------------------------
cap prog drop epostreport
prog def epostreport
  syntax , tag(string) rc(string)
  dis "`tag'_RC=`rc'"
  dis "`tag'_ECMD=" e(cmd)
  local kb = .
  cap local kb = colsof(e(b))
  dis "`tag'_KB=`kb'"
  dis "`tag'_EN=" e(N)
  local ns = .
  cap count if e(sample)
  cap local ns = r(N)
  dis "`tag'_ESAMP=`ns'"
  local b11 = .
  cap local b11 = e(b)[1,1]
  dis "`tag'_B11=" %21x `b11'
  local cn
  cap local cn : colfullnames e(b)
  dis "`tag'_COL1=" word("`cn'", 1)
  local ed
  cap qui est dir
  cap local ed "`r(names)'"
  dis "`tag'_ESTDIR=[`ed']"
  qui frame
  dis "`tag'_FRAME=" r(currentframe)
  dis "`tag'_LASTRC=[${WOOLDID_LASTRC}]"
end

* ---- T1: region-by-year FE where the POLICY VARIES AT THE REGION LEVEL.
*      Entirely plausible spec (state-by-year FE) that absorbs the treatment.
qui est restore good
dis "=============== T1 fe(regyear) ==============="
cap noi wooldid y unit year onset, att fe(regyear)
epostreport , tag(T1) rc(`=_rc')

* ---- T2: poisson on an outcome that contains negative values.
qui est restore good
dis "=============== T2 poisson on negative y ==============="
cap noi wooldid ynet unit year onset, att poisson
epostreport , tag(T2) rc(`=_rc')

* ---- T3: subgroups() on a variable whose levels include wholly-untreated groups.
qui est restore good
dis "=============== T3 subgroups(region) ==============="
cap noi wooldid y unit year onset, att subgroups(region)
epostreport , tag(T3) rc(`=_rc')

* ---- T4: the recon-suggested trigger, on MY data.
qui est restore good
dis "=============== T4 contreat+untreatedZero ==============="
cap noi wooldid y unit year onset, att contreat(dose) contreatcontrols(untreatedZero)
epostreport , tag(T4) rc(`=_rc')

* ---- T5: a PARSE-stage failure (control).
qui est restore good
dis "=============== T5 parse typo ==============="
cap noi wooldid y unit year onset, attt
epostreport , tag(T5) rc(`=_rc')

*---------------------------------------------------------------------
* 3.  Consequences.  Take the first trigger that actually replaced e().
*---------------------------------------------------------------------
dis "=============== SEC3 CONSEQUENCES ==============="
qui est restore good
qui cap noi wooldid y unit year onset, att fe(regyear)
dis "SEC3_RC=" _rc
* a scripted harvest that guards with capture but reads e() unconditionally
dis "SEC3_HARVEST_B=" %21x e(b)[1,1]
dis "SEC3_HARVEST_N=" e(N)
dis "SEC3_TRUTH_B=" %21x `goodb'
* does lincom on an internal regressor succeed?
local cn : colfullnames e(b)
local c1 = word("`cn'", 1)
cap noi lincom _b[`c1']
dis "SEC3_LINCOM_RC=" _rc
cap noi lincom _b[main_att]
dis "SEC3_LINCOMATT_RC=" _rc
cap noi predict double zz, xb
dis "SEC3_PREDICT_RC=" _rc
cap drop zz

*---------------------------------------------------------------------
* 4.  HV-EPOST-02: bare `wooldid' replay after each failure class.
*---------------------------------------------------------------------
dis "=============== SEC4 REPLAY AFTER POST-REGRESSION FAILURE ==============="
* e() currently = the wrecked one from SEC3
cap noi wooldid
dis "SEC4_POSTREG_REPLAY_RC=" _rc

dis "=============== SEC4b REPLAY AFTER PARSE FAILURE ==============="
qui est restore good
qui cap wooldid y unit year onset, attt
dis "SEC4b_TRIGRC=" _rc
cap noi wooldid
dis "SEC4b_PARSE_REPLAY_RC=" _rc

dis "=============== SEC4c REPLAY AFTER est restore OF THE GOOD RUN ==============="
qui est restore good
qui cap wooldid y unit year onset, att fe(regyear)
dis "SEC4c_TRIGRC=" _rc
qui est restore good
cap noi wooldid
dis "SEC4c_RESTORED_REPLAY_RC=" _rc

dis "=============== SEC4d BARE wooldid IN A VIRGIN e() STATE ==============="
* what does bare wooldid say when there never was a wooldid result?
macro drop WOOLDID_LASTRC
qui reg y x1
cap noi wooldid
dis "SEC4d_VIRGIN_RC=" _rc
dis "SEC4d_ECMD=" e(cmd)

* and: would gating on WOOLDID_LASTRC alone (the prober's proposed remedy) be right?
* here e(cmd)=regress and the marker is set -> a LASTRC-only gate would print
* "the stored results predate that failed call" about a regress model.
global WOOLDID_LASTRC 322
dis "SEC4e_ECMD=" e(cmd) " MARKER=[${WOOLDID_LASTRC}]"
macro drop WOOLDID_LASTRC

*---------------------------------------------------------------------
* 5.  My OWN baseline for the "expected" half: does Stata protect e()?
*---------------------------------------------------------------------
dis "=============== SEC5 BASELINE ==============="
sysuse auto, clear
qui regress price mpg
dis "SEC5_START=" e(cmdline)

* (a) a built-in command that fails at PARSE
cap noi regress price nosuchvar
dis "SEC5a_RC=" _rc
dis "SEC5a_CMDLINE=" e(cmdline)

* (b) an eclass ado that estimates then errors -- wooldid's shape
qui regress price mpg
cap prog drop mybad
prog def mybad, eclass
  qui regress price weight foreign
  error 322
end
cap noi mybad
dis "SEC5b_RC=" _rc
dis "SEC5b_CMDLINE=" e(cmdline)

* (c) the same ado WITH the standard idiom
qui regress price mpg
cap prog drop mygood
prog def mygood, eclass
  tempname h
  _estimates hold `h', restore nullok
  qui regress price weight foreign
  error 322
end
cap noi mygood
dis "SEC5c_RC=" _rc
dis "SEC5c_CMDLINE=" e(cmdline)

* (d) how common is the idiom in Stata's own base ado?
qui regress price mpg
dis "SEC5_DONE"

dis "=============== PROBE A COMPLETE ==============="
