* p2_ctw.do -- SELF-CONTAINED.
* Does contreatwithin + contreat() silently disable the two data-integrity screens
* that sit between the CS-98 safety rewrite (ado:4829) and the old rewrite site (ado:5169)?
*   ado:4959  "i must consist of positive, non-zero integers"      (error 411)
*   ado:4973  float-above-2^24 silent-merge caution               (BS-34)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
dis "BUILDCHECK_DONE"

*-------------------------------------------------------------------
* A. FLOAT-ABOVE-2^24 CAUTION
*    i is an 8-digit administrative code stored as float (the default
*    storage type of `gen'), so adjacent ids collide silently.
*-------------------------------------------------------------------
capture program drop mkfloat
program define mkfloat
  clear
  set obs 60
  gen long unit = _n
  gen double ttre = .
  replace ttre = 5 if unit <= 15
  replace ttre = 7 if unit > 15 & unit <= 30
  replace ttre = 9 if unit > 30 & unit <= 45
  gen double ufe = rnormal()
  gen double dose = 0.5 + runiform()
  expand 12
  bysort unit: gen double t = _n
  gen byte post = (t >= ttre) & ttre < .
  gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
  * 8-digit code held in FLOAT -- exactly the BS-34 hazard the caution exists for
  gen float i = 20000000 + unit
end

dis "=== A1: float i, plain binary, NO contreat  (caution expected) ==="
mkfloat
cap noi wooldid y i t ttre, att
dis "A1_RC=" _rc

dis "=== A2: float i, contreat, NO contreatwithin (caution expected) ==="
mkfloat
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "A2_RC=" _rc

dis "=== A3: float i, contreat + contreatwithin (caution SUPPRESSED?) ==="
mkfloat
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "A3_RC=" _rc

dis "=== A4: same as A3 but user EXPLICITLY asks for safety(on) ==="
mkfloat
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin safety(on)
dis "A4_RC=" _rc

*-------------------------------------------------------------------
* B. i-MUST-BE-POSITIVE-INTEGER SCREEN  (ado:4959, error 411)
*    one treated cohort carries i == 0.  _i is built as
*    "i if treated, 0 otherwise", so i==0 collides with the control code.
*-------------------------------------------------------------------
capture program drop mkzero
program define mkzero
  args shift
  clear
  set obs 60
  gen long unit = _n
  gen double ttre = .
  replace ttre = 5 if unit <= 15
  replace ttre = 7 if unit > 15 & unit <= 30
  replace ttre = 9 if unit > 30 & unit <= 45
  gen double ufe = rnormal()
  gen double dose = 0.5 + runiform()
  expand 12
  bysort unit: gen double t = _n
  gen byte post = (t >= ttre) & ttre < .
  gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
  * `shift' = 0  ->  ids 0..59 : one TREATED unit carries i == 0
  * `shift' = 1  ->  ids 1..60 : same panel, same cohort structure, no zero id
  gen double i = unit - 1 + `shift'
end

dis "=== B1: i contains 0, contreat, NO contreatwithin (error 411 expected) ==="
mkzero 0
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "B1_RC=" _rc

dis "=== B2: i contains 0, contreat + contreatwithin (screen skipped?) ==="
mkzero 0
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B2_RC=" _rc
cap matrix B2 = e(b)
cap dis "B2_ATT=" %21x el(e(b),1,1)
cap dis "B2_NOBS=" e(nobs)
cap dis "B2_NCOHORT=" e(ncohorts)

dis "=== B3: SAME PANEL, ids shifted to 1..60, contreat + contreatwithin ==="
mkzero 1
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B3_RC=" _rc
cap matrix B3 = e(b)
cap dis "B3_ATT=" %21x el(e(b),1,1)
cap dis "B3_NOBS=" e(nobs)
cap dis "B3_NCOHORT=" e(ncohorts)

dis "=== B4: i contains 0, BINARY (no contreat), safety default -- control ==="
mkzero 0
cap noi wooldid y i t ttre, att
dis "B4_RC=" _rc

dis "=== B5: i NON-INTEGER, contreat + contreatwithin ==="
mkzero 1
qui replace i = i + 0.5 if unit == 3
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B5_RC=" _rc

dis "=== B6: i NON-INTEGER, contreat, NO contreatwithin -- control ==="
mkzero 1
qui replace i = i + 0.5 if unit == 3
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "B6_RC=" _rc

dis "ALLDONE_P2"
