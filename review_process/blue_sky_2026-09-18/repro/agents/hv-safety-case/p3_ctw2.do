* p3_ctw2.do -- SELF-CONTAINED.
* Clean demonstration: contreatwithin + contreat() silently suppresses
*   (A) the float-above-2^24 caution (ado:4973) on a panel that RUNS TO COMPLETION
*   (B) the "i must be positive non-zero integers" refusal (ado:4959), and the
*       resulting run returns a DIFFERENT number than the same panel with ids shifted by 1.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
dis "BUILDCHECK_DONE"

*===================================================================
* A. float 8-digit id, onset assigned FROM the rounded id so the
*    onset-constancy screen (which gates on the USER's safety value)
*    does not fire.  The float caution is then the only warning there is.
*===================================================================
capture program drop mkA
program define mkA
  clear
  set obs 60
  gen long unit = _n
  * 8-digit administrative code, default `gen' storage = float
  gen float i = 20000000 + unit
  gen double ttre = .
  replace ttre = 5  if i <= 20000020
  replace ttre = 8  if i > 20000020 & i <= 20000040
  gen double ufe = rnormal()
  gen double dose = 0.5 + runiform()
  expand 12
  bysort unit: gen double t = _n
  gen byte post = (t >= ttre) & ttre < .
  gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
end

dis "=== A1 float i, contreat, NO contreatwithin ==="
mkA
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "A1_RC=" _rc
cap noi dis "A1_NCOHORTS=" e(ncohorts)

dis "=== A2 float i, contreat + contreatwithin ==="
mkA
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "A2_RC=" _rc
cap noi dis "A2_NCOHORTS=" e(ncohorts)

dis "=== A3 float i, contreat + contreatwithin + EXPLICIT safety(on) ==="
mkA
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin safety(on)
dis "A3_RC=" _rc

dis "=== A4 DOUBLE i (the fix the caution asks for), contreatwithin ==="
mkA
recast double i
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "A4_RC=" _rc
cap noi dis "A4_NCOHORTS=" e(ncohorts)

*===================================================================
* B. i == 0 on one treated cohort.  _i is built as "i if treated, 0
*    otherwise", so i==0 aliases the control code.
*===================================================================
capture program drop mkB
program define mkB
  args shift
  clear
  set obs 60
  gen long unit = _n
  gen double ttre = .
  replace ttre = 5 if unit <= 15
  replace ttre = 8 if unit > 15 & unit <= 30
  gen double ufe = rnormal()
  gen double dose = 0.5 + runiform()
  expand 12
  bysort unit: gen double t = _n
  gen byte post = (t >= ttre) & ttre < .
  gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
  * shift 0 -> ids 0..59 (unit 1, a TREATED unit, carries i == 0)
  * shift 1 -> ids 1..60 (identical panel, identical cohorts, no zero id)
  gen double i = unit - 1 + `shift'
end

dis "=== B1 i has a 0, contreat, NO contreatwithin (expect error 411) ==="
mkB 0
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic)
dis "B1_RC=" _rc

dis "=== B2 i has a 0, contreat + contreatwithin ==="
mkB 0
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B2_RC=" _rc
cap noi dis "B2_ATT_HEX=" %21x el(e(b),1,1)
cap noi dis "B2_ATT=" el(e(b),1,1)
cap noi dis "B2_NCOHORTS=" e(ncohorts)
cap noi dis "B2_NOBS=" e(nobs)

dis "=== B3 same panel, ids 1..60, contreat + contreatwithin ==="
mkB 1
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
dis "B3_RC=" _rc
cap noi dis "B3_ATT_HEX=" %21x el(e(b),1,1)
cap noi dis "B3_ATT=" el(e(b),1,1)
cap noi dis "B3_NCOHORTS=" e(ncohorts)
cap noi dis "B3_NOBS=" e(nobs)

*===================================================================
* C. does anything at runtime say safety was relaxed?  dump e() keys.
*===================================================================
dis "=== C: e() scan for any safety marker after a contreatwithin run ==="
mkB 1
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) contreatwithin
cap noi dis "C_CMDLINE=" e(cmdline)
foreach k in wdid_safety safety wdid_safetystate {
  cap noi dis "C_`k'=[" e(`k') "]"
}

dis "ALLDONE_P3"
