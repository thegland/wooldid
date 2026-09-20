* repro_ctw_floatcaution.do  -- SELF-CONTAINED.  No other file needed.
* HV-CTW-FLOAT: contreatwithin + contreat() silently suppresses the BS-34
* float-above-2^24 caution (ado:4973), because the CS-98 rewrite at ado:4829
* sets the WORKING safety value to off ABOVE that check.
* sthlp:73 promises the caution unconditionally; sthlp:1195 puts "every
* data-integrity assert" on the user's-own-value side of the rewrite.
* Both are false on this build.  Pre-CS-98 (git e26b3dd^) the check sat ABOVE
* the rewrite and did fire, so this is a regression dated 2026-09-09.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case"
which wooldid

* an 8-digit administrative code held in FLOAT -- the default storage type of `gen'.
* onset is assigned FROM the rounded code so the onset-constancy screen (which
* correctly gates on the user's own safety value) does NOT fire and mask the point.
clear
set obs 60
gen long unit = _n
gen float  ibig   = 20000000 + unit
gen double ibig_d = 20000000 + unit
gen double ttre = .
replace ttre = 5 if ibig <= 20000020
replace ttre = 8 if ibig > 20000020 & ibig <= 20000040
gen double ufe = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort unit: gen double t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
save "`WD'/repro_float.dta", replace
qui tab ibig
dis "DISTINCT_IDS_AS_FLOAT=" r(r)     // 31 : 29 of the 60 units were silently merged
qui tab ibig_d
dis "DISTINCT_IDS_AS_DOUBLE=" r(r)    // 60

dis "=== LEG 1: NO contreatwithin -- the caution PRINTS ==="
use "`WD'/repro_float.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i)
dis "LEG1_RC=" _rc

dis "=== LEG 2: + contreatwithin -- the caution is GONE ==="
use "`WD'/repro_float.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) contreatwithin
dis "LEG2_RC=" _rc
cap noi dis "LEG2_SE=" sqrt(el(e(V),1,1))

dis "=== LEG 3: + contreatwithin + EXPLICIT safety(on) -- still gone ==="
use "`WD'/repro_float.dta", clear
gen float i = ibig
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) contreatwithin safety(on)
dis "LEG3_RC=" _rc

dis "=== LEG 4: the panel the user MEANT (double ids), contreatwithin ==="
use "`WD'/repro_float.dta", clear
gen double i = ibig_d
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) contreatwithin
dis "LEG4_RC=" _rc
cap noi dis "LEG4_SE=" sqrt(el(e(V),1,1))

dis "EXPECTED: LEG2 and LEG3 print the same float caution LEG1 prints."
dis "OBSERVED: LEG1 prints it; LEG2 and LEG3 print nothing, and the run completes"
dis "          on 31 silently-merged clusters instead of 60.  LEG2 SE .07147848"
dis "          vs LEG4 SE .0684769 -- 4.4% apart, with no warning in LEG2."
dis "ALLDONE"
