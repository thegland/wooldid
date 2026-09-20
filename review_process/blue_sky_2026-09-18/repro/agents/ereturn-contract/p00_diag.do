clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid

* variant A: recon's proven shape (80 units, ttre 4/6/8, t=1..10)
clear
set obs 80
gen long id = _n
gen byte grp = mod(id,4)
gen double ttre = .
replace ttre = 4 if grp==1
replace ttre = 6 if grp==2
replace ttre = 8 if grp==3
gen double ufe = rnormal()
expand 10
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .3*t + 1.5*post + rnormal()
cap noi wooldid y id t ttre, att
dis "VA_RC=" _rc

* variant B: mine, 120 units, ttre 5/7/9/11, t=1..12
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
replace ttre = 11 if grp==4
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .25*t + 1.2*post + rnormal()
qui tab ttre, missing
cap noi wooldid y id t ttre, att
dis "VB_RC=" _rc

* variant B2: same but drop the latest cohort (ttre 11 -> never treated)
replace ttre = . if ttre==11
cap noi wooldid y id t ttre, att
dis "VB2_RC=" _rc

* variant C: mine, but t=1..12 with ttre 5/7/9 only
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .25*t + 1.2*post + rnormal()
cap noi wooldid y id t ttre, att
dis "VC_RC=" _rc
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3)
dis "VC_ES_RC=" _rc

* variant D: B with an explicit late-cohort truncation of the window
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
replace ttre = 11 if grp==4
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .25*t + 1.2*post + rnormal()
cap noi wooldid y id t ttre, att oostreatedcontrols(3)
dis "VD_RC=" _rc

dis "DIAG_DONE"
