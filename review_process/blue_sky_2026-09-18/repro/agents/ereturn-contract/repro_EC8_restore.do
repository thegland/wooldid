* EC-8: after any failed wooldid, `estimates restore' returns a valid result set that
* wooldid then refuses to replay, saying the stored results "predate that failed call".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
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
gen double y = ufe + 0.25*t + 1.2*post + rnormal()

wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2)
estimates store GOODM
dis "stored: e(cmd)=`e(cmd)'  K=" colsof(e(b))
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "the failed call: rc=" _rc "   WOOLDID_LASTRC=${WOOLDID_LASTRC}"
estimates restore GOODM
dis "after restore: e(cmd)=`e(cmd)'  K=" colsof(e(b)) "   <- a VALID result set"
dis "WOOLDID_LASTRC is still ${WOOLDID_LASTRC}"
cap noi wooldid
dis "bare wooldid -> rc " _rc "   <- refuses, and its message is untrue"
dis ""
dis "--- and the macro-drop defeat (recon F4), for contrast ---"
macro drop WOOLDID_LASTRC
cap noi wooldid
dis "after macro drop WOOLDID_LASTRC: bare wooldid -> rc " _rc
