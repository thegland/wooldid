clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* tiny panel: 24 units in 4 cohort groups (3 treated + never), 8 periods
set obs 24
gen id = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 4 if grp==1
replace ttre = 5 if grp==2
replace ttre = 6 if grp==3
gen double ui = rnormal()
gen double dose = 0.5 + runiform()
expand 8
bysort id: gen t = _n
gen byte post = (ttre<. & t>=ttre)
gen double rel = cond(post==1, t-ttre, 0)
gen double y = 2 + ui + 0.1*t + rnormal() + post*(1+0.3*rel)
gen double w = 1 + runiform()
* cohort-level id (coarse) so runs are fast
egen cid = group(grp)
gen double negw = w
replace negw = -1 in 1
gen idneg = id
replace idneg = -id
save tiny.dta, replace

*** M1: safety(Off) capitalization -- is the value validated before it is consulted?
di "MECH_BEGIN M1a"
cap noisily wooldid y cid t ttre, att_i unconditionalse robust safety(Off)
di "MECH_RC=" _rc
di "MECH_END M1a"

di "MECH_BEGIN M1b"
cap noisily wooldid y cid t ttre, att_i unconditionalse robust safety(off)
di "MECH_RC=" _rc
di "MECH_END M1b"

di "MECH_BEGIN M1c"
cap noisily wooldid y cid t ttre, att safety(Off)
di "MECH_RC=" _rc
di "MECH_END M1c"

*** M2: does contreatwithin+contreat silently disable the i-positive-integer screen?
di "MECH_BEGIN M2a_no_ctw"
cap noisily wooldid y idneg t ttre, att fastnose
di "MECH_RC=" _rc
di "MECH_END M2a_no_ctw"

di "MECH_BEGIN M2b_ctw"
cap noisily wooldid y idneg t ttre, att fastnose contreat(dose) contreatcontrols(byTime) contreatwithin
di "MECH_RC=" _rc
di "MECH_END M2b_ctw"

*** M2c: negative weights -- the sibling that WAS converted to safetyinput
di "MECH_BEGIN M2c_negw"
cap noisily wooldid y cid t ttre [aw=negw], att fastnose contreat(dose) contreatcontrols(byTime) contreatwithin
di "MECH_RC=" _rc
di "MECH_END M2c_negw"

*** M3: caution printed BEFORE the aborting error (permctl low-reps)
di "MECH_BEGIN M3a"
cap noisily wooldid y cid t ttre, att permute(shuffle) permreps(25)
di "MECH_RC=" _rc
di "MECH_END M3a"

di "MECH_BEGIN M3b"
cap noisily wooldid y cid t ttre, att permute(shuffle) permreps(25) cluster(cid) unconditionalse
di "MECH_RC=" _rc
di "MECH_END M3b"

di "MECH_BEGIN M3c"
cap noisily wooldid y cid t ttre, att permute(shuffle) permreps(25) cluster(cid) subgroups(grp) permonsetcohorts
di "MECH_RC=" _rc
di "MECH_END M3c"

*** M4: bare error 459 from the byCohort dose-slope screen -- is it explained?
di "MECH_BEGIN M4"
gen double dosec = grp
cap noisily wooldid y cid t ttre, att fastnose contreat(dosec) contreatcontrols(byCohort)
di "MECH_RC=" _rc
di "MECH_END M4"

*** M5: a message naming an internal local the user cannot type
di "MECH_BEGIN M5"
cap noisily wooldid y cid t ttre, att contreat(dose) contreatcontrols(byTime) lattice(all 1)
di "MECH_RC=" _rc
di "MECH_END M5"

*** M6: bare replay refusal / replay contract
di "MECH_BEGIN M6"
cap noisily wooldid y cid t ttre, att nonsenseoption
di "MECH_RC=" _rc
cap noisily wooldid
di "MECH_RC2=" _rc
di "MECH_END M6"

di "SENTINEL_DONE"
