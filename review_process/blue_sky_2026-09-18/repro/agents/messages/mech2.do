clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
use tiny.dta, clear

*** M2 retry: contreatwithin with a compatible contreatcontrols type
di "MECH_BEGIN M2d_ctw_basic_negi"
cap noisily wooldid y idneg t ttre, att fastnose contreat(dose) contreatcontrols(basic) contreatwithin
di "MECH_RC=" _rc
di "MECH_END M2d_ctw_basic_negi"

di "MECH_BEGIN M2e_noctw_basic_negi"
cap noisily wooldid y idneg t ttre, att fastnose contreat(dose) contreatcontrols(basic)
di "MECH_RC=" _rc
di "MECH_END M2e_noctw_basic_negi"

*** M2f: float-above-2^24 caution, with and without contreatwithin
gen float bigid = 20000000 + id
di "MECH_BEGIN M2f_float_noctw"
cap noisily wooldid y cid t ttre, att fastnose contreat(dose) contreatcontrols(basic)
di "MECH_RC=" _rc
di "MECH_END M2f_float_noctw"

*** M2g: negative weights (the safetyinput-gated sibling) under contreatwithin
di "MECH_BEGIN M2g_negw_ctw"
cap noisily wooldid y cid t ttre [aw=negw], att fastnose contreat(dose) contreatcontrols(basic) contreatwithin
di "MECH_RC=" _rc
di "MECH_END M2g_negw_ctw"

*** M7: subgroup value-length caution
gen long sgbig = 123456789012345 + grp
di "MECH_BEGIN M7_sgnames"
cap noisily wooldid y cid t ttre, att fastnose subgroups(sgbig)
di "MECH_RC=" _rc
di "MECH_END M7_sgnames"

*** M8: bootctl low-reps + few clusters caution stack
di "MECH_BEGIN M8_boot"
cap noisily wooldid y cid t ttre, att cluster(cid) bootstrap(standard) reps(20) dots(0)
di "MECH_RC=" _rc
di "MECH_END M8_boot"

*** M9: caution before abort in bootctl (bayes cell-weighted under safety(off) then 1 cluster)
gen byte onecl = 1
di "MECH_BEGIN M9_bootabort"
cap noisily wooldid y cid t ttre, att_i cluster(onecl) bootstrap(bayes) reps(20) safety(off) dots(0)
di "MECH_RC=" _rc
di "MECH_END M9_bootabort"

*** M10: bare wooldid replay after success
di "MECH_BEGIN M10_replay"
cap noisily wooldid y cid t ttre, att fastnose
di "MECH_RC=" _rc
cap noisily wooldid
di "MECH_RC2=" _rc
di "MECH_END M10_replay"

di "SENTINEL_DONE"
