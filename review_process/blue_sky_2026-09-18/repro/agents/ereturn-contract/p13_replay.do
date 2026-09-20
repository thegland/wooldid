* ===== the replay guard, est save/restore round-trip, and store->fail->restore =====
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
cap prog drop goodrun
prog def goodrun
  use D1.dta, clear
  cap qui wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2)
  if _rc dis "GOODRUN_FAILED rc=" _rc
end

dis "======== 1. REPLAY GUARD ========"
goodrun
cap noi wooldid
dis "R1_REPLAY_AFTER_GOOD_RC=" _rc
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "R2_PARSEFAIL_RC=" _rc " e(cmd)=`e(cmd)' LASTRC=${WOOLDID_LASTRC}"
cap noi wooldid
dis "R3_REPLAY_AFTER_PARSEFAIL_RC=" _rc
macro drop WOOLDID_LASTRC
cap noi wooldid
dis "R4_REPLAY_AFTER_MACRODROP_RC=" _rc
dis "R4_SHOWS_CMDLINE=`e(cmdline)'"

dis "======== 2. does the marker survive clear / discard? ========"
goodrun
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "M1_LASTRC_SET=${WOOLDID_LASTRC}"
clear
dis "M2_LASTRC_AFTER_CLEAR=${WOOLDID_LASTRC}"
discard
dis "M3_LASTRC_AFTER_DISCARD=${WOOLDID_LASTRC}"
macro drop _all
dis "M4_LASTRC_AFTER_MACRODROP_ALL=${WOOLDID_LASTRC}"
cap noi wooldid
dis "M5_REPLAY_AFTER_MACRODROPALL_RC=" _rc

dis "======== 3. store -> fail -> restore -> replay ========"
goodrun
qui est store GOODM
dis "S1_STORED cmd=`e(cmd)' K=" colsof(e(b))
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "S2_FAIL_RC=" _rc " LASTRC=${WOOLDID_LASTRC}"
cap noi estimates restore GOODM
dis "S3_RESTORE_RC=" _rc " e(cmd)=`e(cmd)' K=" colsof(e(b))
dis "S4_LASTRC_STILL=${WOOLDID_LASTRC}"
cap noi wooldid
dis "S5_REPLAY_AFTER_RESTORE_RC=" _rc

dis "======== 4. a post-regression failure, then restore ========"
goodrun
qui est store GOODM2
qui gen double dose2 = runiform()*(t>=ttre & ttre<.)
cap noi wooldid y id t ttre, att contreat(dose2) contreatcontrols(untreatedZero)
dis "P1_FAIL_RC=" _rc " e(cmd)=`e(cmd)'"
cap noi estimates restore GOODM2
dis "P2_RESTORE_RC=" _rc " e(cmd)=`e(cmd)' K=" colsof(e(b))
cap noi wooldid
dis "P3_REPLAY_AFTER_RESTORE_RC=" _rc

dis "======== 5. est save / est use round trip (hidden entries) ========"
goodrun
local sc0 : e(scalars)
local mc0 : e(macros)
local mt0 : e(matrices)
dis "T0 NSCA=" `:word count `sc0'' " NMAC=" `:word count `mc0'' " NMAT=" `:word count `mt0''
cap erase rt.ster
cap noi estimates save rt, replace
dis "T1_ESTSAVE_RC=" _rc
clear
cap noi estimates use rt
dis "T2_ESTUSE_RC=" _rc
local sc1 : e(scalars)
local mc1 : e(macros)
local mt1 : e(matrices)
dis "T3 NSCA=" `:word count `sc1'' " NMAC=" `:word count `mc1'' " NMAT=" `:word count `mt1''
local lostS : list sc0 - sc1
local lostM : list mc0 - mc1
local lostX : list mt0 - mt1
dis "T4_LOST_SCALARS: `lostS'"
dis "T5_LOST_MACROS: `lostM'"
dis "T6_LOST_MATRICES: `lostX'"
dis "T7_RESTORED_CMD=`e(cmd)' K=" colsof(e(b))
cap count if e(sample)
dis "T8_ESAMPLE_AFTER_USE=" r(N) " rc=" _rc
cap noi wooldid
dis "T9_REPLAY_AFTER_ESTUSE_RC=" _rc

dis "REPLAY_DONE"
