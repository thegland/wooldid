clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
include "estate.doi"

mata:
mata clear
real scalar userfn()
{
  return(42)
}
end

* ------------- the GOOD run, re-established before every failure site -------------
cap prog drop goodrun
prog def goodrun
  use D1.dta, clear
  cap qui wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2)
  if _rc dis "GOODRUN_FAILED rc=" _rc
  cap est drop USERM1
  qui est store USERM1
  cap matrix drop USERMAT
  matrix USERMAT = J(2,2,7)
  cap scalar drop userscal
  scalar userscal = 42
end

goodrun
estate GOOD 0

* ================= FAILURE SITE ENUMERATION =================
* each: re-establish good state, then run the failing call, then dump

#delimit ;
local nsite = 0 ;
#delimit cr

* --- F01 parse: unrecognised option ---
goodrun
cap noi wooldid y id t ttre, attt
estate F01_badopt `=_rc'

* --- F02 too few variables ---
goodrun
cap noi wooldid y id t, att
estate F02_fewvars `=_rc'

* --- F03 nonexistent variable ---
goodrun
cap noi wooldid y id t nosuchvar, att
estate F03_novar `=_rc'

* --- F04 bad if expression ---
goodrun
cap noi wooldid y id t ttre if nosuchvar==1, att
estate F04_badif `=_rc'

* --- F05 no never-treated ---
goodrun
qui replace ttre = 5 if missing(ttre)
cap noi wooldid y id t ttre, att
estate F05_nonevertr `=_rc'

* --- F06 negative aweights ---
goodrun
qui gen double negw = awt - 0.75
cap noi wooldid y id t ttre [aw=negw], att
estate F06_negw `=_rc'

* --- F07 fweights (documented refusal) ---
goodrun
qui gen int fw1 = 1 + mod(id,3)
cap noi wooldid y id t ttre [fw=fw1], att
estate F07_fweight `=_rc'

* --- F08 contreat with unit-level i (collinear dose, r322) ---
goodrun
qui gen double dose = runiform()*(t>=ttre & ttre<.)
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(untreatedZero)
estate F08_contreatcollin `=_rc'

* --- F09 contreat without contreatcontrols ---
goodrun
qui gen double dose = runiform()
cap noi wooldid y id t ttre, att contreat(dose)
estate F09_contreatnoctl `=_rc'

* --- F10 name-length overflow: 8-digit subgroup values, all four estimands ---
goodrun
qui gen long sgbig = 10000000 + mod(id,2)
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod subgroups(sgbig)
estate F10_nameoverflow `=_rc'

* --- F11 straddling subgroup ---
goodrun
qui gen byte sgstraddle = mod(_n,2)
cap noi wooldid y id t ttre, att subgroups(sgstraddle)
estate F11_sgstraddle `=_rc'

* --- F12 permreps below the floor ---
goodrun
cap noi wooldid y id t ttre, att permute(shuffle) permreps(10)
estate F12_permfloor `=_rc'

* --- F13 stray companion option ---
goodrun
cap noi wooldid y id t ttre, att reps(40)
estate F13_straycomp `=_rc'

* --- F14 poisson on a negative outcome ---
goodrun
cap noi wooldid y id t ttre, att poisson
estate F14_poisneg `=_rc'

* --- F15 pretrendpower with no pre periods ---
goodrun
cap noi wooldid y id t ttre, att pretrendpower(0.5)
estate F15_ptpnopre `=_rc'

* --- F16 ccc_absorb with unit-level i ---
goodrun
cap noi wooldid y id t ttre, att ccc_absorb(id) esprelength(2) espostlength(2)
estate F16_cccabsorb `=_rc'

* --- F17 subgroups() on a constant variable ---
goodrun
qui gen byte sgconst = 1
cap noi wooldid y id t ttre, att subgroups(sgconst)
estate F17_sgconst `=_rc'

* --- F18 cluster variable with missing values ---
goodrun
qui gen long clum = clu
qui replace clum = . if mod(id,7)==0
cap noi wooldid y id t ttre, att cluster(clum)
estate F18_clumiss `=_rc'

* --- F19 a run whose outcome is all missing inside the if ---
goodrun
cap noi wooldid y id t ttre if t>100, att
estate F19_emptyif `=_rc'

* --- F20 bootstrap whose replications mostly fail (tiny cluster count) ---
goodrun
qui gen byte clu1 = 1
cap noi wooldid y id t ttre, att cluster(clu1) bootstrap(standard) reps(5) rseed(7)
estate F20_bootfail `=_rc'

* ================= WHAT THE CONSUMERS DO AFTER A FAILED RUN =================
dis "======== CONSUMER BEHAVIOUR AFTER F08 (r322 post-regression failure) ========"
goodrun
dis "PRE-FAIL e(cmd)=" "`e(cmd)'" " K=" colsof(e(b))
qui gen double dose = runiform()*(t>=ttre & ttre<.)
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(untreatedZero)
dis "FAILRC=" _rc
estate AFTERFAIL `=_rc'
cap noi ereturn display
dis "EREDISP_RC=" _rc
cap noi predict double yhat
dis "PREDICT_RC=" _rc
cap sum yhat
dis "PREDICT_N=" r(N)
cap drop yhat
cap noi predict double xbh, xb
dis "PREDICTXB_RC=" _rc
cap sum xbh
dis "PREDICTXB_N=" r(N)
cap drop xbh
cap noi estimates store AFTERF
dis "ESTSTORE_RC=" _rc
cap noi test _b[_cons]==0
dis "TEST_RC=" _rc
cap noi estat vce
dis "ESTATVCE_RC=" _rc
cap noi wooldid
dis "BAREREPLAY_RC=" _rc

* ================= REPLAY-STALENESS GUARD =================
dis "======== REPLAY GUARD ========"
goodrun
cap noi wooldid
dis "REPLAY_AFTER_GOOD_RC=" _rc
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "PARSEFAIL_RC=" _rc
dis "e(cmd) after parse fail = `e(cmd)'"
cap noi wooldid
dis "REPLAY_AFTER_PARSEFAIL_RC=" _rc
macro drop WOOLDID_LASTRC
cap noi wooldid
dis "REPLAY_AFTER_MACRODROP_RC=" _rc

* does clear all wipe the marker?
goodrun
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "LASTRC_SET=${WOOLDID_LASTRC}"
clear
dis "LASTRC_AFTER_CLEAR=${WOOLDID_LASTRC}"
discard
dis "LASTRC_AFTER_DISCARD=${WOOLDID_LASTRC}"

* ================= est save / restore ROUND TRIP =================
dis "======== EST SAVE/RESTORE ROUND TRIP ========"
goodrun
local sc0 : e(scalars)
local mc0 : e(macros)
local mt0 : e(matrices)
dis "PRE_NSCA=" `:word count `sc0''  " NMAC=" `:word count `mc0'' " NMAT=" `:word count `mt0''
qui est store RT1
cap erase rt.ster
cap noi estimates save rt, replace
dis "ESTSAVE_RC=" _rc
clear
cap noi estimates use rt
dis "ESTUSE_RC=" _rc
local sc1 : e(scalars)
local mc1 : e(macros)
local mt1 : e(matrices)
dis "POST_NSCA=" `:word count `sc1''  " NMAC=" `:word count `mc1'' " NMAT=" `:word count `mt1''
local lostS : list sc0 - sc1
local lostM : list mc0 - mc1
local lostX : list mt0 - mt1
dis "LOST_SCALARS: `lostS'"
dis "LOST_MACROS: `lostM'"
dis "LOST_MATRICES: `lostX'"
dis "RESTORED_CMD=`e(cmd)'"
cap noi wooldid
dis "REPLAY_AFTER_ESTUSE_RC=" _rc

* ================= store -> fail -> restore -> replay =================
dis "======== STORE / FAIL / RESTORE / REPLAY ========"
goodrun
qui est drop USERM1
qui est store GOODM
dis "STORED e(cmd)=`e(cmd)' K=" colsof(e(b))
cap noi wooldid y id t ttre, att subgroups(nosuchvar)
dis "FAIL_RC=" _rc " lastrc=${WOOLDID_LASTRC}"
dis "AFTER_FAIL e(cmd)=`e(cmd)'"
cap noi estimates restore GOODM
dis "RESTORE_RC=" _rc " e(cmd)=`e(cmd)' K=" colsof(e(b))
dis "LASTRC_STILL=${WOOLDID_LASTRC}"
cap noi wooldid
dis "REPLAY_AFTER_RESTORE_RC=" _rc

* the same with a post-regression failure (F1 shape)
goodrun
qui est drop USERM1
qui est store GOODM2
qui gen double dose2 = runiform()*(t>=ttre & ttre<.)
cap noi wooldid y id t ttre, att contreat(dose2) contreatcontrols(untreatedZero)
dis "F1FAIL_RC=" _rc " e(cmd)=`e(cmd)'"
cap noi estimates restore GOODM2
dis "RESTORE2_RC=" _rc " e(cmd)=`e(cmd)' K=" colsof(e(b))
cap noi wooldid
dis "REPLAY_AFTER_RESTORE2_RC=" _rc

dis "STALE_DONE"
