*! repro_saveplots_rerun.do
*! SELF-CONTAINED repro for two findings a skeptic can check without any notes.
*!
*! NF-SKEP-01  saveplots() gets no pre-flight file-existence check, so simply
*!             RERUNNING A DO-FILE (the .gph from the first run is still there and
*!             `replace' was not given) raises r(602) and loses the whole wooldid
*!             call -- even though the estimation itself succeeded.  permsave() and
*!             bootsave() were given exactly this pre-flight check on 2026-09-05
*!             (ado:5866-5877, ado:6276-6298) because "a name clash cost the whole
*!             run" (ado:6276-6278).  saveplots() never got one, and under
*!             bootstrap() it pays that cost AFTER the entire replication loop.
*!             NOTE: that a plot failure ERRORS the call is deliberate -- ado:277-280
*!             says so.  The finding is the missing pre-flight on the one trigger
*!             that is trivially pre-flightable.  sthlp:696-697 documents a SOFT
*!             failure for the _cfx.dta write and is SILENT about these two.
*!
*! NF-SKEP-02  When that r(602) lands AFTER wooldid's ereturn post (the event-study
*!             .gph at ado:13325), e() holds the COMPLETE, VALID results of that very
*!             call -- the tables were already printed -- yet bare `wooldid' refuses
*!             to replay them saying "the stored results predate that failed call".
*!             The refusal is design (ado:277-280); the MESSAGE TEXT is false here.
*!
*! PART C also shows the same r(602) landing BEFORE the post (histogramcohorteffects
*! draws at ado:1944, .gph saved at ado:12390), which leaves reghdfe's internal
*! saturated model in e() and destroys the user's previous wooldid results.
*!
*! Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_saveplots_rerun.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

local fl : dir "." files "rrstem*"
foreach f of local fl {
  erase "`f'"
}

* ---- an ordinary 40-county x 8-year panel, policy adopted at the region level
clear
set obs 40
gen int unit    = _n
gen byte region = ceil(unit/5)
gen int onset = .
replace onset = 2015 if region <= 2
replace onset = 2017 if region > 2 & region <= 4
gen double ai = rnormal()
expand 8
bysort unit: gen int year = 2011 + _n
gen double x1 = rnormal()
gen byte post = (year >= onset) & !missing(onset)
gen double y  = 2 + ai + 0.12*(year-2011) + 1.5*post + 0.4*x1 + rnormal()
xtset unit year

*=========================================================================
* PART A -- a good run, kept as the "user's previous results"
*=========================================================================
wooldid y unit year onset, att
local truth = e(b)[1,1]
dis "A_TRUTH_main_att=" %21x `truth'
est store userprev

*=========================================================================
* PART B -- NF-SKEP-02: r(602) AFTER the post.  Run the SAME line twice.
*=========================================================================
dis "== B run 1 =="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(rrstem)
dis "B_RUN1_RC=" _rc
dis "== B run 2 (identical line; rrstem_ES_att.gph now exists) =="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(rrstem)
dis "B_RUN2_RC=" _rc
dis "B_RUN2_ECMD=" e(cmd)
dis "B_RUN2_KB=" colsof(e(b))
dis "B_RUN2_N=" e(N)
cap mat list e(wdidmainresults)
dis "B_RUN2_MAINMAT_PRESENT_RC=" _rc
dis "  ^ e() is the COMPLETE result of run 2: e(cmd)=wooldid, all 6 columns,"
dis "    e(wdidmainresults) present, and its table was printed above."
dis "== B replay of those complete results =="
cap noi wooldid
dis "B_REPLAY_RC=" _rc
dis "  ^ refused as 'predating' a call whose own complete results these are."

*=========================================================================
* PART C -- NF-SKEP-01 / HV-EPOST-01: r(602) BEFORE the post.
*=========================================================================
qui est restore userprev
dis "C_PRE_ECMD=" e(cmd) "  main_att=" %21x e(b)[1,1]
dis "== C run 1 =="
cap noi wooldid y unit year onset, att_i histogramcohorteffects saveplots(rrstem)
dis "C_RUN1_RC=" _rc
qui est restore userprev
dis "== C run 2 (identical line; rrstem_hist_att_i.gph now exists) =="
cap noi wooldid y unit year onset, att_i histogramcohorteffects saveplots(rrstem)
dis "C_RUN2_RC=" _rc
dis "C_RUN2_ECMD=" e(cmd)
dis "C_RUN2_KB=" colsof(e(b))
dis "C_RUN2_N=" e(N)
dis "C_RUN2_B11=" %21x e(b)[1,1]
local cn : colfullnames e(b)
dis "C_RUN2_COL1=" word("`cn'", 1)
dis "C_TRUTH_main_att=" %21x `truth'
dis "  ^ the user's previous wooldid results are gone; e() is reghdfe's internal"
dis "    saturated model, e(N) unchanged at 320, e(b)[1,1] an omitted base cell = 0."
cap noi wooldid
dis "C_REPLAY_RC=" _rc

*=========================================================================
* PART D -- control: the identical reruns WITH `replace' both succeed
*=========================================================================
qui est restore userprev
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(rrstem) replace
dis "D_RC_makeplots=" _rc
cap noi wooldid y unit year onset, att_i histogramcohorteffects saveplots(rrstem) replace
dis "D_RC_hist=" _rc
dis "D_LASTRC=[${WOOLDID_LASTRC}]"
dis "=============== REPRO COMPLETE ==============="
