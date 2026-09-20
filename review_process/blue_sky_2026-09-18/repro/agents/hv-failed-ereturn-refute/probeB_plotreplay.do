*! probeB_plotreplay.do -- SKEPTIC's own probes, beyond the seeded claims.
*!  B1: a plot-SAVE failure after a COMPLETED estimation -> is the replay of a
*!      complete, valid e() then refused with "the stored results predate that
*!      failed call"?  (rerunning a do-file with saveplots() and no replace)
*!  B2: the realistic vector for HV-EPOST-01: a loop over specifications with
*!      `cap noi wooldid ... ; est store m`k''.
*! Self-contained. Run:
*!   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probeB_plotreplay.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771113
set sortseed 771113
which wooldid

* clean slate for the saveplots stem (makes this do-file idempotent)
local fl : dir "." files "pbstem*"
foreach f of local fl {
  erase "`f'"
}

*---------------------------------------------------------------------
* data: 40 counties in 8 regions x 8 years; region-level adoption
*---------------------------------------------------------------------
clear
set obs 40
gen int unit   = _n
gen byte region = ceil(unit/5)
gen int onset = .
replace onset = 2015 if region <= 2
replace onset = 2017 if region > 2 & region <= 4
gen double ai = rnormal()
expand 8
bysort unit: gen int year = 2011 + _n
gen long regyear = region*10000 + year
gen double x1 = rnormal()
gen byte post = (year >= onset) & !missing(onset)
gen double y  = 2 + ai + 0.12*(year-2011) + 1.5*post + 0.4*x1 + rnormal()
xtset unit year

*=====================================================================
* B1
*=====================================================================
dis "=============== B1 FIRST RUN (creates pbstem_ES.gph / .dta) ==============="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(pbstem)
dis "B1_FIRST_RC=" _rc
dis "B1_FIRST_ECMD=" e(cmd)
local b1good = e(b)[1,1]
dis "B1_FIRST_B11=" %21x `b1good'
dis "B1_FIRST_KB=" colsof(e(b))
dis "B1_FIRST_LASTRC=[${WOOLDID_LASTRC}]"
cap confirm file "pbstem_ES.gph"
dis "B1_GPH_EXISTS_RC=" _rc
cap confirm file "pbstem_ES.dta"
dis "B1_DTA_EXISTS_RC=" _rc

dis "=============== B1 SECOND RUN (identical line, file now exists) ==============="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(pbstem)
dis "B1_SECOND_RC=" _rc
dis "B1_SECOND_ECMD=" e(cmd)
dis "B1_SECOND_KB=" colsof(e(b))
dis "B1_SECOND_B11=" %21x e(b)[1,1]
dis "B1_SECOND_N=" e(N)
dis "B1_SECOND_LASTRC=[${WOOLDID_LASTRC}]"
* is the posted result COMPLETE and IDENTICAL to a successful run's?
dis "B1_SAME_AS_GOOD=" (e(b)[1,1] == `b1good')
local cn : colfullnames e(b)
dis "B1_SECOND_COLS=`cn'"
cap mat list e(wdidmainresults)
dis "B1_MAINMAT_RC=" _rc

dis "=============== B1 REPLAY OF THAT COMPLETE RESULT ==============="
cap noi wooldid
dis "B1_REPLAY_RC=" _rc

dis "=============== B1 CONTROL: same run WITH replace ==============="
cap noi wooldid y unit year onset, att esprelength(2) espostlength(2) makeplots saveplots(pbstem) replace
dis "B1_THIRD_RC=" _rc
dis "B1_THIRD_LASTRC=[${WOOLDID_LASTRC}]"
cap noi wooldid
dis "B1_THIRD_REPLAY_RC=" _rc

*=====================================================================
* B2 -- the realistic harvest vector for HV-EPOST-01
*=====================================================================
dis "=============== B2 LOOP OVER SPECIFICATIONS ==============="
est clear
local k = 0
foreach sp in "att" "att fe(regyear)" "att controls(x1)" {
  local ++k
  cap noi wooldid y unit year onset, `sp'
  local rc`k' = _rc
  cap qui est store m`k'
  dis "B2_SPEC`k'_RC=`rc`k''  ECMD=" e(cmd) "  KB=" colsof(e(b)) "  B11=" %21x e(b)[1,1] "  N=" e(N)
}
dis "=============== B2 WHAT THE HARVEST CONTAINS ==============="
qui est dir
dis "B2_ESTDIR=[" r(names) "]"
forvalues j = 1/3 {
  qui est restore m`j'
  dis "B2_STORED`j'_CMD=" e(cmd) "  KB=" colsof(e(b)) "  B11=" %21x e(b)[1,1]
}
cap noi estimates table m1 m2 m3, keep(main_att) b(%12.6f)
dis "B2_ESTTAB_RC=" _rc

dis "=============== PROBE B COMPLETE ==============="
