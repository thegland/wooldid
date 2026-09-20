*! probeD_histsave.do -- SKEPTIC's own finding, part 2.
*! `histogramcohorteffects' + `saveplots()' writes its .gph at ado:12390, which is
*! reached from ado:1944 -- i.e. LONG BEFORE wooldid's own `ereturn post' at ado:4051.
*! So merely RERUNNING A DO-FILE (file exists, no `replace') raises r(602) inside the
*! worker before anything is posted, and e() is left holding reghdfe's internal
*! saturated model.  A file-already-exists rerun is therefore an HV-EPOST-01 trigger.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probeD_histsave.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771113
set sortseed 771113
which wooldid

local fl : dir "." files "pestem*"
foreach f of local fl {
  erase "`f'"
}

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
gen double x1 = rnormal()
gen byte post = (year >= onset) & !missing(onset)
gen double y  = 2 + ai + 0.12*(year-2011) + 1.5*post + 0.4*x1 + rnormal()
xtset unit year

dis "=============== D0 a plain good run to occupy e() ==============="
cap noi wooldid y unit year onset, att
dis "D0_RC=" _rc
dis "D0_ECMD=" e(cmd) "  KB=" colsof(e(b))
local dgood = e(b)[1,1]
dis "D0_B11=" %21x `dgood'
est store dgood

dis "=============== D1 RUN 1 (hist + saveplots, writes pestem_hist_*) ==============="
cap noi wooldid y unit year onset, att_i histogramcohorteffects saveplots(pestem)
dis "D1_RC=" _rc
dis "D1_ECMD=" e(cmd) "  KB=" colsof(e(b))
local fl : dir "." files "pestem*"
dis `"D1_FILES=`fl'"'

dis "=============== D2 restore the good run, then RERUN the identical hist line ==============="
qui est restore dgood
dis "D2_PRE_ECMD=" e(cmd) "  B11=" %21x e(b)[1,1]
cap noi wooldid y unit year onset, att_i histogramcohorteffects saveplots(pestem)
dis "D2_RC=" _rc
dis "D2_ECMD=" e(cmd)
local kb = .
cap local kb = colsof(e(b))
dis "D2_KB=`kb'"
local b11 = .
cap local b11 = e(b)[1,1]
dis "D2_B11=" %21x `b11'
dis "D2_EN=" e(N)
local cn
cap local cn : colfullnames e(b)
dis "D2_COL1=" word("`cn'", 1)
dis "D2_TRUTH_B11=" %21x `dgood'
dis "D2_LASTRC=[${WOOLDID_LASTRC}]"
qui est dir
dis "D2_ESTDIR=[" r(names) "]"
qui frame
dis "D2_FRAME=" r(currentframe)
cap confirm variable _ppmlhdfe_d
dis "D2_PPMLD_RC=" _rc
dis "=== D2 replay ==="
cap noi wooldid
dis "D2_REPLAY_RC=" _rc

dis "=============== D3 CONTROL: same rerun WITH replace ==============="
qui est restore dgood
cap noi wooldid y unit year onset, att_i histogramcohorteffects saveplots(pestem) replace
dis "D3_RC=" _rc
dis "D3_ECMD=" e(cmd)
dis "D3_LASTRC=[${WOOLDID_LASTRC}]"

dis "=============== PROBE D COMPLETE ==============="
