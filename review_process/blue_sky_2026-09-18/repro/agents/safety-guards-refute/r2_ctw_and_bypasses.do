* R2 -- SAFETY-02 (contreatwithin working-value rewrite), SAFETY-04 (weight domain),
*       SAFETY-05 (esf pre/post degradation), plus the treated-cohort-0 wrong-number hunt.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*=====================================================================
* BLOCK A -- float 8-digit calendar codes: does contreatwithin kill the caution?
*=====================================================================
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int yr = 2000 + _n - 1
gen double ttreyr = .
replace ttreyr = 2003 + mod(unit,3) if tg==1
* 8-digit calendar codes stored as FLOAT (the sthlp's own named hazard)
gen float tfl    = 20000101 + (yr-2000)*10000
gen float ttrefl = cond(ttreyr<., 20000101 + (ttreyr-2000)*10000, .)
gen double dose  = runiform()*2
gen double ui = rnormal()
bysort unit (yr): replace ui = ui[1]
gen double y = 1 + 0.2*(yr-2000) + ui + rnormal() + 1.2*dose*(yr>=ttreyr & ttreyr<.)

dis "BLOCKA1_START binary, default safety"
cap noi wooldid y unit tfl ttrefl, att cluster(unit)
dis "A1_RC=" _rc
dis "BLOCKA1_END"

dis "BLOCKA2_START contreat WITHOUT contreatwithin, default safety"
cap noi wooldid y unit tfl ttrefl, att cluster(unit) contreat(dose) contreatcontrols(basic)
dis "A2_RC=" _rc
dis "BLOCKA2_END"

dis "BLOCKA3_START contreat WITH contreatwithin, no safety() typed"
cap noi wooldid y unit tfl ttrefl, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "A3_RC=" _rc
dis "BLOCKA3_END"

*=====================================================================
* BLOCK B -- a TREATED cohort coded 0 (0-indexed export). Does a number move?
*=====================================================================
clear
set obs 120
gen long unit = _n
gen byte tg = unit <= 90
* three treated cohorts + a never-treated group
gen long coh0 = cond(tg==1, mod(unit,3), 9)
gen long coh1 = coh0 + 1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2003 + mod(unit,3) if tg==1
gen double dose = runiform()*2
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t>=ttre & ttre<.)

dis "BLOCKB1_START binary default safety, 0-indexed cohort"
cap noi wooldid y coh0 t ttre, att cluster(unit)
dis "B1_RC=" _rc
dis "BLOCKB1_END"

dis "BLOCKB2_START binary safety(off), 0-indexed cohort"
cap noi wooldid y coh0 t ttre, att att_i cluster(unit) safety(off)
dis "B2_RC=" _rc
if _rc==0 {
  mat b2 = e(b)
  dis "B2_N=" e(N)
  local nm2 : colnames b2
  forvalues j = 1/`=colsof(b2)' {
    local thisn : word `j' of `nm2'
    dis "B2_COEF" `j' "=" %21x b2[1,`j'] " name=`thisn'"
  }
}
dis "BLOCKB2_END"

dis "BLOCKB3_START binary safety(off), 1-indexed cohort (control)"
cap noi wooldid y coh1 t ttre, att att_i cluster(unit) safety(off)
dis "B3_RC=" _rc
if _rc==0 {
  mat b3 = e(b)
  dis "B3_N=" e(N)
  local nm3 : colnames b3
  forvalues j = 1/`=colsof(b3)' {
    local thisn : word `j' of `nm3'
    dis "B3_COEF" `j' "=" %21x b3[1,`j'] " name=`thisn'"
  }
}
dis "BLOCKB3_END"

dis "BLOCKB4_START contreatwithin, 0-indexed cohort, no safety() typed"
cap noi wooldid y coh0 t ttre, att att_i cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "B4_RC=" _rc
if _rc==0 {
  mat b4 = e(b)
  local nm4 : colnames b4
  forvalues j = 1/`=colsof(b4)' {
    local thisn : word `j' of `nm4'
    dis "B4_COEF" `j' "=" %21x b4[1,`j'] " name=`thisn'"
  }
}
dis "BLOCKB4_END"

dis "BLOCKB5_START contreatwithin, 1-indexed cohort (control)"
cap noi wooldid y coh1 t ttre, att att_i cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "B5_RC=" _rc
if _rc==0 {
  mat b5 = e(b)
  local nm5 : colnames b5
  forvalues j = 1/`=colsof(b5)' {
    local thisn : word `j' of `nm5'
    dis "B5_COEF" `j' "=" %21x b5[1,`j'] " name=`thisn'"
  }
}
dis "BLOCKB5_END"

dis "BLOCKB6_START non-integer cohort under contreatwithin"
gen double cohfrac = coh1 + 0.5
cap noi wooldid y cohfrac t ttre, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "B6_RC=" _rc
dis "BLOCKB6_END"

*=====================================================================
* BLOCK C -- SAFETY-05: esfixedbaseperiod with onset at first observed period
*=====================================================================
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2000 if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t>=ttre & ttre<.)

dis "BLOCKC1_START esf default safety"
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(2) espostlength(2)
dis "C1_RC=" _rc
dis "BLOCKC1_END"

dis "BLOCKC2_START esf safety(off)"
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(2) espostlength(2) safety(off)
dis "C2_RC=" _rc
dis "BLOCKC2_END"

*=====================================================================
* BLOCK D -- SAFETY-04: negative aweight sitting only on rows markout discards
*=====================================================================
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2003 + mod(unit,3) if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t>=ttre & ttre<.)
gen double wgt = 0.5 + runiform()
* the outcome is missing on a handful of rows; the negative weight sits ONLY there
gen double ymiss = y
replace ymiss = . if inlist(_n, 7, 8, 9, 101, 102)
replace wgt = -0.5 if missing(ymiss)
qui count if missing(ymiss)
dis "D_NMISS=" r(N)
qui count if wgt < 0
dis "D_NNEGW=" r(N)
qui count if wgt < 0 & !missing(ymiss)
dis "D_NEGW_ON_USABLE=" r(N)

dis "BLOCKD1_START default safety, neg weight on missing-y rows only"
cap noi wooldid ymiss unit t ttre [aw=wgt], att cluster(unit)
dis "D1_RC=" _rc
dis "BLOCKD1_END"

dis "BLOCKD2_START safety(off), same call"
cap noi wooldid ymiss unit t ttre [aw=wgt], att cluster(unit) safety(off)
dis "D2_RC=" _rc
if _rc==0 {
  mat d2 = e(b)
  dis "D2_N=" e(N)
  dis "D2_ATT_HEX=" %21x d2[1,1]
}
dis "BLOCKD2_END"

dis "BLOCKD3_START default safety + if !missing(ymiss) (the trivial workaround)"
cap noi wooldid ymiss unit t ttre [aw=wgt] if !missing(ymiss), att cluster(unit)
dis "D3_RC=" _rc
if _rc==0 {
  mat d3 = e(b)
  dis "D3_N=" e(N)
  dis "D3_ATT_HEX=" %21x d3[1,1]
}
dis "BLOCKD3_END"

dis "R2_DONE"
