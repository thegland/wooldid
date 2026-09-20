* p13_omit.do -- FORCE an omitted (o.) treated cohort-by-period coefficient and see
* (a) whether default safety catches it, (b) what e(gteffects) does under safety(off),
* (c) by how much the documented weighted-mean identity is then wrong.
* M1/M2 use a specification a real user writes: state-by-year fixed effects, with the
* treated cohort living entirely inside one state.  M3 is the deliberate control case.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime

clear
set obs 80
gen long unit = _n
gen byte st = 1 + mod(_n, 4)
gen byte cohp = cond(st == 4, 2, 1)
gen double ttre = cond(cohp == 2, 2004, .)
gen double ai = rnormal()
qui expand 9
bys unit: gen int year = 1999 + _n
tempvar ts
bys year: gen double `ts' = rnormal() if _n==1
bys year: replace `ts' = `ts'[1]
gen double te = 0
qui replace te = 1 + 0.25*(year-ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`ts' + te + rnormal()*0.3
* z is an exact indicator of ONE treated cohort-by-period cell
gen byte z = (cohp == 2 & year == 2006)
tempfile md
qui save `md'

dis "--- M1: state-by-year FE, default safety ---"
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) fe(i.st#i.year)
dis "M1_RC=" _rc

dis "--- M2: state-by-year FE, safety(off) ---"
use `md', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) fe(i.st#i.year) safety(off)
dis "M2_RC=" _rc
gtcheck , tag(M2_styr_off) pre es
cap noi matrix list e(gteffects), format(%12.6g)

dis "--- M3: a control that is exactly one treated cell, default safety ---"
use `md', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) controls(z)
dis "M3_RC=" _rc

dis "--- M4: same, safety(off): does a cell vanish and the identity break? ---"
use `md', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) controls(z) safety(off) verbose
dis "M4_RC=" _rc
gtcheck , tag(M4_zctrl_off) pre es
cap noi matrix list e(gteffects), format(%12.6g)
qui est restore ___wooldidfullmodel
local bn : colnames e(b)
local nom 0
foreach v of local bn {
  if strpos("`v'", "o.") != 0 & strpos("`v'", "0b.") == 0 {
    local ++nom
    dis "M4_OMIT: `v'"
  }
}
dis "M4_OMITTED_TERMS=" `nom'

dis "--- M5: baseline (no z), for the reference numbers ---"
use `md', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2)
dis "M5_RC=" _rc
gtcheck , tag(M5_base) pre es

dis "--- M6: permute() with i = unit and cluster(unit) ---"
mkpanel , n(60) t(9) onsets("2004 2006") ysd(0.4) seed(808)
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects esprelength(2) espostlength(2) permute(shuffle) permreps(20) cluster(unit)
dis "M6_RC=" _rc
dis "M6_perm_gt_blank=`e(perm_gt_blank)'"
gtcheck , tag(M6_permute) pre es
cap noi matrix list e(gteffects), format(%12.6g)

dis "DONE_P13"
