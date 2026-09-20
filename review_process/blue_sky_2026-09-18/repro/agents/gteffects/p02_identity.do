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
local ESW esprelength(3) espostlength(3)

*========================= C1 balanced / i=unit ==========================
mkpanel , n(60) t(10) onsets("2004 2006 2008") ysd(0.3)
tempfile base
qui save `base'
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects `ESW'
dis "C1_RC=" _rc
gtcheck , tag(C1_bal_unit) pre es
gtdefcheck , tag(C1_bal_unit)

*========================= C2 balanced / i=cohort ========================
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW'
dis "C2_RC=" _rc
gtcheck , tag(C2_bal_coh) pre es
gtdefcheck , tag(C2_bal_coh)

*========================= C3 interior gaps ==============================
mkpanel , n(60) t(10) onsets("2004 2006 2008") ysd(0.3) gapprob(0.12) seed(31)
tempfile gapd
qui save `gapd'
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects `ESW'
dis "C3_RC=" _rc
gtcheck , tag(C3_gaps_unit) pre es
gtdefcheck , tag(C3_gaps_unit)
use `gapd', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW'
dis "C3b_RC=" _rc
gtcheck , tag(C3b_gaps_coh) pre es
gtdefcheck , tag(C3b_gaps_coh)

*========================= C4 entry / exit ===============================
mkpanel , n(80) t(10) onsets("2004 2006 2008") ysd(0.3) entryprob(0.3) exitprob(0.3) seed(77)
tempfile eed
qui save `eed'
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW'
dis "C4_RC=" _rc
gtcheck , tag(C4_entryexit_coh) pre es
gtdefcheck , tag(C4_entryexit_coh)

*========================= C5 aweights ===================================
use `base', clear
bys unit (year): gen double wt = 0.2 + 3*runiform() if _n == 1
bys unit (year): replace wt = wt[1]
tempfile wd
qui save `wd'
ereturn clear
cap noi wooldid y unit year ttre [aw=wt], `EST' gteffects `ESW'
dis "C5_RC=" _rc
gtcheck , tag(C5_aw_unit) pre es
gtdefcheck , tag(C5_aw_unit)
use `wd', clear
ereturn clear
cap noi wooldid y cohp year ttre [aw=wt], `EST' gteffects `ESW'
dis "C5b_RC=" _rc
gtcheck , tag(C5b_aw_coh) pre es
gtdefcheck , tag(C5b_aw_coh)

*========================= C6 unbalanced + aw + multi-row cells ==========
mkpanel , n(60) t(9) onsets("2004 2006") ysd(0.3) gapprob(0.08) entryprob(0.25) percell(3) seed(1234)
gen double wt = 0.3 + 2*runiform()
tempfile ub
qui save `ub'
ereturn clear
cap noi wooldid y cohp year ttre [aw=wt], `EST' gteffects `ESW'
dis "C6_RC=" _rc
gtcheck , tag(C6_ub_aw_rcs) pre es
gtdefcheck , tag(C6_ub_aw_rcs)

*========================= C7 collapse oracle ============================
use `ub', clear
ereturn clear
cap noi wooldid y cohp year ttre [aw=wt], `EST' gteffects `ESW' nofastcollapse
dis "C7_RC=" _rc
gtcheck , tag(C7_nofastcollapse) pre es

*========================= C8 timetrends =================================
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW' timetrends
dis "C8_RC=" _rc
gtcheck , tag(C8_timetrends) pre es

*========================= C9 esfixedbaseperiod ==========================
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW' esfixedbaseperiod
dis "C9_RC=" _rc
gtcheck , tag(C9_esfixedbase) pre es

*========================= C10 anticipation ==============================
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW' anticipation(1)
dis "C10_RC=" _rc
gtcheck , tag(C10_antic1) pre es

*========================= C11 level(90) =================================
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW' level(90)
dis "C11_RC=" _rc
gtcheck , tag(C11_level90) pre es

*========================= C12 biennial grid =============================
mkpanel , n(60) t(9) onsets("2004 2008 2012") ysd(0.3) tstep(2) t0(2000) seed(999)
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(6) espostlength(6)
dis "C12_RC=" _rc
gtcheck , tag(C12_biennial) pre es

dis "DONE_P02"
