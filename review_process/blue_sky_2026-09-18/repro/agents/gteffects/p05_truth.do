* p05_truth.do -- EXACT-TRUTH oracle for the cell VALUES (not just the identity).
* With ysd=0 the DGP lies exactly in the span of (unit FE, period FE, cell dummies), so
* every cohort-by-relative-period coefficient must equal the cell's true average effect
* to machine precision, and every pre-treatment cell must be exactly 0.
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

* truthcheck: compare each table row's estimate against the weighted mean of the true
* effect te over the rows of that cell, and check the onset/caltime/relyr/nobs labels.
cap prog drop truthcheck
prog def truthcheck
  syntax , TAG(string) IV(varname) [WV(string)]
  cap confirm matrix e(gteffects)
  if _rc != 0 {
    dis as result "TR|`tag'|NOTABLE"
    exit
  }
  tempname GT
  matrix `GT' = e(gteffects)
  tempvar es
  qui gen byte `es' = e(sample)
  local nr = rowsof(`GT')
  local worst 0
  local worstrow 0
  local nbadn 0
  local nbadct 0
  forval r = 1/`nr' {
    local gi = `GT'[`r',1]
    local gon = `GT'[`r',2]
    local gct = `GT'[`r',3]
    local gry = `GT'[`r',4]
    local ges = `GT'[`r',5]
    local gnb = `GT'[`r',11]
    qui count if `es' & `iv' == `gi' & year == `gct'
    if r(N) != `gnb' local ++nbadn
    if "`wv'" == "" qui sum te if `es' & `iv' == `gi' & year == `gct', meanonly
    else {
      tempvar num
      qui sum te if `es' & `iv' == `gi' & year == `gct', meanonly
    }
    local truth = r(mean)
    if "`wv'" != "" {
      qui sum `wv' if `es' & `iv' == `gi' & year == `gct', meanonly
      local sw = r(sum)
      tempvar prod
      qui gen double `prod' = te*`wv' if `es' & `iv' == `gi' & year == `gct'
      qui sum `prod', meanonly
      local truth = r(sum)/`sw'
      drop `prod'
    }
    local d = abs(`ges' - `truth')
    if `d' > `worst' {
      local worst `d'
      local worstrow `r'
    }
    if abs(`gct' - `gon' - `gry') > 1e-9 local ++nbadct
  }
  dis as result "TR|`tag'|ROWS|`nr'|MAXABSDEV|" %10.3e `worst' "|worstrow|`worstrow'|NOBSMISMATCH|`nbadn'|CALTIMEBAD|`nbadct'"
end

*---- T1: i = unit, balanced, zero noise -------------------------------------
mkpanel , n(30) t(9) onsets("2004 2006") ysd(0) gb1(0) seed(5150)
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects esprelength(3) espostlength(3)
dis "T1_RC=" _rc
gtcheck , tag(T1_truth_unit) pre es
truthcheck , tag(T1_truth_unit) iv(unit)

*---- T2: i = cohort, no unit FE, balanced, zero noise ------------------------
mkpanel , n(60) t(9) onsets("2004 2006") ysd(0) aisd(0) gb1(0) seed(5151)
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(3) espostlength(3)
dis "T2_RC=" _rc
gtcheck , tag(T2_truth_coh) pre es
truthcheck , tag(T2_truth_coh) iv(cohp)

*---- T3: i = cohort, weights, unbalanced, zero noise -------------------------
mkpanel , n(80) t(9) onsets("2004 2006") ysd(0) aisd(0) gb1(0) gapprob(0.1) seed(5152)
bys unit (year): gen double wt = 0.4 + 2*runiform() if _n==1
bys unit (year): replace wt = wt[1]
ereturn clear
cap noi wooldid y cohp year ttre [aw=wt], `EST' gteffects esprelength(3) espostlength(3)
dis "T3_RC=" _rc
gtcheck , tag(T3_truth_aw) pre es
truthcheck , tag(T3_truth_aw) iv(cohp) wv(wt)

*---- T4: out-of-window post periods (short window, long realized horizon) ----
mkpanel , n(60) t(12) onsets("2003") ysd(0) aisd(0) gb1(0) seed(5153)
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(1) espostlength(1)
dis "T4_RC=" _rc
gtcheck , tag(T4_shortwin) pre es
truthcheck , tag(T4_shortwin) iv(cohp)
cap noi matrix list e(gteffects), format(%12.8g)

dis "DONE_P05"
