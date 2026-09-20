* =============================================================================
* SELF-CONTAINED REPRO -- skeptic finding A
*
* CLAIM: with contreat() + contreatcontrols(untreatedZero), the `regress'-based
* estimation path -- customvce(anything), unconditionalse -- and the poisson
* path silently estimate the contreatcontrols(BASIC) model instead.  Every
* post-treatment att and ame row (main, and every ES horizon) comes back as a
* different number, here about 45% of the correct one, with NO caution, the
* SAME e(N), the SAME e(r2), the SAME e(df_r), and IDENTICAL pre-treatment rows.
*
* CAUSE (wooldid.ado:1234): the block that inserts the pooled baseline dose
* control fires on
*     strpos(contreatcontrols,"basic") | unconditionalse | customvce() | poisson
* The untreatedZero branch twenty lines later (ado:1313) rewrites `contreatinj'
* but never clears `contreatcontrolinserterFE'/`ABS', so the baseline dose slope
* that untreatedZero exists to omit is put back on those paths.  Under an
* exposure x post dose it is exactly collinear with the treated-cell dose terms,
* so the fit is unchanged (same r2/rss) but the ame/att decomposition is not.
*
* GROUND TRUTH: the design's true average marginal effect of the dose is 0.8 at
* every horizon; an independent hand-built saturated regression (ORACLE below)
* gives the cell-weighted answer the DEFAULT path reproduces to 10 decimals.
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do finding_A_untreatedzero_regresspath.do
* Read: the VERDICT_* lines at the end of the log.
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
dis "=== build under test ==="
which wooldid

* ------------------------------------------------------------------- DGP ---
* hospital-year administrative panel.  i = state (12), 60 hospitals per state,
* T = 2008..2019.  Medicaid-expansion-style staggered onset 2012/2014/2016,
* three never-expanding states.  dose = hospital exposure, switched on at that
* state's own expansion year -- the canonical exposure x post continuous dose,
* which is exactly the premise contreatcontrols(untreatedZero) documents.
clear
set obs 12
gen byte state = _n
gen int expyr = .
replace expyr = 2012 if inlist(state,1,2,3,4)
replace expyr = 2014 if inlist(state,5,6,7)
replace expyr = 2016 if inlist(state,8,9)
expand 60
bys state: gen int hosp = _n
gen long hid = state*1000 + hosp
gen double expo = 0.2 + 1.6*runiform()
gen double hfe  = rnormal(0,1)
expand 12
bys hid: gen int year = 2007 + _n
gen double dose = expo*(expyr<. & year>=expyr)
gen double y    = 4 + hfe + 0.03*(year-2007) + 0.8*dose + rnormal(0,0.05)
gen byte trt    = (expyr<. & year>=expyr)
gen double ypois = exp(0.2 + 0.05*hfe + 0.01*(year-2007) + 0.25*dose)
compress
save "fa_panel.dta", replace
dis "N_rows=" _N

* ========================== ORACLE (no wooldid) =============================
* the saturated untreatedZero regression written out by hand:
*   y = a_state + d_year + sum_treatedcells [ g_c*1{cell} + b_c*1{cell}*dose ]
* ame = obs-weighted mean of b_c ; att = obs-weighted mean of g_c
egen long cellg = group(state year) if trt == 1
replace cellg = 0 if missing(cellg)
qui levelsof cellg if cellg > 0, local(cells)
qui reghdfe y ib0.cellg ib0.cellg#c.dose, absorb(state year)
local sb 0
local sg 0
local sw 0
foreach c of local cells {
  qui count if cellg == `c'
  local w = r(N)
  local sb = `sb' + `w'*_b[`c'.cellg#c.dose]
  local sg = `sg' + `w'*_b[`c'.cellg]
  local sw = `sw' + `w'
}
scalar ORACLE_AME = `sb'/`sw'
scalar ORACLE_ATT = `sg'/`sw'
dis "ORACLE_AME=" %16.10f ORACLE_AME
dis "ORACLE_ATT=" %16.10f ORACLE_ATT

* ============================== wooldid runs ================================
cap prog drop grab
prog def grab
  syntax , tag(string)
  matrix BB = e(b)
  local bn : colfullnames BB
  local j 0
  foreach nm of local bn {
    local ++j
    if "`nm'" == "main_att"      scalar A_`tag' = BB[1,`j']
    if "`nm'" == "main_ame"      scalar M_`tag' = BB[1,`j']
    if "`nm'" == "ES_att_post0"  scalar E0_`tag' = BB[1,`j']
    if "`nm'" == "ES_ame_post0"  scalar F0_`tag' = BB[1,`j']
  }
  scalar N_`tag' = e(N)
  scalar R_`tag' = e(r2)
  scalar D_`tag' = e(df_r)
end

use "fa_panel.dta", clear
dis "###### 1. DEFAULT path, contreatcontrols(untreatedZero) -- the correct answer ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero)
dis "RC=" _rc
grab , tag(UZ_DEF)

use "fa_panel.dta", clear
dis "###### 2. customvce(ols), contreatcontrols(untreatedZero) ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) customvce(ols)
dis "RC=" _rc
grab , tag(UZ_CV)

use "fa_panel.dta", clear
dis "###### 3. unconditionalse + cluster(hid), contreatcontrols(untreatedZero) ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) unconditionalse cluster(hid)
dis "RC=" _rc
grab , tag(UZ_UN)

use "fa_panel.dta", clear
dis "###### 4. DEFAULT path, contreatcontrols(basic) -- the model the others actually fit ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(basic)
dis "RC=" _rc
grab , tag(BA_DEF)

use "fa_panel.dta", clear
dis "###### 5. robust (reghdfe path) -- control: must equal run 1 ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) robust
dis "RC=" _rc
grab , tag(UZ_RB)

use "fa_panel.dta", clear
dis "###### 6. BINARY treatment, default vs customvce(ols) -- control: must agree ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3)
grab , tag(BIN_DEF)
use "fa_panel.dta", clear
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) customvce(ols)
grab , tag(BIN_CV)

* ---- poisson arm: ado:1234 lists poisson in the same condition --------------
use "fa_panel.dta", clear
dis "###### 7. poisson + untreatedZero ######"
cap noisily wooldid ypois state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) poisson
dis "RC7=" _rc
if _rc == 0 grab , tag(PO_UZ)
use "fa_panel.dta", clear
dis "###### 8. poisson + basic ######"
cap noisily wooldid ypois state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(basic) poisson
dis "RC8=" _rc
if _rc == 0 grab , tag(PO_BA)

* ================================ VERDICT ==================================
dis ""
dis "==================== VERDICT ===================="
dis "VERDICT_oracle_ame      = " %16.10f ORACLE_AME
dis "VERDICT_oracle_att      = " %16.10f ORACLE_ATT
dis "VERDICT_default_uz_ame  = " %16.10f M_UZ_DEF  "   (matches oracle)"
dis "VERDICT_default_uz_att  = " %16.10f A_UZ_DEF
dis "VERDICT_robust_uz_ame   = " %16.10f M_UZ_RB   "   (reghdfe path: same)"
dis "VERDICT_customvce_uz_ame= " %16.10f M_UZ_CV
dis "VERDICT_customvce_uz_att= " %16.10f A_UZ_CV
dis "VERDICT_uncond_uz_ame   = " %16.10f M_UZ_UN
dis "VERDICT_default_BASIC_ame=" %16.10f M_BA_DEF
dis "VERDICT_default_BASIC_att=" %16.10f A_BA_DEF
dis ""
dis "VERDICT_cv_equals_BASIC_ame = " (reldif(M_UZ_CV, M_BA_DEF) < 1e-9)
dis "VERDICT_cv_equals_BASIC_att = " (reldif(A_UZ_CV, A_BA_DEF) < 1e-9)
dis "VERDICT_un_equals_BASIC_ame = " (reldif(M_UZ_UN, M_BA_DEF) < 1e-9)
dis "VERDICT_cv_differs_from_DEF = " (reldif(M_UZ_CV, M_UZ_DEF) > 1e-6)
dis "VERDICT_ratio_ame_cv_over_def = " %8.5f (M_UZ_CV/M_UZ_DEF)
dis "VERDICT_ratio_att_cv_over_def = " %8.5f (A_UZ_CV/A_UZ_DEF)
dis "VERDICT_same_N   = " (N_UZ_CV == N_UZ_DEF) "  (" N_UZ_DEF ")"
dis "VERDICT_same_r2  = " (reldif(R_UZ_CV, R_UZ_DEF) < 1e-12) "  (" %12.9f R_UZ_DEF ")"
dis "VERDICT_same_dfr = " (D_UZ_CV == D_UZ_DEF)
dis "VERDICT_binary_paths_agree = " (reldif(A_BIN_CV, A_BIN_DEF) < 1e-12)
cap dis "VERDICT_poisson_uz_equals_basic = " (reldif(M_PO_UZ, M_PO_BA) < 1e-9) "   uz=" %14.10f M_PO_UZ "  basic=" %14.10f M_PO_BA
dis "================================================"
dis "FA_DONE=1"
