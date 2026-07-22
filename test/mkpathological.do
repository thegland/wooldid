* Build the pathological test datasets referenced by review_process/assessment_grid.md (D12).
* Run from test/:  cd test && "<stata>" -e do mkpathological.do
* Deterministic; regenerating overwrites data/dta_path_*.dta in place.
*
* Shared schema (all files, so any config can run on any shape):
*   i        unit id                     t      calendar time
*   ttre     treatment onset (.=never)   y      continuous outcome (true effect 2 from onset)
*   ycount   count outcome               w      weights in [1,2)
*   sg       subgroup (0/1/2; .a on one unit per file where noted)
*   intensity     continuous dose (all units)     intensity_uz  dose, 0 for never-treated
clear all
set more off
local OUT "data"
cap mkdir "`OUT'"

cap program drop addcommon
program define addcommon
  * builds y/ycount/w/sg/intensity given i t ttre already defined
  gen te = 2*(t >= ttre & ttre < .)
  gen y = te + 0.05*i + rnormal()
  gen ycount = max(0, round(exp(0.25*te + 0.02*mod(i,7) + 0.3*rnormal())))
  gen w = 1 + runiform()
  gen sg = mod(i, 3)
  gen intensity = 1 + mod(i, 5) + 0.5*runiform()
  gen intensity_uz = cond(ttre < ., intensity, 0)
  drop te
end

* ---- dta_path_nonever: staggered onsets, NO never-treated units --------------
* wooldid INTENTIONALLY rejects this design (r459) unless safety(off) is used with a
* continuous treatment; for binary designs, truncate the sample window so the
* latest-treated cohorts become never-treated within the window. this dataset exists
* to verify the rejection stays informative, not to be made to run.
clear
set seed 20260721
set sortseed 20260721
set obs 30
gen i = _n
gen ttre = cond(_n <= 10, 4, cond(_n <= 20, 7, 9))
expand 10
bysort i: gen t = _n
addcommon
save "`OUT'/dta_path_nonever", replace

* ---- dta_path_nocontrols: every unit treated at the SAME time ----------------
* (no controls of any kind; used for error-path testing - must fail informatively)
clear
set seed 20260722
set sortseed 20260722
set obs 20
gen i = _n
gen ttre = 4
expand 8
bysort i: gen t = _n
addcommon
save "`OUT'/dta_path_nocontrols", replace

* ---- dta_path_nopre: one cohort treated at the first sample period -----------
clear
set seed 20260723
set sortseed 20260723
set obs 24
gen i = _n
gen ttre = cond(_n <= 8, 1, cond(_n <= 16, 5, .))
expand 8
bysort i: gen t = _n
addcommon
replace sg = .a if i == 24
save "`OUT'/dta_path_nopre", replace

* ---- dta_path_nopost: one cohort's onset falls after the sample ends ---------
* (ttre=9 with T=8: treated out-of-sample; default handling = control-like;
*  interacts with oostreatedcontrols())
clear
set seed 20260724
set sortseed 20260724
set obs 24
gen i = _n
gen ttre = cond(_n <= 8, 9, cond(_n <= 16, 5, .))
expand 8
bysort i: gen t = _n
addcommon
save "`OUT'/dta_path_nopost", replace

* ---- dta_path_relyrgap: cohorts missing specific relative years --------------
* cohort A (ttre=5): relyr -1 (t=4) absent for ALL its units - the default
*   reference period does not exist for this cohort
* cohort B (ttre=4): relyr +2 (t=6) absent for ALL its units - a mid-window
*   post-period hole
clear
set seed 20260725
set sortseed 20260725
set obs 30
gen i = _n
gen ttre = cond(_n <= 10, 5, cond(_n <= 20, 4, .))
expand 10
bysort i: gen t = _n
addcommon
drop if ttre == 5 & t == 4
drop if ttre == 4 & t == 6
save "`OUT'/dta_path_relyrgap", replace

* ---- dta_path_onepercell: exactly ONE observation per ttre-year cell ---------
* (single-unit cohorts: 4 treated cohorts of one unit each + 2 never-treated)
clear
set seed 20260726
set sortseed 20260726
set obs 6
gen i = _n
gen ttre = cond(_n == 1, 4, cond(_n == 2, 5, cond(_n == 3, 6, cond(_n == 4, 7, .))))
expand 10
bysort i: gen t = _n
addcommon
save "`OUT'/dta_path_onepercell", replace

* ---- dta_path_imbalanced: units span 2-12 periods; + fe() singleton inducer --
* extraFE is constant except unique values on three observations, so fe(extraFE)
* makes reghdfe drop those rows as singletons (e(sample) shrinkage, B10 class)
clear
set seed 20260727
set sortseed 20260727
set obs 30
gen i = _n
gen ttre = cond(mod(_n,3) == 0, 4, cond(mod(_n,3) == 1, 8, .))
gen span = 2 + mod(_n*7, 11)
expand 12
bysort i: gen t = _n
drop if t > span
drop span
addcommon
gen extraFE = 1
replace extraFE = 100 + i if inlist(i, 3, 6) & t == 2
replace sg = .a if i == 29
save "`OUT'/dta_path_imbalanced", replace

* ---- dta_path_sep: poisson separation material -------------------------------
* one never-treated unit with all-zero counts (unit-level separation) and one
* treated cohort-period cell forced to zero
clear
set seed 20260728
set sortseed 20260728
set obs 24
gen i = _n
gen ttre = cond(_n <= 8, 4, cond(_n <= 16, 6, .))
expand 8
bysort i: gen t = _n
addcommon
replace ycount = 0 if i == 20
replace ycount = 0 if ttre == 6 & t == 7
save "`OUT'/dta_path_sep", replace

dis "MKPATHOLOGICAL_DONE"
exit, clear
