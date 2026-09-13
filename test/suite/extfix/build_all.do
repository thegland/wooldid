* Copied from the colleagues' validation bundle
* /Users/hegland/Documents/Research/wooldid_dev/wooldid_external/wooldid_test_comparison_code/tests/fixtures/build_all.do
* on 2026-08-18; pure Stata, deterministic, no RNG.
* MINIMAL ADAPTATION (path plumbing only; fixture content untouched):
*   - runs standalone with an output-dir argument:  do build_all.do <outdir>
*     (default outdir: ./generated). No workspace layout assumed.
*   - expects its two siblings config.do and fixture_common.do in the CURRENT
*     working directory: run it with cwd = this extfix directory.
*   - the original resolved everything off a <workspace> root
*     (<ws>/tests/config.do, <ws>/tests/lib/fixture_common.do,
*      <ws>/tests/fixtures/generated).
* version 17.0 pin kept: verified to run under Stata 19 MP.
version 17.0
set more off

args outdir
local outdir = subinstr(`"`outdir'"', "\", "/", .)
if `"`outdir'"' == "" local outdir "generated"
capture mkdir `"`outdir'"'
local generated `"`outdir'"'

do config.do
do fixture_common.do

local master_keep cohort unit_in_cohort t person_id ///
    repeated_cross_section_id wdbt_obs_id ///
    ttre ttre_single ttre_oos ///
    ever_treated treated relative_time ///
    ever_treated_single treated_single relative_time_single ///
    treated_oos relative_time_oos ///
    sg_one sg_three sg_ttre ///
    fine_i mixed_i mixed_i_rich ///
    categorical_control x1 x2 interactive_x ///
    dose dose_vary dose_zero custom_lattice ///
    dose_common dose_common_vary dose_common_zero custom_lattice_common ///
    anticipation anticipation_single ///
    binary_effect binary_effect_single ///
    continuous_effect continuous_effect_single deterministic_noise ///
    y y_single y_cont_common y_cont_single y_cont_vary y_cont_zero ///
    y_interactive y_oos y_trend y_summary ///
    y_positive y_positive_common y_count y_binary ///
    y_count_zero y_count_oracle ///
    analytic_weight probability_weight custom_weight_multiplier

wdbt_make_master_fixture, reps(3)
preserve
keep `master_keep'
char _dta[wdbt_fixture_id] "binary_balanced"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/binary_balanced.dta"', replace
restore

preserve
drop if (cohort == 2 & inlist(t, 1, 4, 9)) | ///
    (cohort == 5 & inlist(t, 2, 3, 11)) | ///
    (cohort == 8 & inlist(t, 1, 6)) | ///
    (cohort == 11 & inlist(t, 5, 10))
keep `master_keep'
char _dta[wdbt_fixture_id] "binary_mixed_imbalanced"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/binary_mixed_imbalanced.dta"', replace
restore

preserve
replace cohort = fine_i
isid cohort t
keep `master_keep'
char _dta[wdbt_fixture_id] "binary_fine"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/binary_fine.dta"', replace
restore

wdbt_make_master_fixture, reps(8)
preserve
tempvar mixed_cell_n mixed_cell_tag
bysort mixed_i_rich t: generate byte `mixed_cell_n' = _N
egen byte `mixed_cell_tag' = tag(mixed_i_rich t)
assert inlist(`mixed_cell_n', 1, 2, 3, 8)
quietly count if `mixed_cell_tag' & `mixed_cell_n' == 1
assert r(N) > 0
quietly count if `mixed_cell_tag' & `mixed_cell_n' == 2
assert r(N) > 0
quietly count if `mixed_cell_tag' & `mixed_cell_n' == 3
assert r(N) > 0
quietly count if `mixed_cell_tag' & `mixed_cell_n' == 8
assert r(N) > 0
drop `mixed_cell_n' `mixed_cell_tag'
keep `master_keep'
char _dta[wdbt_fixture_id] "binary_mixed_grain"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/binary_mixed_grain.dta"', replace
restore

wdbt_make_master_fixture, reps(4)
preserve
keep `master_keep'
char _dta[wdbt_fixture_id] "continuous_control_compact"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/continuous_control_compact.dta"', replace
restore

preserve
keep if ttre < .
keep `master_keep'
char _dta[wdbt_fixture_id] "continuous_within_compact"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/continuous_within_compact.dta"', replace
restore


wdbt_make_master_fixture, reps(12)
preserve
keep `master_keep'
char _dta[wdbt_fixture_id] "continuous_control_dense"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/continuous_control_dense.dta"', replace
restore

preserve
keep if ttre < .
keep `master_keep'
char _dta[wdbt_fixture_id] "continuous_within_dense"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/continuous_within_dense.dta"', replace
restore

preserve
keep `master_keep'
char _dta[wdbt_fixture_id] "continuous_untreated_zero"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/continuous_untreated_zero.dta"', replace
restore

preserve
drop if (cohort == 2 & inlist(t, 1, 9)) | ///
    (cohort == 5 & inlist(t, 2, 11)) | ///
    (cohort == 8 & inlist(t, 1, 10)) | ///
    (cohort == 11 & inlist(t, 5, 10))
bysort cohort t: assert _N == 12
tempvar reference_n
bysort cohort: egen byte `reference_n' = total(t == ttre - 2)
assert `reference_n' == 12 if ttre < .
drop `reference_n'
keep `master_keep'
char _dta[wdbt_fixture_id] "continuous_imbalanced_dense"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/continuous_imbalanced_dense.dta"', replace
restore

wdbt_make_survey_rcs_fixture
keep state group t state_group group_time ///
    ttre ever_treated treated relative_time ///
    respondent_id repeated_cross_section_id respondent_slot ///
    survey_weight survey_weight10 state_control respondent_control ///
    rcs_cell_shock rcs_treatment_effect y_rcs wdbt_obs_id
char _dta[wdbt_fixture_id] "survey_rcs"
char _dta[wdbt_formula_version] "wdbt-fixture-v2"
save `"`generated'/survey_rcs.dta"', replace

display as result "WDBT_FIXTURES_BUILT"
exit 0
