* Copied from the colleagues' validation bundle
* /Users/hegland/Documents/Research/wooldid_dev/wooldid_external/wooldid_test_comparison_code/tests/lib/fixture_common.do
* on 2026-08-18; pure Stata, deterministic, no RNG (config.do sets seeds but the
* fixture formulas draw nothing). Content unmodified except this header.
version 17.0

capture program drop wdbt_make_master_fixture
program define wdbt_make_master_fixture
    version 17.0
    syntax, Reps(integer)

    if `reps' < 3 {
        display as error "wdbt_make_master_fixture requires at least 3 observations per cohort-period."
        error 198
    }

    local cohort_count 12
    local period_count 12
    local observation_count = `cohort_count' * `period_count' * `reps'

    clear
    quietly set obs `observation_count'
    generate int cohort = floor((_n - 1) / (`period_count' * `reps')) + 1
    generate int unit_in_cohort = floor(mod(_n - 1, `period_count' * `reps') / `period_count') + 1
    generate int t = mod(_n - 1, `period_count') + 1
    generate long person_id = 100 * cohort + unit_in_cohort
    generate long repeated_cross_section_id = 100000 * cohort + 100 * t + unit_in_cohort
    generate long wdbt_obs_id = 1000000 * cohort + 1000 * unit_in_cohort + t

    generate double ttre = cond(cohort <= 3, 6, ///
        cond(cohort <= 6, 7, cond(cohort <= 9, 8, .)))
    generate double ttre_single = cond(cohort <= 9, 7, .)
    generate byte ever_treated = ttre < .
    generate byte treated = ever_treated & t >= ttre
    generate double relative_time = t - ttre if ever_treated
    generate byte ever_treated_single = ttre_single < .
    generate byte treated_single = ever_treated_single & t >= ttre_single
    generate double relative_time_single = t - ttre_single if ever_treated_single
    generate double ttre_oos = cond(cohort <= 3, 6, ///
        cond(cohort <= 6, 8, cond(cohort <= 8, 13, ///
        cond(cohort == 9, 14, .))))
    generate byte treated_oos = ttre_oos < . & t >= ttre_oos
    generate double relative_time_oos = t - ttre_oos if ttre_oos < .

    // Subgroups are cohort-level, as required by the authoritative help file.
    generate byte sg_three = mod(cohort - 1, 3) + 1
    generate byte sg_one = 1 if inlist(cohort, 1, 4, 7, 10)
    generate byte sg_ttre = cond(ttre == 6, 1, cond(ttre == 7, 2, ///
        cond(ttre == 8, 3, .)))

    generate long fine_i = person_id
    generate long mixed_i = cond(cohort <= 4, 10000 + person_id, cohort)
    // The rich mixed-grain identifier has cell sizes 1, 2, 3, and 8 when
    // the fixture is built with eight observations per cohort-period.
    generate long mixed_i_rich = .
    replace mixed_i_rich = 10000 * cohort + unit_in_cohort ///
        if inrange(cohort, 1, 3)
    replace mixed_i_rich = 10000 * cohort + ceil(unit_in_cohort / 2) ///
        if inrange(cohort, 4, 6)
    replace mixed_i_rich = 10000 * cohort + ///
        cond(unit_in_cohort <= 3, 1, cond(unit_in_cohort <= 6, 2, 3)) ///
        if inrange(cohort, 7, 9)
    replace mixed_i_rich = 10000 * cohort + 1 if cohort >= 10
    generate byte categorical_control = mod(cohort + unit_in_cohort, 3)
    generate double x1 = 0.15 * cohort + 0.04 * t + ///
        0.03 * unit_in_cohort + 0.02 * sin(cohort + t)
    generate double x2 = cos(0.7 * cohort + 0.4 * t + 0.2 * unit_in_cohort)
    generate double interactive_x = 0.25 * cohort + 0.07 * unit_in_cohort + ///
        0.03 * t + 0.05 * sin(unit_in_cohort * t)

    // dose is constant within a panel unit but dense within each cohort-period.
    generate double dose = 0.35 + 0.055 * cohort + 0.021 * unit_in_cohort + ///
        0.006 * mod(cohort * unit_in_cohort, 7)
    generate double dose_vary = dose + 0.012 * (t - 6) * ///
        cond(mod(unit_in_cohort, 2), 1, -1)
    generate double dose_zero = cond(treated, dose, 0)
    generate byte custom_lattice = cond(dose < 0.75, 1, cond(dose < 1.05, 2, 3))
    // Common-support dose variables deliberately have the same support in
    // every cohort-period cell.  Use these for higher-order and CFX tests.
    generate double dose_common = 0.30 + 0.05 * unit_in_cohort
    generate double dose_common_vary = dose_common + ///
        0.008 * (t - 6) * cond(mod(unit_in_cohort, 2), 1, -1)
    generate double dose_common_zero = cond(treated, dose_common, 0)
    generate byte custom_lattice_common = ///
        ceil(3 * unit_in_cohort / `reps')

    generate double anticipation = cond(relative_time == -3, 0.06, ///
        cond(relative_time == -2, -0.09, cond(relative_time == -1, 0.03, 0)))
    replace anticipation = 0 if !ever_treated
    generate double anticipation_single = cond(relative_time_single == -3, 0.06, ///
        cond(relative_time_single == -2, -0.09, ///
        cond(relative_time_single == -1, 0.03, 0)))
    replace anticipation_single = 0 if !ever_treated_single
    generate double binary_effect = cond(treated, ///
        0.75 + 0.10 * relative_time + 0.045 * cohort, 0)
    generate double binary_effect_single = cond(treated_single, ///
        0.75 + 0.10 * relative_time_single + 0.045 * cohort, 0)
    generate double continuous_effect = cond(treated, ///
        0.18 * dose + 0.035 * dose^2 + 0.012 * dose^3, 0)
    generate double continuous_effect_single = cond(treated_single, ///
        0.18 * dose + 0.035 * dose^2 + 0.012 * dose^3, 0)
    generate double deterministic_noise = ///
        0.11 * sin(1.7 * cohort + 0.9 * unit_in_cohort + 0.8 * t) + ///
        0.04 * cos(0.3 * cohort * t + 0.5 * unit_in_cohort)

    generate double y = 4.5 + 0.17 * cohort + 0.035 * unit_in_cohort + ///
        0.09 * t + 0.004 * t^2 + 0.22 * x1 - 0.08 * x2 + ///
        0.16 * dose + anticipation + binary_effect + continuous_effect + ///
        deterministic_noise
    generate double y_single = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + 0.16 * dose + anticipation_single + ///
        binary_effect_single + continuous_effect_single + deterministic_noise
    generate double y_cont_common = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + 0.16 * dose_common + anticipation + ///
        binary_effect + cond(treated, 0.18 * dose_common + ///
        0.035 * dose_common^2 + 0.012 * dose_common^3, 0) + ///
        deterministic_noise
    generate double y_cont_single = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + 0.16 * dose_common + ///
        anticipation_single + binary_effect_single + ///
        cond(treated_single, 0.18 * dose_common + ///
        0.035 * dose_common^2 + 0.012 * dose_common^3, 0) + ///
        deterministic_noise
    generate double y_cont_vary = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + 0.16 * dose_common_vary + ///
        anticipation + binary_effect + cond(treated, ///
        0.18 * dose_common_vary + 0.035 * dose_common_vary^2 + ///
        0.012 * dose_common_vary^3, 0) + deterministic_noise
    generate double y_cont_zero = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + 0.16 * dose_common_zero + ///
        anticipation + binary_effect + cond(treated, ///
        0.18 * dose_common_zero + 0.035 * dose_common_zero^2 + ///
        0.012 * dose_common_zero^3, 0) + deterministic_noise
    generate double y_interactive = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + 0.12 * interactive_x + ///
        0.16 * dose_common + anticipation + binary_effect + ///
        cond(treated, 0.18 * dose_common + 0.04 * dose_common * ///
        interactive_x, 0) + deterministic_noise
    generate double y_oos = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + 0.004 * t^2 + ///
        0.22 * x1 - 0.08 * x2 + cond(treated_oos, ///
        0.70 + 0.09 * relative_time_oos, 0) + deterministic_noise
    generate double y_trend = 4.5 + 0.17 * cohort + ///
        0.035 * unit_in_cohort + 0.09 * t + ///
        0.0015 * cohort * t + anticipation + binary_effect + ///
        deterministic_noise
    generate double y_summary = 2.0 + 0.08 * cohort + ///
        0.025 * unit_in_cohort + 0.05 * t + anticipation + ///
        0.65 * binary_effect + 0.20 * deterministic_noise
    generate double y_positive = exp(0.35 + 0.025 * cohort + 0.012 * t + ///
        0.03 * x1 + 0.08 * treated + ///
        0.025 * cond(treated, relative_time, 0) + ///
        0.035 * dose * treated + 0.02 * deterministic_noise)
    generate double y_positive_common = exp(0.35 + 0.025 * cohort + ///
        0.012 * t + 0.03 * x1 + 0.08 * treated + ///
        0.025 * cond(treated, relative_time, 0) + ///
        0.035 * dose_common * treated + 0.02 * deterministic_noise)
    generate long y_count = max(1, round(8 * y_positive))
    generate byte y_binary = mod(unit_in_cohort + cohort + t, 4) == 0
    replace y_binary = 1 if treated & ///
        mod(unit_in_cohort + cohort + t, 4) == 1
    generate byte y_count_zero = mod(unit_in_cohort + 2 * cohort + t, 4)
    replace y_count_zero = y_count_zero + 1 if treated & ///
        mod(unit_in_cohort + cohort + t, 3) == 0
    // With twelve observations per cell the untreated/treated means are
    // exactly .25/.50, making count-family oracle checks deterministic.
    generate byte y_count_oracle = cond(treated, ///
        unit_in_cohort <= 6, unit_in_cohort <= 3)
    generate double analytic_weight = 0.8 + 0.05 * mod(cohort + unit_in_cohort, 5)
    generate double probability_weight = 1 + 0.10 * mod(2 * cohort + unit_in_cohort, 4)
    generate double custom_weight_multiplier = 0.9 + 0.04 * mod(cohort + t, 4)

    label variable cohort "Coarse cohort identifier"
    label variable t "Time period"
    label variable ttre "Staggered treatment-onset period"
    label variable ttre_single "Common treatment-onset period"
    label variable ttre_oos "Treatment onset including out-of-sample onsets"
    label variable dose "Continuous treatment assigned at panel-unit level"
    label variable dose_zero "Continuous treatment equal to zero outside treatment"
    label variable dose_common "Continuous treatment with common cell support"
    label variable dose_common_zero "Common-support dose equal to zero untreated"
    label variable y "Deterministic nondegenerate outcome"
    label variable y_single "Outcome matched to the common-onset assignment"
    label define subgroup_three 1 "SG A" 2 "SG B" 3 "SG C"
    label values sg_three subgroup_three

    sort cohort unit_in_cohort t
    tempvar lattice_min lattice_max
    bysort cohort t: egen byte `lattice_min' = min(custom_lattice_common)
    bysort cohort t: egen byte `lattice_max' = max(custom_lattice_common)
    assert `lattice_min' == 1
    assert `lattice_max' == 3
    drop `lattice_min' `lattice_max'
    if `reps' == 12 {
        tempvar count_cell_mean
        bysort cohort t: egen double `count_cell_mean' = mean(y_count_oracle)
        assert abs(`count_cell_mean' - 0.25) < 1e-12 if !treated
        assert abs(`count_cell_mean' - 0.50) < 1e-12 if treated
        drop `count_cell_mean'
    }
    compress
    isid wdbt_obs_id
    assert ttre == ttre[1] if cohort == 1
end

capture program drop wdbt_make_survey_rcs_fixture
program define wdbt_make_survey_rcs_fixture
    version 17.0

    local state_count 24
    local group_count 2
    local wave_count 8
    local slots_per_cell 6
    local observation_count = `state_count' * `group_count' * ///
        `wave_count' * `slots_per_cell'

    clear
    quietly set obs `observation_count'
    generate byte state = floor((_n - 1) / ///
        (`group_count' * `wave_count' * `slots_per_cell')) + 1
    generate byte group = floor(mod(_n - 1, ///
        `group_count' * `wave_count' * `slots_per_cell') / ///
        (`wave_count' * `slots_per_cell')) + 1
    generate byte t = floor(mod(_n - 1, ///
        `wave_count' * `slots_per_cell') / `slots_per_cell') + 1
    generate byte respondent_slot = mod(_n - 1, `slots_per_cell') + 1
    generate int state_group = 10 * state + group
    generate byte group_time = 10 * group + t
    generate long respondent_id = 100000 * state + 10000 * group + ///
        100 * t + respondent_slot
    generate long repeated_cross_section_id = respondent_id
    generate long wdbt_obs_id = respondent_id

    generate double ttre = cond(state <= 6, 5, ///
        cond(state <= 12, 6, cond(state <= 18, 7, .)))
    generate byte ever_treated = ttre < .
    generate byte treated = ever_treated & t >= ttre
    generate double relative_time = t - ttre if ever_treated

    generate double state_control = 0.08 * state + 0.03 * group
    generate double respondent_control = ///
        0.12 * (respondent_slot - 3.5) + 0.02 * mod(state + t, 3)
    generate double survey_weight = 1 + ///
        0.10 * min(respondent_slot, 7 - respondent_slot) + ///
        0.02 * mod(state, 3)
    generate double survey_weight10 = 10 * survey_weight
    generate double paired_noise = 0.04 * (respondent_slot - 3.5)
    // This common state-by-wave shock is deliberately not additive in the
    // absorbed state-group and group-wave effects.  It yields nonzero cluster
    // scores while cancelling exactly in the between-group effect contrast.
    generate double rcs_cell_shock = ///
        0.07 * sin(0.7 * state + 1.1 * t) + ///
        0.03 * cos(0.2 * state * t)
    generate double rcs_treatment_effect = cond(treated, ///
        0.65 + 0.025 * mod(state, 5) + 0.04 * relative_time + group - 1, 0)
    generate double y_rcs = 3.0 + 0.06 * state + 0.12 * group + ///
        0.05 * t + 0.15 * state_control + paired_noise + rcs_cell_shock + ///
        rcs_treatment_effect

    // Missingness removes a weight-symmetric respondent pair and leaves four
    // observations in every affected state-group-period cell.
    replace y_rcs = . if treated & inlist(respondent_slot, 1, 6) & ///
        mod(state + t, 3) == 0

    sort state group t respondent_slot
    tempvar cell_n onset_min onset_max
    bysort state_group t: egen byte `cell_n' = total(!missing(y_rcs))
    bysort state_group: egen double `onset_min' = min(ttre)
    bysort state_group: egen double `onset_max' = max(ttre)
    assert `cell_n' >= 4
    assert `onset_min' == `onset_max' | ///
        (missing(`onset_min') & missing(`onset_max'))
    assert survey_weight10 == 10 * survey_weight
    drop `cell_n' `onset_min' `onset_max'
    compress
    isid respondent_id
    isid state_group t respondent_slot
end
