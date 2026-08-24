# Preregistration — AME standard-error diagnosis grid (2026-08-24)

Committed while the MC runs and BEFORE any MC result is read (the per-rep CSVs are
unread at commit time; the deterministic structure pass below was run first, by
design, to validate case specs). Purpose (author directive): classify every case of
the continuous-treatment ("ame") SE problem into exactly one of four buckets —
(i) correct SEs with normal desirable properties; (ii) correct SEs with totally
undesirable properties (degenerate/one value); (iii) a correct SE is possible but the
code posts something incorrect; (iv) no correct SE is estimable in principle yet
something is posted — using MC truth, locked theory predictions, and multi-angle
adjudication. wooldid.ado is NOT modified in this exercise.

## Grid

DGPs (30 cohorts, 6 units/cohort, t=1..12, unitid distinct from i, y = cohort FE +
.05t + noise, dose NEVER enters y): d1 unit-level dose (varies across units within
cells, constant within unit); d2 within-unit time-varying dose; d3 cohort-constant
dose; d1c = d1 + common cohort-by-period error shock; d4 all-treated staggered
(onsets 5/7/9) with time-varying dose (the TRUE contreatwithin shape); d5 = d4 with
unit-constant dose. Options: contreatcontrols {basic, byCohort, byTime},
contreatwithin arm = basic+contreatwithin+safety(off); VCE {cluster(i),
cluster(unitid), robust}; window esprelength(4) espostlength(4) esfixedbaseperiod.

MC cases (T=256 each, seeds recorded per rep, driver mc_amese.do):
d1 x 3cc x 3vce, d2 x 3cc x 3vce, d4 x within x 3vce, d5 x within x 3vce,
d1c x basic x {ci, rb} — 26 cases. Structure-only (deterministically degenerate or
duplicate, demoted per the structure pass): d3 x basic x 3vce (all SEs missing +
identification cautions under every VCE); d1/d2 x within (bitwise identical to
d1/d2 x basic — contreatwithin is estimation-inert when controls exist, itself a
recorded observation); d1c_basic_ci total-SE-wipeout channel (contreat-specific:
binary att SEs are fine on the same data) is assigned to the code-auditor.

## Structure-pass facts (recorded pre-MC; log: jobs tmp amese/wooldid.log + struct2.log)

- rank-1 identical-SE collapse (silent): d1_basic_ci, d2_basic_ci, d5_within_ci.
- TOTAL SE wipeout WITH honest caution ("standard errors were not identified for 20
  of 20 ... no e(V) posted"): d1_bycoh_ci, d2_bycoh_ci, d4_within_ci, d1c_basic_ci.
- Full-rank horizon-varying SEs: all cu and rb cases; d1_bytime_ci and d2_bytime_ci
  (cluster(i) + byTime posts FULL-rank V — whether those SEs are CORRECT is exactly
  what the MC adjudicates).

## Preregistered criteria

Primary metric per case and per quantity (main ame, pretreatment ame, each of the 8
ES horizons): ratio = mean(reported SE) / SD_true, SD_true = cross-rep SD of the
point estimate over the rc-0 reps. With T=256 the SD estimate carries ~4.4% noise, so:
ratio in [0.85, 1.15] = "approximately correct"; outside [0.70, 1.40] = "wrong";
between = FLAGGED (theory adjudicates). A case is classified on the pattern across
its 10 quantities, not a single ratio. Secondary: SE dispersion across horizons
(identical-to-14-digits = collapse signature), V_es_ame eigen range, rc/missing
rates, and the cross-rep SD of the POINT estimates themselves (a degenerate-estimand
case is a different diagnosis from an untrackable-variance case).

## Theory protocol

Two econometrician agents (incidental-parameters lens; design-based lens) produce
LOCKED per-case predictions + bucket classifications + remedy prescriptions before
seeing any MC output (they cannot read the CSVs; their reports are committed on
arrival). At adjudication, a third theorist reads the pooled results BLIND to the
predictions; a code-path auditor traces exactly what wooldid computes per case
(bucket (iii) vs (iv) hinges on it, incl. the d1c wipeout channel and the
byTime+cluster(i) partial-collapse question); an adversarial referee attacks the
final classification. Deviations from this protocol are documented in the output.

## Output

wooldid_standing_issues.md family-A section updated with: per-case bucket table,
the cross-case mechanism narrative, and detailed fix recommendations mapped to the
author's three response types (say-not-estimable / user-misspecification-steer /
code-fix), each with its re-bless blast radius. No code changes.
