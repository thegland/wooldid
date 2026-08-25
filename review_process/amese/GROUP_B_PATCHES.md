# Group B Patches — the byCohort point-identification screen (author-ruled 2026-08-25)

Separate revert pool from the Group A Patches (which the author is testing). Group B
implements the ruled hardened variant of screen option 1 from the byCohort
solver-normalization diagnosis:

- EXACT case (a treated cohort whose untreated/reference-pool observations carry ZERO
  variation in the continuous treatment): hard refusal, r(459), naming the smallest
  offending cohort, advising examination of that cohort's dose values and dropping it
  from the estimation sample — the pre-treatment-period sufficiency errors' sibling,
  one layer down (rows can exist while the identifying VARIATION does not).
- THIN case (uncovered-pool dose sd below 1/1000 of the sample dose sd): a
  computational-fragility caution naming the count and smallest cohort, suggesting a
  larger regtol() (tighter convergence tolerance) or dropping the cohort(s).
- Always-on when contreatcontrols(byCohort)+contreat are active (any VCE, fastnose,
  impvalues): POINT identification is at stake, not variances. Runs pre-estimation.
- Implementation: wooldid_amescreen mode(bycoh); grouping keys on the ORIGINAL cohort
  variable among treated rows (the `_i' cell coding pools all reference rows into
  cohort code 0 — the first-cut icell grouping false-fired universally and was
  corrected before commit); one gegen pass, O(N).

| sha | subject | touches |
|---|---|---|
| (this commit) | the screen + sthlp + cert B1-B3 | wooldidinstall/wooldid.ado, wooldidinstall/wooldid.sthlp, test/test_amescreen.do |

Verification: cert AMESCREEN_ALL_PASS (13 scenarios: A1-A10 regression-clean incl.
A5's byCohort refit, B1 refuses-459-naming-cohort, B2 thin-caution-completes,
B3 healthy-silent); battery 34/34 rtol-0 unchanged; suite unreachable except
ext_cont_ext_controls_cohort_time (probed directly on its fixture: screen silent,
rc 0 — zero suite surface change, no re-bless, no ledger rows).

Revert: git revert of this single commit; no baselines or ledger rows involved.
