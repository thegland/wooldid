# fastnose-refute — skeptic's notes

Build asserted in every log: `which wooldid` ->
`/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`, `*! version 1.4 2026-08-19`.
All data generated in-file; no project test suite, no project data, no prober do-file run.

## Probes (all self-contained; each generates its own panel)
| file | what | result |
|---|---|---|
| `r1_fn03.do` / `.log` | MY panel (105 units, T=2010-2019, onsets 2013/2016/2019 + never, dose varies within cohort, subgroup = late-adopter region). contreat x subgroups x ES, analytic vs fastnose, plus a no-subgroups control | FN-03 reproduced on my own geometry |
| `r2_fn01_fn02.do` / `.log` | MY binary panel. permute+permanalyticse (A), +fastnose (B), plain permute (C); gteffects analytic vs fastnose | FN-01 reproduced; FN-02 behaviour reproduced but its doc claim is false |
| `r3_escalate.do` / `.log` | same panel as r1 under permute() default, permute+permanalyticse, bootstrap() | escalation: P rc0 with p-values on zeros; PA rc **459**; BT rc **459** |
| `r4_controls.do` / `.log` | control matrix pinning the interaction | B1 nosg rc0, B2 no-ES rc0, B3 BINARY rc0, B4 contreat+sg+ES rc **459**, B5 same w/ robust rc **459** |
| `r5_mechanism.do` / `.log` | first (failed) attempt at a margins-level marker demo | inconclusive — design was fully estimable; kept as a negative record |
| `r6_omarker.do` / `.log` | margins-level demo, no wooldid | **mechanism proven**: `df()` -> colnames `1bno.t 2o.t 3o.t 4.t 5.t 6.t` (3 markers); `nose` -> `1bn.t 2.t 3.t 4.t 5.t 6.t` (0 markers) |

## Verdicts
- **FASTNOSE-03 CONFIRMED at high (NOT raised).** Reproduced on my own data
  in both `contreatcontrols(basic)` and `untreatedZero`; control without `subgroups()` is
  identical in both modes. Mechanism the prober could only infer is now PROVEN at the
  margins level by r6. Escalations they never reached: the `permute()` FAST-MODE DEFAULT
  (no `fastnose` typed) publishes a whole subgroup ame event study the analytic run skips,
  WITH permutation p-values and shift intervals, and `bootstrap()` / `permute()
  permanalyticse` ABORT r(459) on the same spec. Filed as SKEP-01/SKEP-02 so the shared
  root cause is not double-counted. I checked and REJECTED the critical framing: on the
  POOLED path the analytic run also reports pre3/pre2/pre1 = 0, flagged by e(es_ndegen) —
  documented design (sthlp:1325). The subgroup row set is the defect.
- **FASTNOSE-01 CONFIRMED at medium**, with two corrections in the prober's favour that
  they missed: sthlp:1234's engine-compatibility table says under `permute()` fastnose
  "*is* the default ... `permanalyticse` opts back in", so this is also a doc-vs-code
  divergence; and the run does NOT print the fast-mode header line, so nothing anywhere
  says fast mode is in force.
- **FASTNOSE-02 REFUTED.** sthlp:1149-1150 says verbatim "Fully available under
  {it:fastnose} (the regression's own VCE is intact even though the {it:margins} stage
  skips variances), making {it:gteffects} the one inference surface that survives fast
  mode." The claim "the gteffects paragraph never mentions fastnose" is false. Behaviour
  reproduced (32x15 e(gteffects), byte-identical in %21x) but it is documented, not a gap.

## Free side-checks done while in here
- known_bydesign §7 item 11: `fastnose` DOES post `e(b_es_att)` (cols pre3..post1) and
  `e(numpre_att)`=3, and prints the Caution sthlp:1325 describes. sthlp:1325 keeping
  `fastnose` at the head of that list is ACCURATE. **No divergence — negative result.**
