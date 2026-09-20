# Predicted battery classes, non-number-moving batch — written BEFORE the compare was opened

Target: `test/runB2` against `test/baseline`, which `git log` confirms was re-blessed at
`1c09efa` (the R1-R19 slate). So every group below is **this batch's**, not R-batch residue.

Config inventory in `runbattery.do` (63 `dumpres` calls): `contreatwithin` 0, `gteffects` 5,
`boottstat` 13, `unconditionalse` 1, `pretrendpower` 2.

| class | what | where | mechanism |
|---|---|---|---|
| D1 | `macro vce` robust -> unconditional; `macro vcetype` Robust -> Unconditional; `macro vce_stage1` ADDED (= robust) | the 1 config at runbattery.do:127 (`att robust unconditionalse`) | VCV-3 |
| D2 | `scalar gt_ncontained` ADDED | the `gteffects` configs whose VCE can see containment: `gt_basic` (cluster(i)) and `boot_ts_gt` (cluster(i)); NOT `fn_gt`, which is `fastnose` and gated out | GTE-4 |
| D3 | `matnames`/`mat boot_jointWobs` ADDED on the studentized path, and `boot_advWobs` ADDED where an adversarial test ran | the `boottstat` configs that reach a joint test | BSAVE-1 |
| D4 | `scalar roth_cellmax_*` ADDED | at most the 2 `pretrendpower` configs, and only if the pre-period block is rank-deficient there | PTP Proposal 2 |

## Correction made before the compare was opened

D4 first said this class **could not appear**, on the reasoning that `e(scalars)` does not
enumerate `ereturn hidden` results. **That reasoning is wrong**, and it was checked against
dumps already on disk rather than assumed: `test/baseline` carries 16 `scalar|df_r_analytic|`
lines and 16 `matnames|boot_ptails|` lines, both of which are posted `ereturn hidden`. So
Stata's `e(scalars)`/`e(matrices)` do enumerate hidden results and the dumper sees them. D4 is
therefore a live class, not an expected absence. Recorded here so that either outcome —
present or absent — is pre-explained rather than rationalised afterwards.

## Expected ABSENCES (recorded because absence is information)

* **G2** cannot fire: no battery config specifies `contreatwithin`, so the `safety` ->
  `safetyinput` rekey of the identifier-domain assert has no arm here.
* **LAT L1**, **SMP-2/BSAVE-3**, **G1b**, **FIX-PERM-3** and **G4** change displayed text, an
  error path, or a display format only. None touches `e()`, so none may appear.
* **DOC-1 / IMP P2 / GTE-5 / BAND-3 / BAND-1b / L4** are help-file only.

## The falsifier

Anything outside D1-D4 is unexplained and must be adjudicated per the project rule: name a new
cell CLASS in the comparer, pristine-twice, never blanket the matrix. In particular **no matrix
cell, no existing scalar, no `rc` and no `varlist` may move anywhere.**
