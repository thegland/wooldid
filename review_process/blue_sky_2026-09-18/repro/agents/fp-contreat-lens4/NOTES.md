# fp-contreat-D — LENS 4 (documented-behaviour) adjudication

## Verdict (set from primary sources BEFORE the confirmatory run): REFUTED, severity none

### Three primary sources, all verified on the pinned build
1. **The behaviour is documented verbatim, with the identical example.**
   `build/wooldid.sthlp:1052` (a `{pmore}` inside the `fe()` entry):
   "This holds on the default estimation path. Under customvce() or unconditionalse,
   estimation runs through regress ... and a BARE token there is read as a CONTINUOUS
   variable: fe(state) fits a linear term in the state code rather than a set of state
   indicators. Write the prefix explicitly - fe(i.state) - whenever those options are in
   use, or on any specification you may later run both ways."
2. **It was ruled, and the ruling was to document rather than error.**
   `review_process/rulings_2026-09-07.md:42`, section **B. Documentation-only rulings**:
   | BS-1 | **Document the difference** | state that `fe()` tokens are categorical on the
   reghdfe path and continuous under `customvce()`/`unconditionalse`, so a bare token there
   fits a linear term in the identifier. |
   The prober's EXPECTED ("wooldid knows at parse time the combination is invalid -- an
   error would cost nothing") is precisely the alternative BS-1 rejected.
3. **The ruling was BUILT, not merely decided** (the coordinator note warns that "ruled"
   != "implemented"; that gap is closed here):
   - `git log -S"BARE token" -- wooldidinstall/wooldid.sthlp` -> `e26b3dd` 2026-09-09
     "Ruled batches one and two: the documentation set, and the refusals and screening order"
   - md5 of `wooldidinstall/wooldid.sthlp` == md5 of the pinned `build/wooldid.sthlp`
     (3cda877961315e9dde2a00415499da14); same for the .ado (444892481be8d50ec962a102d6abef51).
   The runtime caution is `ado:4637-4638` and predates all of this (`fc17051`, 2023-08-15).

### The one thing neither the prober nor the skeptic showed, and the pass criterion
Both showed only that the two spellings DIFFER (R2 0.96 vs 0.73; ame -0.178 vs 1.370).
Neither showed that the DOCUMENTED REMEDY RECOVERS the absorbed answer. That is the only
thing that could still flip this finding.

**Pass criterion, fixed before reading the log** (`L4_fe_remedy.do`, prober's own probe08
fixture, y exactly linear so the absorbed fit is exact and TRUTH_AME is closed-form):
- PASS (verdict stays REFUTED): `HELPFIX_ame` (fe(i.unit) + customvce), `CTRLFIX_ame`
  (controls(i.unit), the runtime caution's literal remedy) and `HELPFIXUNC_ame`
  (unconditionalse spelling) all match `REF_ame` (fe(unit) + cluster, reghdfe-absorbed)
  to rtol 1e-10. Point estimates only -- the plain-reg SE non-reproducibility is by design.
- FAIL (verdict becomes PARTIAL): any of those legs fails to recover `REF_ame`. A documented
  remedy that does not work is a live residual.
Artifact to IGNORE: y is deterministic, so REF fits exactly and the covariance is degenerate;
any all-zero-covariance caution or scrubbed SE here is the DGP, not known_bydesign sec.5.

### Residual (cosmetic, does not move the verdict)
The runtime caution "You must specify anything in fe as controls" reads as though fe()
contents are DROPPED on this path, when the code in fact passes them straight through as a
bare (continuous) token. The following clause "manually specifying any needed i. prefixes"
already flags the actual trap, and the help states it precisely. Wording nit.

## RESULTS — PASS on every leg. Verdict stands: REFUTED, severity none.

`L4_fe_remedy.do` / `L4_fe_remedy.log` (build asserted: `which wooldid` resolves to the
pinned .../bluesky/build/wooldid.ado). My TRUTH_AME scalar (2.270000000000 = mean of 0.9 + betatrue over
treated post rows) is NOT the quantity wooldid reports, and the gap is arithmetic, not a
defect: contreatcontrols(basic) puts a pooled c.dose main effect in the model, so the 0.9
baseline slope is a control and the reported ame is the average of betatrue alone. Checked
by hand: betatrue = 2.0 - 0.3*(t-ttre); cohort i=1 has t=5..10 (mean 7.5/6 = 1.25) over 60
rows, cohort i=2 has t=7..10 (mean 6.2/4 = 1.55) over 40 rows, so the row-weighted mean is
(60*1.25 + 40*1.55)/100 = 1.37 exactly -- which is REF_ame to all twelve printed decimals.
REF is therefore the exact-fit absorbed truth and is the comparison target for every leg.

| leg | spelling | att | ame | r2 |
|---|---|---|---|---|
| REF | `fe(unit) cluster(unit)` (reghdfe, absorbed) | 4.169863808849 | 1.370000000000 | 1.00000000 |
| BARE | `fe(unit) customvce(cluster unit)` | 4.085936283577 | **-0.178105713914** | 0.33211 |
| HELPFIX | `fe(i.unit) customvce(cluster unit)` | 4.169863808849 | **1.370000000000** | 1.0 |
| CTRLFIX | `controls(i.unit) customvce(cluster unit)` | 4.169863808849 | **1.370000000000** | 1.0 |
| HELPFIX_UNC | `fe(i.unit) unconditionalse cluster(unit)` | 4.169863808849 | **1.370000000000** | 1.0 |

- The BARE leg reproduces the prober's numbers exactly (ame -0.178105713914, r2 0.332).
- **The documented remedy recovers the absorbed answer bit-for-bit at 12 decimals**, on all
  three spellings: the help's `fe(i.unit)`, the runtime caution's `controls(i.unit)`, and the
  `unconditionalse` spelling. Pass criterion (rtol 1e-10) met with room to spare — these are
  exact ties, not approximations.
- Therefore the finding's EXPECTED ("an error would cost nothing") is false twice over: the
  combination is a working, documented capability, and erroring at parse would delete it.

Caution placement (awk over the log): the runtime caution fires on BARE, on HELPFIX and on
HELPFIX_UNC — i.e. it fires even when the user has ALREADY written the i. prefix — and does
not fire on CTRLFIX (no fe()). Mildly noisy, not wrong: absorption genuinely is not happening
on that path, the terms are ordinary regressors.

## Additional documentation found, beyond what the skeptic cited
The behaviour is documented in THREE places, not one:
- `sthlp:1052` — the fe() entry, the BARE-token paragraph quoted above.
- `sthlp:268` — the `unconditionalse` entry: "...it cautions with aweights and with fe() or
  ccc_absorb(), which the regress path cannot absorb - specify those as controls() with your
  own i. prefixes."
- `sthlp:273` — the `customvce()` entry: same sentence.
Both option entries also cover the `ccc_absorb()` half of the same runtime caution, so the
adjacent lead (does ccc_absorb() under customvce get the same bare-token treatment, and is
that documented?) is documented too. Separately, `sthlp:532` shows the author tracking this
exact distinction elsewhere: under bootstrap(), "The check compares bare tokens, so fe(state)
is noticed where fe(i.state) is not."

## Where I differ from the skeptic (marginally)
Same verdict, stronger basis. The skeptic showed only that the two spellings differ (R2 0.96
vs 0.73), which does not establish that the documented remedy WORKS — a documented remedy
that failed to recover the absorbed answer would have been a live residual and would have
made this PARTIAL. It does work, exactly. I also add the ruling record (BS-1, documentation-
only, 2026-09-07), the proof it was BUILT (e26b3dd + md5 identity of pinned build vs HEAD),
and the two further help entries.
