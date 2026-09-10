# Open defects ledger

Running list of defects found while working the ruled slates, held for a fix pass at the
end rather than fixed inline. Opened 2026-09-10 at the author's instruction ("keep the
defects found in a separate to do file which we will address at the end").

Rules for this file. A defect goes in OPEN when it is a genuine bug or a real deviation
that (a) is not what the item under work was about, or (b) would move a reported number
and so wants its own verified commit. Anything that is a judgment call for the author
goes to `held_for_ruling_2026-09-09.md` instead, not here. Anything fixed gets moved to
CLOSED with its commit, so this file is also the record of what the sweeps bought.

---

## OPEN — to address at the end

### O-1. PS-46: measured 2026-09-10, and the answer needs a ruling before code moves
`permute()`'s fast-mode default posts no `e(V)`. The ruling is conditional: *"First
establish whether the permutation fast-mode default actually breaks esttab-class table
tools. If it does, post an ARTIFICIAL covariance that lets those tools run and that is
transparently false to anyone who inspects it (an all-zero matrix being the example the
author gave), rather than leaving nothing posted. If it does not break them, leave as
is."* **MEASURED: yes, it breaks them, silently.** `esttab` prints an empty standard-error row
under a footer asserting standard errors are present, and writes that into `.tex`/`.csv`.
Report: `review_process/research_2026-09-10/permute_fastmode_esttab.md`. The all-zero
covariance the ruling names fails the ruling's own "transparently false" test on the
surfaces it was meant to protect, and a better route (post `e(se)`/`e(p)`, which estout
reads ahead of `e(V)`) was measured working. **The choice is now the author's** — the five
options are tabulated in `held_for_ruling_2026-09-09.md`. No code change here until then.

The one piece that moves regardless, because it is a false statement about current
behaviour rather than a judgment call: the help says `esttab` "will refuse" under fast
mode. It does not. Corrected in the P8 sweep.

### O-2. `contreatwithin` silently disables the whole collinearity-screen family, and one of them reports a structural zero as an estimate
`contreatwithin` with `contreat()` rewrites the working `safety` macro to `off` (ado:4595).
`wooldid_ccscreen` is gated only on that rewritten value (ado:1599), so a user who typed **no**
`safety()` at all loses the entire family — including the SI-A4/R6 ame-void screen at
ado:10703-10722, which exists to report an unidentified dose slope as **missing** rather than as
the regression's structural zero. Measured on `dta_cont.dta` with a cohort-constant dose: without
`contreatwithin` the table prints `ame | . . . . .` under the named caution; **with**
`contreatwithin` the same run prints `ame | 0 . . . .` and says nothing. rc 0 both times.
**Severity: this is a reported number that is wrong, printed in silence.** Narrow fix: gate
ado:10703-10722 (and the ado:1599 family gate) on `safetyinput` rather than on `safety` — the
same one-token change P6 ruled for the maxvar preflight at ado:1247, at a site the ruling did not
name. It moves a printed cell from `0` to `.`, so it is the author's call rather than a silent
correction, and it wants its own verified commit.

### O-3. `e(df_r_analytic)` is kept only on the `bootdf(normal)` branch
`wooldid_bootwire` (ado:9163-9179) reads the analytic residual df into a local and then posts the
hidden twin **only** in the `else` branch — the `bootdf(normal)` case. Under `bootdf(clusters)`,
which is the default whenever `cluster()` is given, the analytic value is dropped. Measured:
`bootstrap(standard) cluster(i)` posts `e(V_analytic)` and `e(V_es_analytic_att)` but no
`e(df_r_analytic)`; `bootstrap(bayes) cluster(i) reps(20)` gives `.` for it with `e(df_r)` = 79.
Three help sites promise it unconditionally, and the phrase "alongside the analytic covariances"
is exactly the promise that fails — the covariances are kept unconditionally, the df is not. Fix:
hoist the post out of the `else`. Pure addition to `e()`; no printed number moves.

### O-4. The `unconditionalse` prerequisite screen accepts any `customvce()`, robust or not
ado:4442 tests only that `customvce()` is non-empty. The screen's own comment says it exists to
catch a bare r(322) early, and `wooldid y i t ttre, att unconditionalse customvce(ols)` walks
straight through it into that r(322) — `option vce(unconditional) is allowed only after
estimation with robust, semirobust, or linearized standard errors`. Fix: test the `customvce()`
text for a robust/cluster/semirobust/linearized form rather than for non-emptiness.

### O-5. `makecfxplotsbysg` is a silent no-op in two situations
Without `subgroups()` the subgroup list is never populated, so its loop runs zero times; without
`makecfxplots` the whole block is unreachable. Both return rc 0 and produce nothing. This is
against the package's own precedent: `makecfxplots` without `contreat()` (ado:5012-5017) and
`lattice()` without `contreat()` (ado:4581-4586) were both turned from silent no-ops into errors,
and P6's ruled policy is that no accepted option is a silent no-op.

### O-6. A caution points at another caution that was never printed
The all-zero-covariance disclosure at ado:1453 ends by referring the reader to "the
cell-containment caution above", and prints unconditionally — but the caution it names is
produced only when `contreat()` is given (screen gated at ado:1420). On a binary-outcome run the
reader is sent to a message that does not exist. Fix: make the pointer conditional, or inline the
one clause it is pointing at.

### O-7. UNRESOLVED — does a studentized subgroup contrast disclose its own degeneracy?
Under `bootstrap() boottstat cluster(i) subgroups()` the analytic pairwise-degeneracy screen never
runs (ado:4256 gates it off the engine paths), and the bootstrap-side guard
(`wooldid_bootsgw`, ado:13864) is a zero/missing test, not a rank test — the comment at
ado:9905-9906 says as much. On a design with pairwise residue se .0018 and rho .9995 nothing
should fire, which would leave the studentized contrast built silently on the residue. **Not
confirmed by a run.** What would settle it: that exact call on that fixture, checking whether any
caution prints and what the contrast's p-value is.

---

## HELD FOR THE AUTHOR (not defects — judgment calls)

These live in full in `held_for_ruling_2026-09-09.md`; listed here only so one file shows
everything outstanding.

- **PS-2** — ruled but not built.
- **PS-62 sentinel diagonal** — the ruling said zero the sentinel row and column; the
  off-diagonals are zeroed and the `-1e100` diagonal is KEPT, because zeroing it would
  print a fabricated `0` standard error where an honest `.` prints today. Two-line
  reversal recipe recorded.
- **BS-15 degeneracy cut** — moved from the spec's `1-1e-3` to `1-1e-2` on measured
  evidence.
- **PS-45 bootstrap band studentization** — the band scales by the coefficient-draw
  standard deviation even under `boottstat`, where every other reported quantity uses each
  replication's own standard error. Cost of the alternative is written up.

---

## CLOSED

### C-1..C-4. The four defects the batch-5 shared assessment surfaced (8280a2e, 2026-09-10)
1. `cleanmatrices` runs `clear mata` inside the worker, dropping wooldid's own Mata
   kernels; PS-62 had put a Mata call on the mainline analytic path, so a plain
   `att cleanmatrices` fit died `r(3499)`. Guarded. **No cert suite caught this** — only
   the 501-config dual-build comparison did. `test_stored` S9 now pins it.
2. `wooldid_bootwire` relabelled `e(wdid_vcvstate)` as `analytic-diagonal` on the
   repost-failure fallback, even when the matrix left standing was the spliced one.
3. The PS-62 splice ran inside every bootstrap/permute replication, where nothing reads
   `e(V)`'s off-diagonals. Gated on the replication marker.
4. `e(wdid_vblocks)`/`e(wdid_vblockn)` describe the hidden `e(V_analytic)` under
   `bootstrap()`, not the reported `e(V)`. The help now says so.
