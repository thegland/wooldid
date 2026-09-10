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

### BUILT 2026-09-10 on the author's rulings

- **O-1 / PS-46 — ruled: post `e(p)` and the interval, nothing else.** `wooldid_permwire`
  now posts `e(p)`, `e(ci_l)` and `e(ci_u)` under the fast-mode default, so `esttab` and
  `estout` carry this run's real permutation p-values and shift intervals. No covariance and
  no `e(se)`: `test`, `lincom`, `_se[]` and `estat vce` go on refusing, and nothing is
  invented, because a permutation run has no standard error to report.
- **O-3 — the analytic residual df is preserved on every `bootdf()` branch**, not only the
  one that blanks `e(df_r)`.
- **O-6 — the all-zero-covariance caution** no longer points at a caution that was never
  printed; off the `contreat()` path it names the alternatives directly.
- **PS-45 — ruled: percentile-t.** Under `boottstat` the simultaneous band now uses each
  replication's own standard error in the denominator and the observed analytic ones in the
  endpoints, so it sits on the same scale as the pointwise intervals printed beside it.
  `e(esband_engine)` gains a third value, `bootstrap-studentized`.

### CLOSED BY RULING 2026-09-10 — not defects, not to be reopened

- **O-2 `contreatwithin` and the collinearity screens.** Put to the author with the two
  console tables side by side (`ame | .` with the caution, versus `ame | 0` without it).
  **Ruled: "Completely ignore this. It's fine. It's basically already documented. No
  reasonable observer will look at a pile of all 0s with no SEs and mistake them as real
  results."** The safety inventory now describes the relaxation as the named exception it is
  rather than as a gap. Do not raise this again.
- **O-4 the `unconditionalse` prerequisite screen** (accepts any `customvce()` text).
  **Ruled unimportant; not to be considered again.**
- **O-5 `makecfxplotsbysg` as a silent no-op** without `subgroups()` or `makecfxplots`.
  **Ruled unimportant; not to be considered again.**

### O-7 — RESOLVED 2026-09-10 by running it: the disclosure already exists

The worry was that `bootstrap() boottstat cluster(i) subgroups()` might build a studentized
contrast silently on a near-degenerate analytic pairwise standard error, because the analytic
screen does not run on that path and the bootstrap-side guard is a zero/missing test rather
than a rank test. Run: rc 0, contrasts studentized (p = .118 and .083 on the fixture), and a
caution DOES print, on exactly this point:

> Caution: under clustering at or above the i level the analytic pairwise-contrast standard
> errors can be degenerate (the same limitation behind the refusal of jointtests with
> analytic clustered SEs), and boottstat studentizes each contrast by that standard error in
> every replication, so the studentized contrasts and the across-subgroup equality test can
> be absent or numerically unstable here; the non-studentized bootstrap (bootstrap() without
> boottstat) forms the contrasts from the coefficient draws instead.

It names the mechanism, says the studentization consumes the degenerate quantity, and names
the remedy. Nothing to fix. Closed.

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
