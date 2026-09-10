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
