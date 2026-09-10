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

### O-8. PS-3 (per-column table formats) — attempted 2026-09-10, reverted, still open
The other eleven items the audit found are closed (see the BUILT section). This one is not.

**Ruling.** *"Per-column formats. Keep the table structure, give each column an explicit
format."*

**What was tried.** A `wooldid_mlfmt` helper building a `matlist cspec()` from each matrix's
own column names - general format for estimates and interval bounds, two decimals for
statistics, three for p-values, integer for `relyr` - so no site hard-codes a column count.

**Why it was reverted.** `matlist` refuses `cspec()` without a matching `rspec()`, and these
tables carry row-equation groups whose separator rows `matlist` inserts itself, so the row
count is not `rowsof()` and not any simple function of it. Three row-spec forms
(`rowsof+1`, `3*rowsof+8`, and a generous fixed pad) all failed with
`too few specifications in rspec()` - on the summary-statistics table, before the main one
was reached. The tables were left on their working renderer and re-verified rendering.

**What would settle it.** Either determine `matlist`'s separator-row rule for
equation-grouped matrices empirically (count the rules it emits for a known matrix and derive
the formula), or take the sweep's other named route and build the tables through
`_coef_table` instead of `matlist`, which formats by convention and needs no spec at all.
The second is the larger change and the one that would actually make the output look like
every comparator's, which is what the ruling is after.

**Priority: cosmetic.** No number is wrong; a p prints as `1.12e-35` where `regress` would
print `0.000`.

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

### The eleven audit items closed 2026-09-10
From `review_process/research_2026-09-10/sweep_audit.md`. PS-3 is the twelfth and is above.

- **PS-20** the Ibragimov-Muller banner: the header now names IM inference and its reference
  distribution instead of the analytic cluster specification the run did not use.
- **PS-64** the joint pre-trend caution now fires under `pretrendpower()` as well as
  `jointtests`, and names whichever entry point asked for it - "one statistic, one rule, in
  both entry points", as the ruling says.
- **CS-98** resolved by NAMING rather than by moving: the two screens that run before the
  `contreatwithin` resolution now read `safetyuser`, snapshotted once at the top of the
  validator, so no reader has to know a line number to know which value a gate means. Hoisting
  the resolution instead would have extended the relaxation into a screen it has no
  relationship with, and those two screens terminate the run, so they were never in CS-98's
  hazard class - a message overturned by a later rewrite. No behaviour moves.
- **CS-60** the subgroup joint p is a scalar at all THREE posting sites (the audit found two;
  the permutation engine had the same pair as the bootstrap one).
- **CS-76/CS-115** the bootdisplay() entry no longer claims every family is computed on every
  run - a claim it contradicted five sentences later.
- **CS-118** `bootdisplay(bca)` added to the `e(boot_stat)` value list; `e(boot_display)`'s
  legal values enumerated.
- **CS-120** `e(wdid_fcncell)` and `e(wdid_fcskip)` documented, and the Replay paragraph no
  longer disowns the two cell-collapse names the `fastcollapse` entry tells users to read.
- **P8(ii)** the Standard Errors synopsis split into labelled sub-runs (analytic /
  shared-by-both-engines / permutation-only / bootstrap-only); `ties` states that it requires
  `bootstrap()`; `permside` and `permrawp` carry their abbreviation markers; `regtol` shows
  its argument.
- **P7 accident 3** both inline compatibility lists reshaped: reasons on every item, and the
  two refusals `permute()`'s list omitted (multi-way `cluster()`, `customvce()`) added.
- **CS-25/40** the control-group refusal added to the inventory's class rows, and the "gates
  fifty-eight sites" claim replaced with what the section actually delivers.
- **SI-A1** (author ruling 2026-09-10, a note and only a note): the `cluster()` entry states
  that clustering below the level of `i` is the intended practice, what goes wrong at or above
  it (pointwise standard errors understated to roughly 0.66-0.75, joint tests oversized), and
  why there is no plug-in correction.


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
