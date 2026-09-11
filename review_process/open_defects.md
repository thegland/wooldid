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

Nothing. Every item recorded in this file has been built, ruled closed, or resolved by
measurement. PS-3, the last one, was built on 2026-09-10 (see CLOSED below).

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

### PS-3 — built 2026-09-10, after a wrong diagnosis
**Ruling.** *"Per-column formats. Keep the table structure, give each column an explicit
format."* **Author ruling 2026-09-10 on the four options put to him: the conventional
formats, as ruled** - accepting that a p of 1.84e-21 prints as 0.000, which is what every
comparator does. The magnitude is not lost from the data, only from the printed cell;
`e(wdidmainresults)` still carries the full value.

**Built.** `wooldid_mlfmt` builds matlist's `cspec()` and `rspec()` from each matrix's own
names, so no call site hard-codes a column count: statistics two decimals, p-values three,
`relyr` an integer, everything else a general format so no estimate loses precision. The
p-value rule is EXACT-match, because the summary-statistics table has `p25` and `p75`
columns that are percentiles. The rules are placed where the default put them - none above
the header, one under it, one at each equation-group boundary, none at the foot - so the
structure is unchanged and only the number formatting moves.

**The earlier failure was mine, not matlist's.** The first attempt was abandoned on a probe
that passed a `nooutput` option matlist does not have, so every row-spec length errored
regardless, and this file recorded "matlist's contract fights these tables" as a finding. It
does not. `rspec()` takes one character per horizontal rule and the count is rows + 2,
independent of equation groups - measured on all three tables. **A single failing probe is
not a property of the tool.**

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
