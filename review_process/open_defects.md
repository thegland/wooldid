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

## RULINGS ON THE DOCUMENTATION REVIEW'S ITEMS (author, 2026-09-11)

All fifteen were put to the author as multiple choice and all fifteen were ruled. Build
order below is by risk, not by item number: the harness first (so everything after it is
judged by the stricter gates), then the changes that move no number, then the three that do.

| id | ruling |
|---|---|
| F033 cells | **Apply the add-one.** `2*(min+1)/(n+1)` in the studentized cell branch, matching `wooldid_boottp` and `wooldid_bootpwt`. One convention for every recentred rank-based p. Moves reported cell p-values; needs a re-bless. |
| D-1 | **Carve-out like `fastnose`.** Keep the 6-column cfx shape under `permute()` even with `cfxplotnose`, so the cells get real permutation inference. |
| D-2 | **Fix the wording.** "no permutation inference", naming which columns survive; the analytic SE stays, as on every other permutation surface. |
| D-8 | **Keep the refusal, say why in the caution.** The bootstrap critical value is a maximum over the family taken replication by replication; the caution must name the horizon that caused the refusal. |
| D-7 | **Caution naming the families** that got no band, replacing the all-or-nothing one. |
| D-10 | **Caution on both paths** in `wooldid_permsaveasg` (no assignment block written / N cluster columns omitted). |
| D-6 | **Split the caution** so `wdptsmiss` is spent: too few completed replications versus no usable observed analytic SE. |
| D-3 | **Rename to a `wooldid_` name** at all four `mata: M = ...` sites. |
| D-5 | **Delete all four inert sites.** All behaviour-neutral; the explanatory comments go with them. |
| D-4 | **Honour `dots()`** on the physical-expansion leg, so the option means one thing everywhere. |
| D-9 | **Reseed per family** in the band loop, so a family's band depends only on that family. Moves band endpoints; needs a re-bless. |
| D-11 | **Widen the regex and add both savers** to `mkrun.py`'s file pin. |
| D-12 | **Remove `compare.py`'s absolute floor at all rtols**, matching `compare_oldnew.py`. |
| D-13 | **Make UNRECOGNIZED fatal and add the `fw` bucket** in `coverage_check.py`. |
| D-14 | **Tighten `check_validate_sync.py`** to a real assignment. |

---

## OPEN — to address at the end

Nothing. Every item recorded in this file has been built, ruled closed, or resolved by
measurement. The fifteen items the documentation review of 2026-09-10 turned up were ruled
on 2026-09-11 and all fifteen are built - see CLOSED below.

---

## HELD FOR THE AUTHOR (not defects — judgment calls)

Nothing. This list was stale: all five items it carried were ruled on 2026-09-10 and are
recorded as built or confirmed in the table at the head of
`held_for_ruling_2026-09-09.md` — PS-2 (name the exporters only) built, PS-62's sentinel
diagonal confirmed as built, BS-15's cut confirmed at 1 - 1e-2, PS-45's band built as a
percentile-t under boottstat, PS-46 built. The F033 cell p-value entry added on 2026-09-10
was ruled on 2026-09-11 and built the same day (see CLOSED).

---

## CLOSED

### The documentation review's fifteen items — ALL BUILT 2026-09-11
Ruled by the author on 2026-09-11 (multiple choice, all fifteen answered; the ruling table is
above). Built in three risk-ordered commits, each gated.

**`ed95674` — the harness (D-11..D-14), first, so everything after it is judged more strictly.**
- D-12 `test/compare.py`: the 1e-12 absolute floor is gone at every tolerance, matching
  `compare_oldnew.py`. The selftest now pins the property with the two cases that used to
  pass: `close(0.0, 1e-300)` and `close(1e-13, 5e-13)`.
- D-13 `coverage_check.py`: an UNRECOGNIZED option token FAILS the run; `fweight` has its own
  bucket.
- D-14 `check_validate_sync.py`: the assignment regex requires the name to BE the thing
  assigned. Tightening it immediately surfaced `item`, which the worker sets fourteen times
  as a `foreach` loop variable and the old pattern had been suppressing for the wrong reason;
  `foreach`/`forvalues` are recognized now, so the script is OK on its own terms rather than
  by accident.
- D-11 `mkrun.py`: the file pin covers every file-writing option at every legal abbreviation,
  and `dumpres` takes a list of prefixes. `save(...)` used to slip it silently; `permsave()`
  and `bootsave()` were outside it altogether.

**`9bb24c6` — the seven that move no number (D-2..D-8, D-10).**
- D-7/D-8: `wooldid_bootband` reports WHY it refused, naming the horizon; `wooldid_esband`
  collects the declined families and names them. The refuse-the-whole-family policy stays, as
  ruled, and the caution now carries its reason.
- D-2: the drift caution says what is true — no permutation inference, p and interval blank,
  the se column analytic, the t blank unless `permtstat`. The replay line matches.
- D-6: the permutation blanking caution is split, so `wdptsmiss` is finally read: too few
  completed replications (raise `permreps()`) versus no usable observed analytic SE (a
  property of the specification).
- D-10: both silent exits in `wooldid_permsaveasg` now report.
- D-3: the four one-liner Mata temporaries are `wooldid_Mtmp`.
- D-4: `dots()` is honoured on the physical-expansion leg, using the documented cadence.
- D-5: the four inert sites are deleted, each leaving a note saying what went and why it
  could not fire.

**`82934b3` — the three that move numbers, each with before/after evidence.**
- **F033 on the studentized CELL p-value.** Now `2*(min+1)/(n+1)`, one convention for every
  recentred rank-based p. The effect was larger than the estimate in the finding: on the cert
  fixture at `reps(40)`, the OLD build produced **22 of 50 cell p-values at exactly 0** with a
  minimum of 0; the new build produces none at zero and none below `2/(completed+1)`.
- **D-1, `permute()` + `cfxplotnose`.** `cfxplotnose` now takes the same `permute()` carve-out
  `fastnose` has. Old build: the first cfx matrix had `colsof = 1`, no `e(permdist_*)`, and a
  missing cell p. New: `colsof = 6`, `permdist_` posted, a real permutation p.
- **D-9, the band's seed.** Derived per family from the family's NAME rather than set once
  before the loop. Same fixture and spec, `e(esband_crit_att_it)` with `att` banded before it
  versus `att_it` drawing first: before `2.474072413106145` vs `2.492870068598342`; after
  `2.47972450315949` in both.

**Pinned.** None of the 504 dual-build configs reaches these three paths, so the changes would
have been invisible to the suite. Three cert legs were added and each was verified to FAIL on
the pre-change build: `test_esband.do` B10 (the band's critical value is keyed to the family,
not the draw order), `test_permute3.do` Q33 (`cfxplotnose` no longer truncates under
`permute()`), `test_bootstrap_tstat.do` T12 (no studentized cell p is zero or below the
add-one floor).


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
