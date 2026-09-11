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

Every item from the ruled slates has been built, ruled closed, or resolved by measurement;
PS-3, the last of those, was built on 2026-09-10 (see CLOSED below). The items below are
new, and all come from one source: the **documentation review of 2026-09-10** (a 12-agent
read of the thin regions of the ado, the Mata block and the Python harness). Each is
recorded as a comment at the site where it bites, so nothing here is invisible to the next
reader; they are listed here because each is a real deviation rather than a gap in the
prose. None was fixed inline — the review's own commits are comment-only apart from two
text corrections named in them.

**D-1. `permute()` + `cfxplotnose` silently leaves the plot cells analytic.**
`wooldid_permcellsx` skips any cfx matrix with fewer than six columns, and `cfxplotnose`
sets `cfxmatendpoint` to 1 unconditionally (`fastnose` has an explicit `permute()` carve-out;
`cfxplotnose` does not). So every cfx matrix is skipped, `ndone` stays 0, the drift caution
never prints, `e(perm_cfx_blank)` is never posted, and the user gets analytic plot cells
under a permutation header with no message at all. This is the "no silent no-op" policy's
own failure mode. Fix is either a carve-out matching `fastnose`'s or a caution on the skip.

**D-2. The drift caution's wording does not match what the code blanks.**
`wooldid_permcellsx`'s drift branch blanks columns 4-6 always and column 3 only when the
statistic is not the t-ratio. Column 2 keeps the ANALYTIC standard error, and under
`permtstat` column 3 keeps the analytic b/se. The caution says the cells' "inference columns
are blank (points only)", and the replay line repeats it. Either the wording should say
"no permutation inference" (which is what is true), or the blanking should extend to the se
and the `permtstat` t. Number-moving in the second form, so it is the author's call.

**D-3. Four `mata: M = ...` one-liners leave a global Mata `M` in the user's session.**
`wooldid_bootcellsx` (one), `wooldid_permwire`'s gteffects blanking (one) and
`wooldid_permcellsx` (two). Mata has no statement-local scope at the `mata:` prompt, so a
user's own matrix named `M` is silently overwritten, once per family per cfx matrix. Costs
wooldid's results nothing; damages interactive state with no error. Fix is a `wooldid_`
prefix or a `mata drop` tail — decide once for all four.

**D-4. `dots()` is not honoured on the physical-expansion leg of `wooldid_bootrun`.**
That leg prints one dot per replication and never the 50-mark running count, while the
count path, the re-identification leg and the jackknife loop all gate on `mod(r, dots)`.
Cosmetic, but a test that counts marks on this leg must expect `reps` of them.

**D-5. Inert code, now documented as inert rather than deleted.** Four sites, none harmful,
each of which reads as a bug or invites a wrong "repair": the two `if level != 95` interval
renames in `wooldid_permwire` (they can never fire — the worker already names the columns
from the run's own level); the pre-fill loop in `wooldid_bootjwr` (overwritten wholesale
five lines later); `matrix rownames `summarystats' = `rownamessum'` in the worker (expands
to an empty namelist, a measured silent no-op — the real labelling is inside
`wooldid_sumstat`); and the initialisation loop in `wooldid_sgdom` (locals that die at its
`end`). Deletion is safe for all four; it is not a documentation matter and was left alone.

**D-6. `wdptsmiss` is counted and never read.** `wooldid_permwire` increments it once per
coefficient blanked under `permtstat` for want of a usable observed analytic se, but the
caution states both causes in one sentence without being able to say which applies. The two
are different fixes for the user (raise `permreps()` versus "this element has no identified
standard error to studentize by"), so the counter is worth spending on a caution split.

**D-7. Per-family simultaneous-band drops are silent.** `wooldid_esbandv` exits without a
message when a family has no stored band, while `wooldid_esband`'s caution fires only when
NO family got one. Under `simultaneousbands` a user can get a dashed band on the ATT figure,
none on the AME_it figure, and be told nothing.

**D-8. The two band kernels take opposite policies on a degenerate horizon.**
`wooldid_esbanda` drops the horizon and bands the rest; `wooldid_bootband` refuses the whole
family. So an 8-horizon analytic band can come with no bootstrap band at all when horizon 9
is structurally unidentified. Both are defensible; the asymmetry is undecided rather than
ruled.

**D-9. `wooldid_esbanda` consumes the caller's RNG stream, seeded once for the family loop.**
A family's critical value therefore depends on how many families were banded before it, so
adding `subgroups()` or an estimand moves an unrelated family's band in its low-order digits
(simulation noise at nd = 10000). Reproducible within a run; not reproducible across runs
that differ in family count.

**D-10. `wooldid_permsaveasg` has two silent drops** — a conformability return that writes no
assignment block at all (indistinguishable from a run that never asked for one), and a skip
of any cluster whose `strtoname()`d column name is already taken, which leaves the saved file
one cluster short with nothing marking the gap.

**Harness items** (same review, `test/` files):

**D-11. `mkrun.py`'s saveplots pin matches one spelling.** The option's legal minimum
abbreviation is `save`, which resolves uniquely, so `save(svp_x)` is valid Stata the regex
misses — re-opening the hole SI-21 closed. `permsave()` and `bootsave()` are outside the pin
altogether.

**D-12. `test/compare.py`'s 1e-12 absolute floor is live at every nonzero rtol**, including
the 1e-9 default the README's own invocation uses; the B21 fix exempted `--rtol 0` only.
`compare_oldnew.py` deliberately has no absolute floor, so the two comparers disagree on
near-zero quantities.

**D-13. `coverage_check.py` can overstate coverage without failing.** An UNRECOGNIZED option
token never sets `fail`, and the abbreviation matcher resolves an ambiguity to whichever
option comes first in the syntax line where Stata would refuse it (measured: no ambiguous
pair exists today, so this is right by luck of naming). Separately the weight classifier
files `fweight` — which the syntax line accepts — as unweighted.

**D-14. `check_validate_sync.py`'s assignment regex counts a READ as an assignment.** Any
`local X = `nsg'` line registers as assigning `nsg` and suppresses the finding, in the script
whose whole job is to catch stranded syntax-parsed values. Benign today (only `item` is
affected, via `tempvar `item'`), and the mechanism is ordinary.

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
- **F033 propagation to the studentized CELL p-value** (found by the documentation review,
  2026-09-10; recorded as an OPEN QUESTION at the site in `wooldid_bootcellwire`). The cell
  branch ranks recentred draws — `tr = (b_r - b)/se_r` against `to = b/se_obs` — with the
  plain `2*min(m0,m1)/n` rule, while `wooldid_boottp` and `wooldid_bootpwt` use the add-one
  form `2*(min(m0,m1)+1)/(n+1)` that F033 mandates for exactly that construction. As written,
  under `boottstat` a cell p and the coefficient p for the same quantity differ by roughly
  `2/(n+1)` — about 0.02 at `reps(99)`, the size defect F033 was measured to remove — and the
  cell p can still print as exactly 0 where the coefficient p cannot. The non-studentized
  cell branch correctly has no add-one (its draws are not recentred). Moving reported numbers,
  so held.

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
