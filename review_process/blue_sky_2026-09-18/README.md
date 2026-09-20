# Blue Sky Review of wooldid — 2026-09-16 to 2026-09-18

**Status line, dated: written 2026-09-18; fix proposals and their adversarial review added 2026-09-19.**
**Nothing in the package is fixed except the three edits recorded in `shipped_edits.md`.**
This is an identification-and-investigation record, per the author's instruction. It will go stale
the moment anything is patched; what would make it false is any commit touching the sites it cites.

Pinned build: `c1d8979`. The existing testing-and-assessment suite was **not used** and was
off-limits by design — it is green, and the author asked for probing outside its boundary.
Mechanically verified over the final corpus: of **1,168** agent do-files, **exactly one** names
anything under `test/`, and it does not run the suite — `fineclus-cellwt-refute/n05_stalecomment.do`
greps two suite *prose* files for a stale comment while checking documentation drift. No probe
executed `runbattery`, a `test_*.do`, a baseline comparison, or a suite dataset. The repo working
tree was never modified during the review phase.

## How this was produced

| wave | agents | what |
|---|---|---|
| Recon | 10 | 11 attack catalogs: option inventory, ~110 metamorphic invariants with measured tolerance classes, e() contract, message inventory, interaction matrix, help-claim extract, data generators, architecture, and a QA'd "already ruled" suppression list |
| Structure/UX | 16 | 8 static reviewers + 8 cross-critics on the program as a whole |
| Probing | 160 | 27 method-assigned adversarial cells + 12 verification-only cells, each pipelined through an independent refuter and a third lens on high-stakes items |
| Coordinator | — | 5 independent Stata verifications of the highest-stakes claims |

Each probe cell was assigned a **method**, not just a target, to force different tacks:
first-principles reimplementation oracle, exact enumeration oracle, metamorphic invariance,
internal algebraic identity, Monte Carlo size/coverage, numerical perturbation, state/leakage
differencing, realistic messy data, help-conformance, and randomized option fuzzing.
**83 of 149 probe findings were self-found rather than from seeded leads.**

Verification outcome: 134 verdicts — 75 CONFIRMED, 55 PARTIAL, 4 REFUTED — plus **74 findings the
refuters discovered themselves**. 37 of 39 refuter cells completed.

**Read the high tier with this in mind — corrected 2026-09-19.** The paragraph that stood here
described 45 critical/high findings of which 18 were single-source. That is no longer the state, and
the count it gave was itself produced by a defective join (see below). After three further
adjudication rounds **every one of the 223 findings carries at least two independent reviewer
opinions** — 118 with two, 94 with three, 11 with five. Severity is now the **median** of all
opinions rather than the most severe, which is why the high tier reads 48 rather than 66.

Final tier: **4 critical, 48 high, 87 medium, 72 low, 12 reduced to none.**

The reviewers are still not interchangeable, and `mechanisms.md` and `all_findings.md` name which
kind each finding got: 118 findings were checked by a prober and an independent refuter that wrote
its own repro, 16 of those also got a third lens, and **89 were checked only in the later
adjudication rounds**, because their refuter or lens cells were lost to the weekly usage limit.
"3 reviewers" does not always mean prober+refuter+lens.

**A defect in this report's own pipeline, and what it cost.** Refuters independently chose generic
finding ids (`SK-01`, `SKEP-01`, `SK-1`, `SK-2`), and the join that attached verdicts to findings was
keyed on that id alone, so distinct findings collapsed into one another and some verdicts landed on
the wrong rows. It was caught when a critical finding still showed one reviewer after a round that
should have given it three. `mechanisms.md`, `all_findings.md` and `refuted.md` were regenerated on
2026-09-19 from a dataset keyed on unique ids; the mechanism groupings were re-derived at the same
time, which is why `SG` no longer carries `permute()`'s name. The four criticals survived the
correction unchanged and all four are confirmed by multiple reviewers.

The four criticals and the headline are the best-evidenced items here; the criticals were
additionally reproduced by the coordinator on independent designs.

## What is NOT wrong — read this before spending effort

- **The core estimand machinery is clean.** A from-scratch numpy oracle reproduces `att`, `att_it`,
  `att_i`, `att_itime`, both pretreatment rows and every realized event-study horizon to
  **<= 2.7e-13 across 14 option regimes on 5 designs** (unbalanced, `aw`, `pw`, `esfixedbaseperiod`,
  `customweightmultiplier`). The defects in this report are in the inference and peripheral layers,
  not in the estimator.
- **Recon's top lead was refuted.** `customweightmultiplier()` entering all four weight numerators
  and no denominator is coherent and authorised by sthlp:1107. Two oracles were built — one from the
  code's formula, one from the help's literal formulas — and only the help-literal one disagreed.
  What survives is a LOW documentation gap (sthlp:220-236 omits the `cwm` factor while still calling
  the result "the simple, unweighted mean").
- **The caution corpus is not too big.** 84 caution blocks over 17,021 lines; MEASURED, a bare
  `att` run prints **one** caution and the common `cluster(i) jointtests` run prints **two**. Of 141
  Caution/Note sites read in full, exactly **five** carry textbook-econometrics content.
- **91.6% of the help is on-topic reference.** The (a)/(b)/(c) material the two-track split targets is
  8.4% of the file. Deleting all of it would leave the help essentially as long. That file needs
  **structure**, not deletion: `permute()` is 2,481 words in a single paragraph, `bootstrap()` 1,668
  with no internal structure.
- **Twelve findings were reduced to severity `none`** by review, and 25 REFUTED verdicts were
  returned across the rounds. They are listed in `refuted.md` — the process cutting its own false
  positives. The later rounds killed single-source findings at **21%**, against ~3% for findings that
  had already survived a prober-refuter pair: a measure of how much the pairing was worth.

## Corrections to interim statements made during the review

Stated once, not re-litigated:
1. sthlp:1551's worked example is the **lattice** arm, not a count outcome and not `poisson`.
2. sthlp:1097's "impvalues ... remains valid" sits inside the contreat cell-containment paragraph,
   i.e. the **safe** regime. That sentence is correct where it stands; the `impvalues` realism case
   rests on sthlp:67/255/257/375 instead.
3. "The bootstrap `e(V)` is rank-deficient at every `reps()`" is **reframed, not confirmed** — see
   `unadjudicated.md`.

## Where to read next

- `headline.md` — the one result that should drive what happens next
- `critical.md` — the four critical findings, with verdict chains and repro paths
- `mechanisms.md` — the high tier grouped by **fix site**, not by symptom
- `refuted.md` — what the verification stage killed
- `unadjudicated.md` — claims I could not settle, stated as open questions
- `all_findings.md` — the full 223-record appendix by mechanism
- `shipped_edits.md` — **the only package changes this review made** (two help edits and one
  `impvalues` caution), including the predicate defect the adversarial round found in the caution
  *after* it was written, and the measurement that fixed it
- `fix_proposal_status.md` — 56 fix proposals after two adversarial rounds: the author's ruling
  table for the 13 that move reported numbers, the 6 not to implement, and the per-proposal
  corrections that make up each stable specification
- `fixprop/UZ_FINAL.md`, `fixprop/FIX2_FINAL.md` — the two finalized fix specifications
- `fixprop/reviewer_verdicts.json` — the raw reviewer returns behind that status file
- `closes_held_item_5.md` — the one diagnosis run held-for-ruling item 5 (2026-09-13) asked for, and
  its answer: neither documented containment nor a lattice-construction defect, but exact fit of
  saturated lattice pieces
- `fixprop/proposals_*.md` — fix proposals by mechanism, under adversarial review
- `repro/` — the six coordinator runs, and `repro/agents/` with 1,168 agent do-files by cell
  (8 MB; delete it if you do not want it in the repo — nothing else depends on it)

## findings_raw.json

The full, untruncated record of all 223 findings - claim, expected, observed, refuter note
and lens detail - which the published reports abbreviate. **Join on `uid`, not on `id`.**
Nine `id` values were assigned independently by different agents and name unrelated
findings on different surfaces; `SK-01`, for instance, covers a subgroup-degeneracy guard,
a semielasticity averaging defect and a cell-collapse read. Keying on `id` silently
conflates them, which is the defect the published reports were regenerated to remove.
`uid` is the `id` where it is unique and `id__1`, `id__2`, ... where it is not, so prose
references to an unambiguous id still resolve.
