# PERMSEM-01 — LENS 0 (independent reproduction). Running notes.

## Assignment
Third, independent lens on PERMSEM-01. Prober said [low/certain]; skeptic said PARTIAL/low.
My job: build the case a third time from the claim alone and say whether the world behaves
as claimed.

## Source read (pinned build, done before any Stata)
- `wooldid_permqf` (ado:16193) and `wooldid_permadv` (ado:16386): both open with
  `Wr = J(rows(A),1,.)` and `st_matrix(outW, Wr)`, i.e. the W vector is allocated at the
  ACCUMULATOR length and failed/incomplete rows stay missing. The posted p uses
  `ok = rowmissing(A) :== 0`, `n = sum(ok)` and `nex = sum(select(Wr, Wr :< .) :>= Wobs)`
  — so the PRINTED p correctly ignores the missing rows.
- `wooldid_bootdist` (ado:14644) compacts: `ok = rowsum(A :< .) :> 0; select(A, ok)`.
  Called at ado:8191/8195 to build `e(permdist_b)`.
- `wooldid_bootsavemk` (ado:14871) writes `st_addobs(R+1)`, `_rep = 0::R` off the SAME
  accumulator `A`. permsave calls it at ado:8859 with the FIVE-argument form, so the
  permsave file carries `_rep`, the coefficient columns (and `se_*` under permtstat) and
  the `onset_*` block — and **no W_r columns**. bootsave (ado:6546) uses the 11-arg form
  and does write `jw`/`aw` columns under boottstat. That asymmetry is real.
- The `jn` local (`outn`, the number of replications behind a joint/adversarial test) is
  computed inside permqf/permadv and **never ereturned** — grep finds only `jp` posted.
  `e(perm_n)` is per-COEFFICIENT (ado:8146), not per-test.
- `e(perm_addone)` (ado:8204) is built from the per-coefficient kernel only (ado:16147) —
  it does not cover the joint/adversarial tests.

## Help read
- sthlp:1348-1349 (permutation stored results): lists `e(permdist_b)` and calls
  `e(perm_jointW)`/`e(perm_sgW)`/`e(perm_advW)` "the replication statistics behind the
  joint, equality and adversarial p-values ... the p is reconstructable from the pair".
  No statement anywhere that the W matrices are at the REQUESTED length with missing rows
  while `e(permdist_b)` beside them is compacted.
- sthlp:495 (permsave()): DOES state "a failed replication keeps its row, all missing" for
  the FILE, and closes "Every reported permutation quantity is otherwise reconstructable
  from the file." So the prober's "the help ... says nothing for the permutation side" is
  too strong: the file convention IS stated; what is missing is the e()-side bridge.
- sthlp:1340 (bootstrap block) describes `e(boot_jointW)` with no missing-row statement
  either — the gap is engine-agnostic (source-read only; not measured here).

## My probe (independent mechanism)
`L0_permrows.do` — ROLLING ADMINISTRATIVE FRAME, deliberately different from both earlier
repros (prober: permprobs with a zero-treated draw; skeptic: early leavers/attrition):
30 districts 2010-2019; 12 report the whole window, 10 enter the frame in 2015, 8 enter in
2017; three adopt in 2017, all three on full-window districts in the OBSERVED assignment.
A shuffle that lands every onset on a 2015 entrant leaves relyr -3 with no contributing
unit → the replication COMPLETES but its -3 column is missing (a PARTIAL row, which
neither earlier repro produced); a shuffle landing every onset on a 2017 entrant leaves no
pre-period at all → outright failure.
esprelength(3) espostlength(2) jointtests adversarial esfixedbaseperiod
permute(shuffle) permreps(120) rseed(4041) permsave(L0_save).

## Status
- [x] source + help read
- [x] Stata run A (`L0_permrows.do` / `.log`) — data + shapes measured; its mata block died
      r(3201) on a non-conformable `select(dd, dd :< .)` before the p-value arithmetic
- [ ] Stata run B (`L0b_permrows.do`) — same data, same seeds, arithmetic fixed
- [ ] verdict

## Results — run A (measured)
```
L0_NOBS=194   (30 districts, unbalanced by frame entry)
L0_reps=120 completed=120 failed=0          <-- ZERO outright failures
L0_DIM e(permdist_b) = 120 x 7
L0_DIM e(permdist)   = 120 x 1
L0_DIM e(perm_jointW)= 120 x 1
L0_DIM e(perm_advW)  = 120 x 1
L0_PERMN min=97 max=120 completed=120       <-- 23 replications lost >=1 element
L0_SAVE_ROWS=121  (rep 0 + 120)  ONSETCOLS=30  COEFCOLS=7
L0_FILE_FAILEDROWS=0
L0_FILE_PARTIALROWS=23
printed: "120 of 120 replications completed"; joint p = 0.8367; adversarial p = 0.9694
NO caution about the lost elements (e(perm_incomplete)/e(perm_addone) are gated on
permexact, ado:16130-16131, so a SAMPLING run records the loss only in e(perm_n))
```

**This is a shape neither earlier repro produced, and it inverts the finding's premise.**
Because no replication failed OUTRIGHT, `wooldid_bootdist`'s `rowsum(A :< .) :> 0` mask
drops nothing: `e(permdist_b)` is 120 x 7, the SAME row count as `e(perm_jointW)`. So on
this panel the two objects share a row basis and the "different lengths" tell is ABSENT —
yet `e(perm_jointW)` still carries missing rows (the replications that lost a tested pre
horizon), and `e(perm_failed)` is 0, which tells the user there is nothing to drop.
Run B measures how many and what the naive reconstruction costs.

## Results — run B (`L0b_permrows.do` / `.log`), arithmetic complete
Same data, same seeds, reproduces run A.
```
LB_reps=120 completed=120 failed=0        e(perm_incomplete) and e(perm_addone) EMPTY
LB_DIMS permdist_b=120x7 permdist=120x1 jointW=120x1 advW=120x1   <-- all the SAME length
LB_FILE_ALLMISS=0   LB_FILE_PARTIAL=23
LB_COLMISS ES_att_pre4 = 23 ; ES_att_pre3 = 23 ; every other column 0
LB_C1 maxabsdiff(permdist_b, permsave coefficient block) = 0 over 794 finite cells
LB_C1 missing cells: file=46  permdist_b=46      (23 x 2 lost horizons)
LB_C2 jointW joint_horiz: rows=120 MISSING=23, while e(perm_failed)=0
LB_C3 rows missing in W but KEPT in the file = 23 ; present in W but all-missing in file = 0
LB_C4 correct  n=97  nex=81  p=0.83673469 (+1.ac687d6343eb2X-001)
LB_C4 naive    N=120 nex=104 p=0.86776860 (+1.bc4c2a50658dcX-001)  overcount=23
      e(perm_failed)=0   reps - n_joint = 23
LB_C5 POSTED joint_horiz = 0.83673469 (+1.ac687d6343eb2X-001)  matches_correct=1 matches_naive=0
LB_ADV adversarial_pval_att: 23 missing rows; n=97; correct 0.96938776 = POSTED;
      naive 0.97520661; overcount 23
```

## Verdict — PARTIAL (severity low, realistic true)

| claim | my measurement |
|---|---|
| permdist_b/permdist compacted, jointW/advW full-length | TRUE in source (`bootdist` ado:14644 `select(A, ok)` vs `Wr = J(rows(A),1,.)` ado:16208/16402) and in the two earlier runs. **But not a reliable tell:** my run has zero outright failures, so compaction drops nothing and all four matrices are 120 rows, while the hazard is still live and larger. |
| "row r of one is a different replication from row r of the other" | only when a replication fails OUTRIGHT. Not so here. |
| "neither is indexed by permsave's _rep" | **FALSE for perm_jointW/perm_advW** — measured: their 23 missing rows are exactly the 23 file `_rep`s that lost the tested horizons, 0 violations. Also false for `permdist_b` here (maxabsdiff 0 against the whole file block). It is true of `permdist_b` only in the outright-failure regime. |
| naive reconstruction over-counts | **TRUE, and it is the real content of the finding.** Mata `:>=` ranks missing above every finite value: 104 vs 81 exceedances, p 0.8678 vs 0.8367. |
| "over-counts by exactly e(perm_failed)" | **FALSE.** Measured over-count 23 with `e(perm_failed) = 0`. The quantity is `reps - n_test`, and `n_test` (permqf/permadv's `outn`) is **never ereturned** — `e(perm_n)` is per COEFFICIENT, not per test. |
| "no printed number is wrong" | **CONFIRMED.** Posted `joint_horiz` and `adversarial_pval_att` equal the correct hand computation bit-for-bit (hex match) and differ from the naive one. |
| "the help states this for the bootstrap side and says nothing for the permutation side" | **Overstated.** sthlp:493 does say "a failed replication keeps its row, all missing" for permsave. What is genuinely missing: (a) the e()-side bridge sentence bootsave's entry carries (sthlp:557 "so _rep indexes the draws exactly as e(bootdist) does after its missing rows are dropped"); (b) any statement that the W matrices carry missing rows at all; (c) permsave writes NO `jointW_*`/`advW` columns (5-arg `bootsavemk` at ado:8859) although bootsave does under `boottstat` (11-arg at ado:6546), while the permsave entry closes "Every reported permutation quantity is otherwise reconstructable from the file." |

### Mine, not seeded
1. The **partial-replication regime**: `e(perm_failed) = 0`, all four replication matrices the same length, the run prints "120 of 120 replications completed" and **no caution at all** — yet 23 W rows are missing. Every tell the finding names is absent there, and the error is the same size.
2. `e(perm_incomplete)` / `e(perm_addone)` are gated on `M < .` (`permexact`) at ado:16130-16131, so a SAMPLING permute() run has no per-element incompleteness record beyond `e(perm_n)`, and none at all for a test.
3. `e(perm_sgW)` (built by the same `permqf` at ado:8535) shares the property; the claim names only jointW/advW.
4. The permsave/bootsave W-column asymmetry above.

### Caveats on my own evidence
- In run B every file row was kept (`allmiss ≡ 0`), so the "0 violations" leg of LB_C3 is
  arithmetically vacuous. What I actually have is a COUNT agreement (23 W-missing rows =
  23 rows missing ES_att_pre4/pre3) plus `permqf`'s `ok = rowmissing(A[., idx]) :== 0`,
  which forces the identity. The row-by-row identity on OUTRIGHT failures is the skeptic's
  measurement, not mine.
- `e(perm_sgW)` sharing the property, and the bootstrap-side `e(boot_jointW)` twin, are
  SOURCE-READ ONLY here (ado:8535 and ado:16351); not measured.
- The 23/120 size of the over-count is design-dependent (3 onsets, 18 of 30 districts
  entering the frame late). The point is that it is PRESENT with `e(perm_failed)=0`, not
  that it is large.
- Adjacent, uncertain, NOT folded into this finding: a sampling run that loses elements on
  COMPLETED replications prints "120 of 120 replications completed", no caution, and no
  `e(perm_incomplete)`; `e(perm_n)` is the only record. That may well be by design, since
  the help documents `e(perm_n)` as "replications carrying the element".

### The sentence that survives
`e(perm_jointW)` / `e(perm_sgW)` / `e(perm_advW)` carry a missing row for any replication
INCOMPLETE ON THE TESTED COLUMNS, not only for a failed one. The help's advertised
reconstruction route ("the p is reconstructable from the pair", sthlp:1349) neither says to
drop those rows nor posts the per-test n, and `e(perm_failed) = 0` with matching row counts
against `e(permdist_b)` does not mean there is nothing to drop.
