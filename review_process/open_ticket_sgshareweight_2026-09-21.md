# TICKET — remove the subgroup share-weighting fallback

**Opened 2026-09-21 on the author's instruction:** *"I'm seeing some references to a share
weighting feature used when a given level of i crosscuts multiple subgroups. Please open a
ticket to delete this feature."*

**Status: OPEN. Not started. One scope question for the author before any edit (s5).**

Line numbers below are as of `1323749` on `bugfixes`.

## What the feature is

When a `subgroups()` variable does not partition units cleanly, wooldid refuses the run under
default safety and offers `safety(off)` as an escape. Under that escape the run proceeds with
**share weighting**: a unit that touches more than one subgroup enters each one weighted by the
share of its sample weight sitting there, instead of with the equal-unit weighting that `att_i`
and `att_itime` are defined to carry. By-subgroup `att` and `att_it` are unaffected.

Two distinct data conditions reach it:

- **Split unit.** The subgroup value changes over time within a level of `i`.
- **Partial labelling.** The subgroup value is present on some of a unit's rows and missing on
  others. The split-unit screen folds missing through `cond(sg >= ., ., sg)`, so this passes that
  screen and needed its own.

## The finding that decides what "delete" means

**There is no share-weighting code to delete.** Share weighting is not an algorithm wooldid
selects. It is what the ordinary weight machinery already does when a unit is split, and the
screens are purely diagnostic.

The estimand weights `w_att_i`, `w_att_it` and `w_att_itime` carry denominators (`w_per_i`,
`w_per_it`, `t_per_i`) built at the weights block over the whole estimation sample, while
`subgroups()` only selects the `margins` subpop and `margins` renormalises within it. A unit whose
labelled rows carry share *s* of its weight therefore enters its own subgroup with total weight
*s* rather than 1, automatically.

Confirming this in the code: the screens create `f13sgimin` and `f13sgimax`, run `assert` against
them, and `drop` them at 12818. They are never read by anything that builds a weight.

So deleting the feature means **hardening three refusals so `safety(off)` can no longer lift
them**, and deleting the three caution branches that currently describe the fallback. No weighting
logic changes, and no reported number moves on any run that is legal today.

## Site inventory

| # | Program | Anchor | Condition | Today under default safety | Today under `safety(off)` |
|---|---|---|---|---|---|
| S1 | `wooldid_sgdom` | 12759-12771 | Split unit | Error r(459) | Caution, share weighting |
| S2 | `wooldid_sgdom` | 12803-12817 | Partial labelling | Error r(459) | Caution, share weighting |
| S3 | `wooldid_permctl` | 6392-6402 | Split unit, over all units | Error r(459) | Caution, share weighting |

S2 is additionally gated on `"`rm'" != "1"`, so it does not run inside replications. That gate is
deliberate and should survive the change: under `permonsetcohorts` a replication that rebuilds `i`
can pool a labelled unit with a wholly unlabelled one, which is a resampling artefact rather than
something to tell the user about.

S3 is the `permute()` mirror of S1 and runs over every in-sample row rather than treated units
only, because any cluster can be assigned treatment in a replication.

## What moves with the deletion

**Tests.** `test/test_sgsplit.do` is the dedicated suite, opened on the author's directive of
2026-07-28. Its G2 leg asserts that a split unit under `safety(off)` **runs** and prints the
share-weighting caution. That leg inverts: the run must now fail. G1, G3 and G4 are unaffected. The
file's header comment states the current contract in its first three lines and must be rewritten.

`test/test_fastcollapse.do` also references the condition and needs a read before the edit.

**Suite artifacts.** `test/suite/out/degeneracy/cautions.csv` carries three share-weighting rows.
Regenerate rather than hand-edit.

**Help file.** Three passages in `wooldidinstall/wooldid.sthlp` describe the fallback: the safety
discussion at 242, the `subgroups()` entry at 622, and the subgroup discussion at 645-649. Each
currently says the run proceeds under `safety(off)`; each becomes an unconditional refusal. Note
the register constraint from the documentation review: new sthlp prose is measured against
`2a2da17` on countable markers, so keep the edits minimal and in voice.

**Agent guide.** `wooldidinstall/wooldid_agent_guide.md` carries the prose at 473, 486 and 494, and
the generated messages chapter carries the four message rows at 980, 1077, 1078, 1223 and 1317.
Regenerate chapter 14 from the edited ado with `scratchpad/gen_messages_chapter.py` rather than
editing those tables by hand.

**Message review page.** Four cards on the review artifact correspond to the messages that change
or disappear. They are reseeded from the ado by `msgreview/consolidate.py`, so this lands there
automatically on the next rebuild.

## Open question for the author

**Q1. Do the two within-cell screens go with it?**

There is a fourth and fifth site that are adjacent but behave differently, and the instruction does
not obviously reach them:

- `wooldid_sgdom` at 12741-12754: the subgroup varies **within** an `i`-by-`t` cell.
- `wooldid_permctl` at 6383-6391: the `permute()` mirror of the same.

These sit entirely inside `if "`safetyinput'" != "off" {`, so `safety(off)` **skips the check**
rather than downgrading it to a share-weighting caution. Their stated rationale is also different:
subgroup contrasts built from mixed cells are structurally unable to reflect within-cell
heterogeneity and can come back identical under strong true heterogeneity. That is a
"the answer is uninformative" problem, not a "the weights are not what the estimand says" problem.

They are worth deciding together because the help file's subgroup discussion at 645-646 bundles all
of them into one sentence, so the prose edit has to take a position either way.

- **Option A.** Delete only the share-weighting fallback, S1 to S3. The within-cell screens keep
  their `safety(off)` bypass.
- **Option B.** Harden all five, making every unclean subgroup partition an unliftable refusal.

Option B gives the simpler user-facing rule, which is that a `subgroups()` variable must partition
units cleanly with no exceptions. Option A is the narrower reading of the instruction as written.

## Verification plan once ruled

1. Edit the sites, rebuild, `check_validate_sync.py`.
2. Invert `test_sgsplit.do` G2, rerun `test_sgsplit`, `test_sgdomain` and `test_fastcollapse`.
3. Confirm the three retired `safety(off)` configs measure rc 0 before and r(459) after, the shape
   the finer-than-`i` tickets used.
4. Battery must come back **IDENTICAL**: no legal run changes, so any movement is a defect in the
   edit.
5. Regenerate the degeneracy CSV and guide chapter 14; rerun the coverage audit.
