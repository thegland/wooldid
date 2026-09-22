# TICKET — make the subgroup containment guards unliftable

**Opened 2026-09-21 on the author's instruction.** First stated as *"I'm seeing some references
to a share weighting feature used when a given level of i crosscuts multiple subgroups. Please
open a ticket to delete this feature."* Rescoped the same day once the survey showed there is no
feature to delete: *"don't worry about removing features, let's just update the ticket to harden
our guards for this behavior to always fire."*

**Status: PARTIALLY BUILT 2026-09-21, and the remaining scope is now in doubt.**

The author ruled on the message-review page rather than on this file, and the ruling is narrower
than this ticket assumed. He unsilenced exactly ONE of the three containment guards, S3, the
permute() split-unit screen: card #180, `silenced_by_safety` set false. **S3 is built** and pinned
by a new G5 leg in `test/test_sgsplit.do`, which asserts r(459) at both safety settings.

He left S1 and S2 alone, and on the same pass he REWROTE both of their messages while deliberately
keeping the escape clause in the text: the worker's split-unit and within-cell errors still end "To
ignore this error, specify safety(off)." and the partial-labelling error still ends "Specify
safety(off) to proceed anyway." So his most recent and most specific signal is that the worker's
screens keep their escape, which contradicts s1 to s3 of this ticket as originally written.

**S1 and S2 are therefore NOT built, pending his word.** Building them would have made those three
retained sentences false, which is the defect class the message review exists to remove. If he does
want them hardened, the escape clauses in all three messages have to go in the same edit.

Line numbers are as of `1323749` on `bugfixes`.

## The rule, and what the share-weighting references actually are

The rule is containment: **each level of `i` must lie wholly inside one subgroup**, so subgroups
are unions of whole cohorts. wooldid enforces it by default and refuses the run when it fails.

Every reference to share weighting is text hung off that refusal. The eight mentions in the ado
are three error messages, three cautions and two code comments. **None is executable.** The
estimand weights never branch on the subgroup variable, and the two temporary variables the
screens build are asserted on and then dropped.

Share weighting is just what the arithmetic does once a split unit is let through. The `att_i`,
`att_it` and `att_itime` denominators are built over the whole estimation sample while
`subgroups()` only selects the `margins` subpop, and `margins` renormalises inside it. A unit whose
labelled rows carry share *s* of its weight lands in its subgroup with total weight *s* rather than
1. Nobody wrote that behaviour and nothing selects it.

So this ticket removes no logic. It removes the `safety(off)` bypass from three guards, so the
containment rule holds unconditionally.

## Site inventory

| # | Program | Anchor | Condition caught | Today, default safety | Today, `safety(off)` | After |
|---|---|---|---|---|---|---|
| S1 | `wooldid_sgdom` | 12759-12771 | Unit's subgroup changes over time | Error r(459) | Caution, run proceeds | Error r(459) always |
| S2 | `wooldid_sgdom` | 12803-12817 | Unit labelled on some rows, missing on others | Error r(459) | Caution, run proceeds | Error r(459) always |
| S3 | `wooldid_permctl` | 6392-6402 | S1's condition, over all units | Error r(459) | Caution, run proceeds | Error r(459) always |

Each site has the same shape: an `if "`safetyinput'" != "off"` arm that errors and an `else` arm
that cautions. The edit drops the condition and deletes the `else` arm.

Two things to preserve:

- **S2 keeps its `"`rm'" != "1"` replication gate.** Under `permonsetcohorts` a replication that
  rebuilds `i` can pool a labelled unit with a wholly unlabelled one, which is a resampling
  artefact rather than something to report. The observed run's screen is the one that matters.
- **Wholly unlabelled units stay legal.** S2 fires only when `f13sgimin < .`, so a unit with no
  label on any row still belongs to no subgroup, exactly as the help promises.

## Message rewording is part of the work, not a follow-up

All three errors currently end by offering the escape that this ticket removes:

- S1: *"To proceed anyway, specify safety(off): by-subgroup att and att_it then weight each
  observation to its own subgroup as usual, while att_i and att_itime enter a split unit into each
  subgroup weighted by its share of sample weight there..."*
- S2: *"To proceed anyway, specify safety(off): the affected rows are then share-weighted, exactly
  as a split unit is."*
- S3: *"To proceed anyway, specify safety(off); share weighting then applies inside the
  replications as it does in the observed run."*

Those clauses become false on the day the guard hardens, so they go with it. A message that names
an option which no longer changes anything is the exact defect class the message review is
cataloguing. What should survive is the **why**: each error should still say that the unit would
otherwise enter its subgroup by weight share rather than with the equal-unit weighting `att_i` and
`att_itime` name. That is the justification for refusing, and it is worth keeping even though the
user can no longer choose it.

The remedies the errors already list are unaffected and stay: define `i` at the finer level with
`egen newi = group(i sg), missing`, blank the label on all of a unit's rows, fill it in, or
restrict with `if()`.

## Precedent for the tier

This is not a new kind of check. The `permonsetcohorts` bijection verification already sits at an
unliftable tier, with the comment at 6281 stating it is *"a hard error at the same tier as the
onset-constant-within-cluster check above and is not reachable with safety(off)."* Its rationale
transfers: without the precondition the downstream object is not well defined, so there is nothing
coherent for `safety(off)` to let the user proceed into. Containment is the same shape. Without it
the estimand `att_i` names is not the estimand the run computes.

## Impact, measured

**No shipped suite config exercises the path.** Four configs pair `safety(off)` with `subgroups()`
(`sg_attit`, `sg_attitime`, `cont_sg`, `sg_missif`), but all carry `safety(off)` for the unrelated
SI-A1 cluster-at-`i` joint-test refusal ruled 2026-08-25, and every one of their subgroup variables
is unit-constant:

| Variable | Definition | Verdict |
|---|---|---|
| `sg` (dta_binary and friends) | `1 + mod(i,2)` family | Constant within `i` |
| `gen_sg2` | `mod(i,2)+1` | Constant within `i` |
| `sgm_missif` | `mod(i,3)`, blanked where `mod(i,11)==0` | Blanks whole units, so legal |

The `sg = 1 + mod(_n,2)` at `mkdata_suite.do:151` is generated **before** `expand 12`, when `_n`
indexes units, so it too is unit-constant.

Consequences:

- **Battery must come back IDENTICAL.** No config changes behaviour, so any movement is a defect in
  the edit.
- **`test/test_sgsplit.do` G2 inverts.** It currently asserts that a split unit under `safety(off)`
  runs and prints the share-weighting caution. It must now assert the run fails. G1, G3 and G4 are
  untouched. The file header states the current contract in its first three lines and needs
  rewriting.
- **`test/suite/out/degeneracy/cautions.csv`** carries three share-weighting rows. Regenerate, do
  not hand-edit.
- **Users relying on the escape lose it.** This is the intended breaking change. Worth a line in
  whatever release note the author keeps.

## Documentation

**Help file**, three passages, each of which currently says the run proceeds under `safety(off)`:
the safety discussion at 242, the `subgroups()` entry at 622, and the subgroup discussion at
645-649. Keep the edits minimal and in voice: new sthlp prose is measured against `2a2da17` on
countable markers per the documentation review's register trap.

**Agent guide**, prose at 473, 486 and 494, plus the generated message rows at 980, 1077, 1078,
1223 and 1317. Regenerate chapter 14 from the edited ado with `scratchpad/gen_messages_chapter.py`
rather than editing tables by hand.

**Message review page**, four cards correspond to messages that change or disappear. They reseed
from the ado via `msgreview/consolidate.py` on the next rebuild.

## Open question for the author

**Q1. Do the two within-cell guards harden too?**

Two further screens catch a subgroup that varies **within** a single `i`-by-`t` cell:
`wooldid_sgdom` at 12741-12754 and its `permute()` mirror `wooldid_permctl` at 6383-6391.

They differ from S1 to S3 in two ways. They sit wholly inside `if "`safetyinput'" != "off"`, so
`safety(off)` **skips the check** rather than downgrading it to a caution. And their rationale is
not share weighting at all: subgroup contrasts built from mixed cells cannot reflect within-cell
heterogeneity, and can come back identical even under strong true heterogeneity. That is an
uninformative-answer problem, not a wrong-weights problem, which is why no share-weighting language
attaches to them.

- **Option A, the narrow reading.** Harden S1 to S3 only. The within-cell guards keep their bypass.
- **Option B.** Harden all five, giving one rule with no exceptions: a `subgroups()` variable must
  partition units cleanly, and `safety(off)` never lifts it.

The instruction says "this behavior", which points at the containment guards just discussed, so
**this ticket is scoped to Option A unless the author says otherwise.** Flagging it because the
help file's subgroup discussion at 645-646 bundles all five into one sentence, so the prose edit
has to take a position either way, and because a half-hardened rule is harder to state than either
end.

## Verification plan

1. Edit the three sites and their message text. Rebuild, run `test/check_validate_sync.py`.
2. Invert `test_sgsplit.do` G2 and rewrite its header. Run `test_sgsplit`, `test_sgdomain` and
   `test_fastcollapse`.
3. Pin the three retired bypasses the way the finer-than-`i` tickets did: each config measured at
   rc 0 before the change and r(459) after.
4. Full battery, expecting **IDENTICAL**.
5. Regenerate the degeneracy CSV and guide chapter 14, then rerun the coverage audit.
