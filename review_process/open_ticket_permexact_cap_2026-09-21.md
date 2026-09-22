# TICKET — let safety(off) lift the permexact enumeration cap

**Opened 2026-09-21 on the author's instruction,** left as a reviewer note on card #240 of the
message review page:

> *"Please update this so that it exists as an error by default, but that when safety(off) is set,
> it appears merely as a caution. Please also include a caution that this means the permutation
> matrix won't post to e() - converting that other guard to a caution also. Please include a 'Set
> safety(off) to disable this rep cap.' statement at the end of the error text here."*

**Status: OPEN. Not built, and it should not be built as written.** One measured fact changes the
design, and four user-facing strings need the author's wording before anything ships, under his own
rule of the same day.

## The fact that changes the instruction

**`safety(off)` cannot disable the cap, because a harder ceiling sits underneath it.** Measured, not
inferred: Stata matrices stop at **65,534 rows**, and a request for one row more fails r(915). The
replication accumulators are Stata matrices sized by the replication count, so a run with more than
65,534 assignments dies with an opaque r(915) **after** the full enumeration has already been built
and paid for.

There is a second unliftable case. When the assignment count overflows a double, the count becomes
missing by design, and then neither the assignment-matrix allocation nor the replication loop can
run at all.

So the current single test fuses three different situations, and only one of them is about the
user's patience:

| Case | Today | Should be |
|---|---|---|
| Count overflows to missing | Refused | **Refused, unliftable** |
| Count above 65,534 | Refused | **Refused, unliftable** |
| Count between 10,001 and 65,534 | Refused | Refused by default, **downgradable by `safety(off)`** |

The sentence the author dictated, "Set safety(off) to disable this rep cap", is therefore accurate
only on the third row. On the first two it would promise an escape that does not exist and send the
user to an r(915) they cannot diagnose.

## What the work is

**Thread the safety flag.** The routine that holds the cap is Mata and cannot see the safety
setting. Four small edits carry it: add the option to the replication runner's call in the shell,
add it to that program's `syntax`, pass it as a new argument on the Mata call, and add the parameter
to the Mata signature. Use the user's own typed value, matching the convention the other gated
screens use. The validation program needs **no** change: it never reaches this routine.

**Split the fused test into three branches,** per the table above, with a second constant beside the
existing one for the Stata ceiling.

**The other guard is two places, not one.** The draw matrices are dropped above ten thousand
completed replications by two separate hard-coded literals, one in the shell for the main matrix and
one in the rewiring program for the two companions. Neither shares storage with the enumeration cap;
they coincide only by hand maintenance. Both are **silent today**. The author asked for a caution,
so both get one, and the two literals should be pulled into one named constant in the same edit so
they cannot drift.

Worth knowing: that posting cap is **already reachable today** through a large `permreps()`, which
already drops the matrices without a word. So this caution fixes an existing silent loss, not only
one created by lifting the enumeration cap.

## The four strings that need the author's wording

Under the rule he set the same day, these are his to write, not mine. Drafts are offered only as a
starting point.

1. **The default refusal, with his dictated sentence appended.** His words: *"Set safety(off) to
   disable this rep cap."* The rest of the sentence already exists.
2. **The caution when `safety(off)` lifts the cap.** Must say the run is now enumerating a large
   number of assignments and will be slow.
3. **The caution that the draw matrices will not post.** Must name what is lost and that
   `permsave()` still writes them, since it reads from the accumulators and is not subject to the cap.
4. **The two unliftable refusals.** These are new cases with no existing text, and they must not
   repeat the "set safety(off)" offer.

## Test position

**No test anywhere asserts the current cap text**, so nothing moves and nothing needs re-blessing.
That also means the new behaviour would ship untested, so the same commit should add: a run just
above the cap refused by default and admitted under `safety(off)`; a run above the Stata ceiling
refused at both settings; and an assertion that the draw-matrix caution fires when the matrices are
dropped.

## Recommendation

Build the three-branch split and the threading, which are unambiguous. Put the four strings to the
author first, and put the 65,534 ceiling to him with them, because it narrows what his dictated
sentence can truthfully promise.
