# TICKET — raise the permexact enumeration cap to the Stata matrix ceiling

**Opened 2026-09-21 on the author's instruction,** left as a reviewer note on card #240 of the
message review page:

> *"Please update this so that it exists as an error by default, but that when safety(off) is set,
> it appears merely as a caution. Please also include a caution that this means the permutation
> matrix won't post to e() - converting that other guard to a caution also. Please include a 'Set
> safety(off) to disable this rep cap.' statement at the end of the error text here."*

**Status: BUILT 2026-09-22, in a narrower form than this ticket proposed. One item remains for
the author (s3).**

He ruled: *"For the rep cap, let's just use the matrix cap that you listed out."* That simplifies
the design away. Instead of threading the safety flag and splitting the test three ways, the cap is
simply moved to where the machinery actually stops:

    capexact = 65534      (was 10000)

That single constant change does everything the three-branch design would have. The `M >= .`
overflow arm of the existing test is untouched and still catches the other unrunnable case. No
safety flag is threaded, no program signature changes, and the two refusals that would have been
"unliftable" are now just the one refusal, which was never liftable.

**The ceiling is measured, not chosen.** `J(65535, 3)` fails r(915) and `J(65534, 3)` does not, and
the replication accumulators are Stata matrices with one row per replication. One assignment more
and the run dies r(915) after paying for the whole enumeration.

**The draw-matrix caution asked for in his earlier note IS built**, because raising the cap is what
makes it matter: an exact run between 10,001 and 65,534 assignments is now legal and silently loses
`e(permdist)`, `e(permdist_b)` and `e(permdist_se)` to the separate 10,000 posting cap. That cap now
appears three times and they share no storage, so `check_validate_sync.py` gained a static check
that all three agree.

**What is verified.** `test_permute.do` P12a calls the enumeration routine directly on a synthetic
8-cluster frame giving 40,320 assignments, above the old cap and below the new one, and asserts it
enumerates rather than refuses. P12b asserts a 12-cluster design is refused with the NEW ceiling
named. P12c refuses end to end through the command. The direct-call approach is deliberate: testing
the opened band through the command would mean actually running 40,320 replications.

## Still for the author

**The wording of the new caution is mine, not his,** and is flagged on the review page. Draft:

> Caution: N completed replications exceeds the 10,000 cap on posting draw matrices, so
> e(permdist), e(permdist_b) and e(permdist_se) are not posted. Use permsave() to write the draws.

His dictated sentence, *"Set safety(off) to disable this rep cap."*, is **not** used, because under
this ruling `safety(off)` has nothing to do: the cap is now a machine limit rather than a policy,
and nothing can lift it. Confirm that is what he intended.

One observation, not built: `e(bootdist)` has the identical 10,000 posting cap on the bootstrap
side and is still silent. Nothing about this ruling forces a change there, but the asymmetry is now
deliberate rather than accidental.

## The original design, kept for the record

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
