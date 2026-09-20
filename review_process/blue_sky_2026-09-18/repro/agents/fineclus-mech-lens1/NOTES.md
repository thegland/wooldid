# FCM-2 — third lens (DOCUMENTED-BEHAVIOUR). Working notes.

## Question I am answering
Is the regime-predicate-on-`smp` behaviour a DEFECT, or is it documented / ruled /
correct-but-surprising design?

## Record check (done, read-only, no suite files used)
- `feature_tickets_2026-09-15.md`: T1/T2/T3 "COMPLETE, merged and blessed" 2026-09-16;
  **"The full assessment suite is deferred to 2026-09-17 by the author; it is the only
  outstanding item."** So this path is one day old and its assessment has NOT run.
  known_bydesign §7 lists it as a hunting ground, not as ruled.
- `git show 456f636` (T2/T3 machinery, 2026-09-16 09:54): states the design choice
  explicitly — "wooldid_bootrun now decides the regime on the same `smp' the resample
  runs on". Nowhere does it consider smp vs e(sample) / markout. Same message states the
  invariant that is violated: "keeping i under nesting would make a twice-drawn cluster's
  cells identical duplicates that the equal-cell weighting counts once, collapsing the
  composition variance of att_i toward zero."
- `git show 60e31c0` (gate retirement): confirms the A10 caution is now conditional on
  the nesting regime; nothing screened in place of the retired refusals.
- A7 (`git show 1ad2f65`, 2026-08-30) is the ruling for the RESAMPLE FRAME's sample
  concept ("user-if sample + nonmissing cluster"), written when no predicate existed.
  It is a ruling about a different object.
- NOT read: `test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md` (off-limits per brief;
  only test/suite/README.md is whitelisted). Noted as unread.

## Help-file check
- sthlp:532 states the rule with NO sample qualification: "what becomes of the unit
  identities on that frame depends on where cluster() sits relative to i. Where the
  cluster nests i, each drawn copy of a cluster becomes a fresh set of units".
- The help uses "the estimation sample" for this engine's own sample concept elsewhere
  ("bootstrap() requires at least two clusters in the estimation sample"; bootsave's
  "_rep = 0 row holds 1/N, N the number of observations in the estimation sample").
  Nowhere does the help say the bootstrap's frame is the if-sample before markout.
- No `{marker}` section documents a bootstrap-specific sample concept; enginediff says
  nothing about which rows.

## Observability check (grep, to be confirmed in Stata)
- `wdb_nest` is `c_local`-returned and used ONLY for the A10 caution, which additionally
  requires an i-keyed covariate in fe()/controls()/... — so with no such covariate the
  regime switch prints NOTHING.
- No e() scalar/macro records the regime. Contrast `e(perm_irebuilt)`, which the author's
  2026-09-16 ruling on the permute display line cited as "the record" for which path a run
  took. The bootstrap leg has no analogue.

## Prober's one overread (skeptic is right)
ado:8975-8985 defines "the sample" as `smp` in its own preceding sentence, so it is
internally correct about raw-frame-vs-smp and never claims smp is the estimation sample.
The correct statement is narrower and still damning: the comment ARTICULATES the principle
("the predicate is computed on the rows bsample actually sees ... a cluster that nests i
within the sample but straddles two i's only through rows the sample excludes would be
misclassified") and the implementation applies it one level too shallow.

## My own probe (p_lens1.do): three arms
A base (nest) / B phantom cluster change on the missing-y rows (keep-i) / C = B with those
rows PHYSICALLY DELETED (nest). A, B, C have identical e(sample), N, b, analytic SE.
If se_boot(A) == se_boot(C) != se_boot(B), the bootstrap SE depends on rows the estimator
never uses, and the help's "where cluster() sits relative to i" is identical in B and C on
every estimated row.

## Status
Record + help check DONE. Stata probe: WRITTEN, RUNNING.
Tentative verdict: CONFIRMED, high, realistic-but-narrow.

## RESULT (p_lens1.log, complete, rc 0 throughout)

Three arms, one panel (96 firms in 24 industries, 8 years):

| arm | rows | e(sample) | NEST on smp | NEST on e(sample) | att_i (hex) | analytic se (hex) | boot se(att_i) |
|---|---|---|---|---|---|---|---|
| A base            | 768 | 764 | 1 | 1 | +1.02b32930bf451X+000 | +1.4a2b00da9e50eX-005 | 0.14073971 |
| B phantom         | 768 | 764 | **0** | 1 | +1.02b32930bf451X+000 | +1.4a2b00da9e50eX-005 | **0.09955304** |
| C phantom deleted | 764 | 764 | 1 | 1 | +1.02b32930bf451X+000 | +1.4a2b00da9e50eX-005 | 0.14073971 |

- **A and C are BIT-IDENTICAL on the bootstrap SE** (+1.203c2471fb1aeX-003 both). The four
  phantom rows are provably inert to the estimator AND to the draws: deleting them changes
  nothing. The ONLY thing that moves the number is the cluster VALUE written on rows the
  regression never sees. B is 0.7074 of A/C — a 29.3% understatement, anti-conservative.
- **Completely silent.** `diff out_A.txt out_B.txt` over the whole console block returns only
  the numbers themselves (se, z, CI, percentile row). No Caution anywhere in the run.
- **No e() record.** `eREGIMEMARKER=0` in all three arms (no scalar or macro whose name
  contains nest/regime/reid/rebuil); `e(scalars)`/`e(macros)` name lists are IDENTICAL between
  A and B (enames_A.txt vs enames_B.txt diff clean). The A10 caution is additionally gated on an
  i-keyed covariate, so even that channel is shut on an ordinary specification.

## One further help sentence, found after the run (the `if` entry, sthlp:61)
"The rule is deliberate ... and it applies identically to the estimation sample, to the
randomization design under permute() and to the bootstrap frame, **so the three always agree**."
Scoped to the `if` rule, so not a blank promise that the three row sets are equal — but it is the
ONLY place the help relates the bootstrap frame to the estimation sample, and it points the wrong
way for FCM-2.

## VERDICT: CONFIRMED. Not documented, not ruled, not correct-but-surprising.
