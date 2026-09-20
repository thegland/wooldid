# FINDINGS -- perm-semantics
## (semantic conformance: every permute() p-value and interval derived from the saved null distribution)

METHOD.  For every reported permutation quantity I re-derived the number myself from the
saved null distribution (`permsave()` file, `e(permdist_b)`, `e(permdist_se)`,
`e(perm_comparator_b/_se)`) using ONLY the help file's definitions, and compared it
bit-for-bit with what wooldid posted (`e(perm_kernel)` rows 1-18 and the stripe surfaces
`e(perm_p)`, `e(perm_ci)`, `e(perm_p_lo/hi)`, `e(perm_p_se)`, `e(perm_p_ci)`,
`e(perm_sd)`, `e(perm_nullband)`, `e(perm_n)`, `e(perm_se_shift)`, `e(perm_tlb_*/tub_*)`,
`e(perm_stepdown*)`, `e(perm_jointW/sgW/advW)` + their `*obs` twins,
`e(perm_pwsgs_kernel)`), plus `e(p)/e(ci_l)/e(ci_u)` and `e(wdidmainresults)`.
p-values were required to match EXACTLY (they are rationals of integer counts);
quantile-derived numbers to rtol 1e-12.

**The engine is conformant.** Across 16 probes, five p-value conventions, both sampling
modes, permtstat, permanalyticse, permtci, permstrata, permonsetcohorts, subgroups,
element-level holes and whole-replication failures, NOT ONE reported p-value, interval,
null band, SD-of-draws, Monte-Carlo standard error, Clopper-Pearson interval, step-down
adjusted p or quadratic-form test p differed from the documented definition. The two
findings below are both about how the REPLICATION-LEVEL objects are stored, not about a
reported number, and both are rated low for that reason.

---

## F1 (low) -- `e(permdist_b)` and `e(perm_jointW)`/`e(perm_advW)` are not row-aligned when replications fail

**Claim.** In a single `permute()` run with failed replications, `e(permdist_b)` and
`e(permdist)` are compacted to the completed replications while `e(perm_jointW)` and
`e(perm_advW)` keep the full requested length with missing rows. Row *r* of one is
therefore a different replication from row *r* of the other, and neither is indexed by
`_rep` in the `permsave()` file. Nothing in `e()` or in the help says so.

**Repro.** `/private/tmp/.../probes/perm-semantics/repro_permdist_rowalign.do`
(self-contained; 12 clusters, `permute(probabilities)` where a draw can leave zero
treated clusters).

**Observed.** requested=120, completed=111, failed=9:
```
e(permdist_b)   : 111 rows x 6 cols
e(permdist)     : 111 rows x 1 cols
e(perm_jointW)  : 120 rows x 1 cols   (111 non-missing)
e(perm_advW)    : 120 rows x 1 cols
permsave file   : 121 rows  (120 replications + _rep = 0)
```

**Expected.** Either the same row basis for every replication-indexed object of one run,
or a statement of the difference. The bootstrap side of the same help file carries
exactly that statement -- "a replication that failed keeps its row with every value
missing, so `_rep` indexes the draws exactly as `e(bootdist)` does after its missing rows
are dropped" (`bootsave()`) and "`e(bootdist_se)` ... on the same row mask as
`e(bootdist)`". The permutation paragraph says only "the draw matrices `e(permdist_b)`
(and `e(permdist_se)` under `permtstat`)" and, for the test statistics, "the replication
statistics behind the joint, equality and adversarial p-values, one named column per
test". The help also advertises this surface as the reconstruction route ("Every reported
permutation quantity is otherwise reconstructable from the file", "the p is
reconstructable from the pair"), so pairing a coefficient draw with its W is exactly the
use the documentation invites.

**Consequence.** No printed number is wrong; this is a stored-results hazard. `select()`
on the mismatched pair errors loudly, so the dangerous shapes are the quiet ones: a
positional loop (`for r=1..rows(permdist_b) { ... jointW[r] ... }`), and merging
`e(permdist_b)` onto the `permsave` file by row number -- both silently read the wrong
replication from the first failure onward. A sharper case: the help's own reconstruction
formula for the joint p, `(1 + #{W_r >= W_obs})/(n+1)`, applied to the full-length
`e(perm_jointW)` without first dropping its missing rows, over-counts by exactly
`e(perm_failed)`, because Mata's `:>=` treats missing as +infinity -- so on this run a
reader gets `(1+9+nex)/121` instead of `(1+nex)/112`.

**Realistic.** Yes. Failures are not pathological: any `permute(probabilities)` design in
which an onset can go undrawn, and any shuffle design where a reassigned onset makes a
draw non-estimable, produces them; the help gives failures their own caution threshold
(>10%) and their own refusal rule, so they are expected traffic.

**Severity** low (documented-vs-actual gap in the stored results; no printed number is
wrong). **Confidence** certain.

---

## F2 (low) -- `e(perm_sd)` silently switches denominator inside one `permexact` table

**Claim.** Under `permexact`, `e(perm_sd)` is the population SD (divisor *n*) for elements
estimable in all *M* enumerated assignments and the sample SD (divisor *n*-1) for elements
that lost assignments -- in the same run, in the same table, with no marker. The help
describes `e(perm_sd)` only as "SD of the permutation distribution".

**Repro.** `/private/tmp/.../probes/perm-semantics/repro_permsd_denominator.do`
(self-contained; 8 clusters, 2 treated, M = 28, a subgroup level carried by two clusters).

**Observed.**
```
main_att          n=28 M=28  posted=+1.85227db25c93dX-001  == population SD (1/n)
main_SG0_att      n=27 M=28  posted=+1.6a5103d115142X-001  == sample SD (1/(n-1))
main_SG1_att      n=13 M=28  posted=+1.771507ea81bbbX+000  == sample SD (1/(n-1))
```
The n=13 rows are inflated by sqrt(13/12) = 4.1% relative to the convention the n=28 rows
in the same table use.

**Expected.** One convention per table, or a documented switch. The switch is also
inconsistent with the reasoning the project adopted next to it: BS-23 (author ruling
2026-09-07) deliberately moved the SHARE convention off completeness and onto MEMBERSHIP,
on the ground that an element whose reference set contains the observed statistic is
exactly enumerated over the assignments that completed -- and this run confirms the
p-values follow that rule (2/27 and 2/13 exactly, `e(perm_p_se)` withheld,
`e(perm_incomplete)` naming the elements, `e(perm_addone)` empty). The SD beside them was
left on the completeness flag (`wooldid_permtail`: `out[6] = (exact ? population :
sample)`, where `exact` is completeness), so the same reference set is treated as a
population for its p and as a sample for its SD.

**Consequence.** Cosmetic-to-low: `e(perm_sd)` is descriptive and feeds no p-value or
interval (the step-down computes its own studentization). It does make SDs
non-comparable across rows of one table.

**Realistic.** Yes, but narrow: `permexact` plus an element that some enumerated
assignments cannot estimate (a small subgroup, a horizon a reassignment removes).

**Severity** low. **Confidence** certain on the behaviour; the code comment shows the
completeness rule is deliberate, so this is a consistency/documentation finding, not a
miscomputation.

---

## NEGATIVE RESULTS (what I attacked and what would have shown up)

See the table in NOTES.md. Highlights, all CLEAN and all bit-exact unless stated:

- **Every convention's p.** default equal-tailed `min(1, 2 min(P_lo,P_hi))`, `permabsp`,
  `permrawp`, `permside(left)`, `permside(right)`, on the same design and seed. Had any
  convention been mis-wired (a swapped tail, a missing add-one, a wrong denominator) the
  exact integer-count comparison would have caught it on the first column.
- **permside direction.** With a true positive effect, `permside(left)` returns 1.000 on
  the post horizons and `permside(right)` 0.010; `e(perm_p_lo)`/`e(perm_p_hi)` are
  identical across all five convention runs and the default p is exactly
  `min(1, 2 min)` of them. A sign error was the thing to find here; there is none.
- **p = 0 hunt.** Zeros appeared only under `permrawp` (ruled design). The route that
  could produce one elsewhere -- plain shares under a complete enumeration whose observed
  replication does NOT carry the observed statistic -- was attacked directly: in all four
  `permexact` regimes (plain, `permabsp`, `permanalyticse`, `permtstat`) exactly one of
  the 28 enumerated replications carries the observed assignment and its coefficients are
  BITWISE equal to the rep-0 comparator, so membership holds and the floor is real.
- **permtci.** Exactly `b_obs +/- q_level(|t_r|) * se_obs` recomputed from the dump
  (agreement to 1 ulp), and it coincides with the `permabsp` table interval to 3e-16.
- **permstepdown.** Post and pre families recomputed from the dump: `t_obs`, `p_raw`,
  `p_singlestep`, `p_stepdown` all exact; monotone in the |t_obs| order; and the
  invariant `p_raw <= p_stepdown <= p_singlestep` holds for every horizon.
- **Joint / equality / adversarial.** All eight test p-values reconstruct exactly as
  `(1 + #{W_r >= W_obs})/(n+1)` from `e(perm_jointW)/e(perm_sgW)/e(perm_advW)` and their
  `*obs` twins.
- **Pairwise subgroup contrasts.** Identified the coefficient pair by its comparator
  difference (no name matching) and reproduced `e(perm_pwsgs_kernel)` exactly.
- **Are the draws what they claim to be?** The strongest check I ran: rebuild `ttre` from
  the dump's `onset_` columns and refit wooldid with no `permute()`. Five replications
  (including rep 0) came back BIT-IDENTICAL, on the ordinary design and again under
  `permonsetcohorts` with a hand-rebuilt `i`.
- **permonsetcohorts** (§7, shipped 2026-09-16, assessment deferred): conformant, and its
  rebuilt-`i` draws are bit-identical to a hand rebuild.
- **permstrata**: `e(perm_npossible)` = 90^4 = 65,610,000 sampled and 36 exact, both
  exactly right; zero within-stratum composition violations in 100 sampled and 36
  enumerated assignments; the 36 enumerated assignments are distinct.
- **Failure accounting**: the 19 replications that drew zero treated clusters are exactly
  the 19 all-missing rows and exactly `e(perm_failed)`; `e(perm_n)` is the completed count
  and every p uses it.
- **Element-level holes**: a subgroup level on 3 of 24 clusters gives n = 193 vs 199, and
  every kernel row is conformant on that element's own n.
- **Fast-mode estout handoff**: `e(p)`, `e(ci_l)`, `e(ci_u)` are byte-equal to
  `e(perm_p)`/`e(perm_ci)`, and `e(wdidmainresults)` columns 4-6 are byte-equal to the
  kernel, on every run.

- **`permsave()` reconstructability, the strong form.** Under `permtstat` I rebuilt every
  reported `e(perm_p)` from the FILE ALONE (its `_rep = 0` row and its `se_` columns) with
  zero mismatches, and `e(perm_ase_obs)` is bitwise equal to `e(perm_comparator_se)`, so
  the t-ratio the file implies is the t-ratio the engine ranked.
- **`e(gteffects)` under `permute()`**: columns `se`, `p`, `lb_95ci`, `ub_95ci` are missing
  on all 80 rows and `e(perm_gt_blank)` is 1, exactly as documented; `estimate`, `relyr`,
  `nobs` and `w_att` are intact.
- **Failure-driven denominators**: `e(perm_n)` is the completed count everywhere, and the
  Monte-Carlo `e(perm_p_se)`/`e(perm_p_ci)` use `cc/n` on that same n.

- **Sampled shuffle uniformity.** 1400 draws on a design with exactly M = 28 assignments:
  all 28 appeared, chi-square 21.24 on 27 df, p = 0.77. The whole exactness argument rests
  on the sampled draws being uniform over the assignment set; they are.
- **Continuous treatment.** `contreat()` + `permute(shuffle)`: the `ame` rows
  (`main_ame`, `pretreatment_ame`, `ES_ame_*`) are permuted like the `att` rows and every
  kernel row reconstructs exactly. `lattice(all 3)` (assignment-invariant) gives cfx cell
  tables with p and interval populated; `lattice(baseTreated 3)` (assignment-dependent)
  blanks columns 3-6, posts `e(perm_cfx_blank) = 1` and writes no `permdist_cfx*` twins --
  exactly as documented.
- **The "fewer than half" caution fires under permexact**: with four SG1 elements at
  n = 13 of 28, the run printed "Caution: 4 reported coefficient(s) were present in fewer
  than half of the completed replications".
