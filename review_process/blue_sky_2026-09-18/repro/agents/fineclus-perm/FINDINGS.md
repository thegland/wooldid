# FINDINGS -- permute() under finer-than-i clustering, and permonsetcohorts
Agent: fineclus-perm.  Method: DIFFERENTIAL + IDENTITY probing against hand-computed
permutations.  Build: .../bluesky/build/wooldid.ado (asserted with `which wooldid`).

--------------------------------------------------------------------------------
## FP-1 (HIGH, certain, realistic, SEEDED)
### permonsetcohorts: a covariate keyed to the OBSERVED cohort is inert for the point
### estimate but rewrites the permutation null distribution

**Claim.** Under `permonsetcohorts`, `wooldid_permrun` substitutes the rebuilt `i`
into the VARLIST only (ado:7639, `local repcall `y' `repi' `t''`), while `repopts`
(ado:125-163) still names the ORIGINAL variables in `fe()`, `controls()`,
`coarsecohortcontrols()`, `ccc_absorb()`, `interactivecontrols()` and
`customweightmultiplier()`. A covariate that is a function of the observed onset -
and under `permonsetcohorts` the user has ASSERTED that `i` is the onset grouping, so
any adoption-timing covariate is one - is therefore perfectly collinear with the `i`
fixed effects in the OBSERVED fit (silently dropped, zero effect on the estimate) and
NOT collinear inside a replication, where it is live and soaks up exactly the
between-observed-cohort variation the randomization test is generating. The observed
statistic and the reference statistics are then not the same function of the data,
which is the one thing a randomization test may not do.

**Measured** (design B: 6 states, ragged windows and caseloads; i = coh, the onset
grouping; cluster(state) strictly finer; `permexact`, M = 90 enumerated assignments;
`cohx = 10*coh + .5*coh^2`, `wave = (coh<=2)`):

| arm                                | observed b (%21x)      | p      | null SD |
|------------------------------------|------------------------|--------|---------|
| baseline                           | +1.962577d3c317fX+000  | .04444 | 0.7465  |
| `controls(cohx)`                   | +1.962577d3c317fX+000  | .02222 | 0.6627  |
| `fe(wave)`                         | +1.962577d3c317fX+000  | .02222 | 0.7295  |
| `coarsecohortcontrols(c.cohx)`     | +1.962577d3c317fX+000  | .68889 | 2.6027  |

The realistic, dangerous arms are `controls(cohx)` and `fe(wave)`: they are
ANTI-conservative, they cross the 5% line, the observed estimate is bitwise unchanged,
and the user gets no feedback at all because reghdfe drops the collinear term under
`qui`. The `coarsecohortcontrols` row is the same mechanism but its null SD of 2.60
against a baseline 0.75 looks like near-singularity in the replication fits rather than a
coherent alternative model, so it should not carry the headline.

Every observed estimate is BITWISE identical to the baseline. The p-value moves by a
factor of 2 (anti-conservative, crossing 5%) and by a factor of **15.5**
(`coarsecohortcontrols`, conservative). On a second panel the same `controls()` arm
moved p from .11111 to .04444 with the null SD shrinking .7713 -> .5888.

**Control arm.** With `cluster(coh)`, i.e. cluster AT `i` so no rebuild happens, the
same inert option leaves the p-value bit-identical (.33333333 both ways). It is the
rebuild, not the option.

**Mechanism, hand-computed.** Build one replication exactly as the ado does (states
1,3 -> onset 4; 2,5 -> onset 6; 4,6 never; `repi` = the drawn onset's observed cohort
code). Then `controls(wave)` moves that replication's estimate by .0218, while the
same covariate cannot move the observed estimate at all.

**No signal to the user.** Nothing is printed. The bootstrap's structurally identical
re-identification leg DOES warn (ado:403-414, the A10 caution "a covariate
specification names the i variable ...", documented at sthlp:532); the permute rebuild
leg has no equivalent, and `fe(coh)` - literally naming `i` - is accepted silently
(p .04444 -> .02222, observed b unchanged to 3 ulps).

**Repro.** `repro_F1_cohortkeyed_covariate.do` (self-contained, generates its own data).
Arm sweep: `probe15_F1arms.do`.  First measurement: `probe04_repoptsB.do`.

**Anticipating the obvious defense.** "Stata's own `permute` prefix holds everything but
the permuted variable fixed, and the bootstrap's idcluster convention already documents
keeping `i`-keyed terms at their original codes." That defense does not reach this case,
for two reasons. (i) wooldid has ALREADY departed from fix-everything-else here: it
rebuilds `i` inside the replication. Rebuilding `i` while leaving `i`-derived terms fixed
is an inconsistency that neither pure convention has. (ii) The bootstrap's A10 convention
is harmless because re-identification REFINES the partition - each drawn copy of a cluster
is a subset of one original cohort, so a cohort-constant covariate stays inside the span
of the new `i` fixed effects and stays dropped. The permonsetcohorts rebuild REGROUPS: a
rebuilt cohort pools clusters that came from different observed cohorts, so the same
covariate leaves the span and goes live. That is exactly why the help documents the
convention where it does no harm and is silent where it does. (This second point is a
structural argument from ado:9297-9349 and sthlp:532, not a measurement - the confirming
run, `bootstrap(standard)` with and without `fe(wave)`, was not made: the project's
8-instance Stata cap was already breached by other agents.)

**Detectability.** The exposure is mechanically detectable with the check permctl already
runs for the onset: a covariate named in fe()/controls()/interactivecontrols()/
coarsecohortcontrols()/ccc_absorb()/customweightmultiplier() that is constant within `i`
on the design sample (`gegen min/max ... by(i)`). Noted as a detection hook, not a
proposed patch.

**Severity.** High. It is a published p-value that is wrong by a factor of 2-15 with a
bitwise-unchanged point estimate, on an option combination the feature exists to serve.
Realistic: `permonsetcohorts` is FOR the case i = adoption cohort, and an
adoption-wave fixed effect or an adoption-timing control is a natural thing to add.

--------------------------------------------------------------------------------
## FP-2 (MEDIUM, certain, doc-vs-code, SELF-FOUND alongside FP-1)
### The A10 convention is documented and cautioned for the bootstrap re-id leg and is
### neither documented nor cautioned for the permonsetcohorts rebuild

sthlp:532 spells out, for `bootstrap()`, that "A term inside fe(), controls(),
interactivecontrols(), coarsecohortcontrols() or ccc_absorb() that names the i
variable itself keeps its original codes across the copies", and ado:412 prints a
caution when it sees the bare `i` variable in one of those options. The
`permonsetcohorts` rebuild does the same kind of re-identification, has the same
exposure (worse, because the term can be inert in the observed fit and live in the
replications), and the help's `permonsetcohorts` text (sthlp:68, 255, 431, 1240, 1352)
says nothing about covariates at all. Repro: the `fe(coh)` arm of `probe12_F1realistic.do`
prints no caution; `probe15_F1arms.do` shows the numeric consequence.

--------------------------------------------------------------------------------
## WHAT CAME BACK CLEAN (negative results, all measured)

1. **The rebuild machinery is exactly faithful.** All 31 rows of a 30-assignment
   enumeration hand-recomputed externally three ways (observed-code map; fresh
   `egen group()` codes; relabelled codes 10/20/30): all agree with the stored
   replication to <= 1e-15 ABSOLUTE (reghdfe convergence noise) and with each other.
   (One rep differed by 114 ulps in relative terms, but its value was 3.5e-3, so that is
   ~9e-17 absolute - the same 1e-16 noise seen everywhere, inflated by a small denominator.)
   No dependence on the `i` CODES. `probe02_identity.do`.
2. **Bitwise identity under `permtstat`, `aweight`s and `customweightmultiplier()`.**
   Hand-recomputed replications reproduce the stored `main_att` AND `se_main_att`
   bit for bit. `probe09_tstat_weights.do`.
3. **rep-0 comparator.** `e(perm_comparator_b) == e(b)` bitwise (K=12, all columns) in
   fast mode and under `permtstat`; 4.4e-16 under `permanalyticse`, which is the
   documented F9 contract. `probe13_comparator_misc.do`.
4. **The randomization space is exactly right.** Hand enumeration: design A gives 30
   distinct assignments (hand count 5!/(2!1!2!) = 30), design B gives 90
   (6!/(2!2!2!) = 90); every drawn assignment is a permutation of the observed onset
   multiset; the observed assignment appears exactly once. `probe14_enumeration.do`.
5. **p-convention.** Hand-computed two-sided equal-tailed exact share reproduces
   `e(perm_p_att)` (2*2/30 = .1333 on design A).
6. **Invariance.** Under `permexact` the sorted reference distribution is invariant to
   row order, to relabelling `i` (including a permutation that moves the never-treated
   code) and to relabelling cluster ids, to <= 3.6e-15. Under sampled `permute` at a
   fixed `rseed`, row order and `i` relabelling reproduce the draw sequence exactly;
   relabelling CLUSTER ids changes it (the map frame is sorted by cluster value) -
   expected, not filed. `probe06_invariance.do`.
7. **permstrata() x permonsetcohorts composes.** Hand count 3!x3! = 36 matches
   `e(perm_npossible)`; 0 strata violations across all 36 drawn assignments.
   `probe08_strata_probs_guards.do`.
8. **Every guard fires.** i splitting the never-treated -> r459; onset varying within i
   -> r459; cluster finer than i without the option -> r459; `subgroups()` ->
   r198; cross-cutting cluster -> r459; `cluster()` at `i` with the option ->
   `e(perm_irebuilt)=0` and a bit-identical run to not giving the option at all.
   `probe08_strata_probs_guards.do`, `probe10_contreat_misc.do`.
9. **The `wdp_ncode` sentinel claim holds.** With a treated cluster's `cluster()` value
   set missing, the observed run, the permute run and a hand-built rep-0 all return the
   same `e(N)` (70) and the same bitwise `main_att`. `probe13_comparator_misc.do`.
10. **Attrition.** A cluster that draws a late onset and has no post-treatment rows does
    not fail the replication; it simply contributes no treated cells. 60/60 completed.
    `probe07_estimands.do`.
11. **The "relabelling i leaves colfullnames identical" claim holds** on a run with four
    estimands and an event-study window. `probe07_estimands.do`.
12. **permute(probabilities) failure censoring is NOT rebuild-specific**: 12/150 failures
    with the rebuild vs 42/150 without it on the same data, the failures being draws with
    no never-treated cluster. Pre-existing accounting design (the >10% caution fired on
    the 42 case). `probe17_probsfail.do`.
