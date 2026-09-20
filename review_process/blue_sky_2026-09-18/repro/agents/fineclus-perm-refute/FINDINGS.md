# SKEPTIC REPORT -- refutation pass on fineclus-perm FP-1 and FP-2
Agent: fineclus-perm-refute.  Build asserted with `which wooldid` in every probe:
`/private/tmp/claude-501/.../scratchpad/bluesky/build/wooldid.ado`.
Method: rebuilt every claim from its own terms on MY OWN data (7 ragged administrative
states, i = onset cohort, cluster(state) strictly finer, M = 210 enumerated assignments
so the p grid is 2/210 = .00952, four times finer than the prober's M = 90).
I did not run the prober's do-files.

================================================================================
## VERDICT ON FP-1 -- **PARTIAL** (mechanism CONFIRMED and it generalises;
## the stated magnitude is WRONG IN BOTH DIRECTIONS)

### What I confirmed, independently

1. **The premise holds.** A covariate constant within the observed `i` is inert for the
   reported point estimate. On my panel, `%21x` of `e(b)[main_att]`:
   ```
   no covariate                    +1.6690b859bb834X-002
   controls(cohbase)               +1.6690b859bb834X-002   BITWISE IDENTICAL
   coarsecohortcontrols(c.cohbase) +1.6690b859bb834X-002   BITWISE IDENTICAL
   fe(wavegrp)                     +1.6690b859bb838X-002   4 ulps (reghdfe noise)
   fe(coh)                         +1.6690b859bb82cX-002   8 ulps
   ```
2. **The reference distribution really is rewritten, and only under the rebuild.**
   Paired assignment by assignment (the enumeration is identical across arms -- I checked
   the `permsave` `onset_*` blocks: 0 mismatching rows of 211), `controls(cohbase)`
   changes **208 of the 210** reference statistics.  The two it does not change include
   the identity assignment, where the difference is **bitwise zero**
   (`+0.0000000000000X-3ff`).  That is exactly FP-1's mechanism: observed statistic and
   reference statistics are not the same function of the data.
3. **The control arm behaves.** With `cluster(coh) == i` (no rebuild, `e(perm_irebuilt)=0`)
   the same covariate leaves p and the null SD **bit-identical** (p = .33333333,
   SD = .2566696171, with and without).
4. **Nothing is printed.** Verified NON-QUIETLY for all three arms
   (`probe07_cautions.log`, `probe10_cccnoisy.log`):
   `permute(shuffle) permonsetcohorts fe(coh)` -- the covariate literally naming `i` --
   prints no covariate caution of any kind; nor does `controls(cohbase)`; and the
   `coarsecohortcontrols(c.cohbase)` arm prints EXACTLY the same three cautions
   (permreps grid, few-treated-clusters, and the permute footer) as the same run without
   the option -- 3 cautions in each, none of them about the covariate.  The no-permute
   CCC run prints an entirely ordinary table with no caution at all and the same
   estimate .35016144.  So nothing in the output distinguishes the arm whose reference
   distribution has been replaced.
5. **It is the inert-in-observed / live-in-replication asymmetry, not the option.**
   `coarsecohortcontrols(c.a_s)` with a STATE-level control (varies within the cohort, so
   it is live in the observed fit too) gives a normal reference distribution
   (p = .4762, null SD = .2169) on the same panel where the cohort-level control gives
   SD = 1.2655.

### The corrections

**(a) "p moves 2x-15x" does not reproduce for the arms FP-1 headlines.**
On two independent panels, with a p grid four times finer than the prober's:

| panel | arm | p | null SD | SD vs BASE |
|---|---|---|---|---|
| A | BASE | .02857143 | .22768479 | -- |
| A | `controls(cohbase)` | .02857143 | .22076719 | **-3.0%** |
| A | `fe(wavegrp)` | .02857143 | .21969891 | -3.5% |
| A | `fe(coh)` | .03809524 | .21818330 | -4.2% |
| A | `coarsecohortcontrols(c.cohbase)` | **.89523810** | **1.26552397** | **+456%** |
| B | BASE | .10476190 | .24825485 | -- |
| B | `controls(cohbase)` | .11428571 | .22942360 | -7.6% |
| B | `fe(wavegrp)` | .08571429 | .23534421 | -5.2% |
| B | `fe(coh)` | .11428571 | .22928515 | -7.6% |
| B | `coarsecohortcontrols(c.cohbase)` | **.80000000** | **1.37229962** | **+453%** |

In the `controls()`/`fe()` arms the p-value moves by **at most one grid step**
(2/210 = .00952) and **in both directions** -- panel B's `controls()` arm made the p
LARGER.  The prober's "p .04444 -> .02222, a factor of 2" is one step on their M = 90
grid (2/90 = .0222), i.e. **one reference statistic crossing the observed value**, not a
doubling of anything.  Their null-SD shrink of 11-24% does not reproduce either: I get
3.0-7.6%.  "Anti-conservative, crosses 5%" is therefore a statement about their
particular panel, not a property of the mechanism.

**(b) The arm FP-1 explicitly demoted is the one that matters.**  FP-1 says the
`coarsecohortcontrols` row "should not carry the headline" because "its null SD of 2.60
against a baseline 0.75 looks like near-singularity in the replication fits rather than a
coherent alternative model."  I reproduce that arm on both my panels, and it is the
largest, most consistent and most consequential effect in the whole probe: the null SD
goes up **5.5x on both panels**, and the reported p goes from **.0286 to .8952**
(31x) and from **.1048 to .8000** (7.6x) while the point estimate is **bitwise
unchanged**.  `e(perm_failed) = 0` and no caution is printed, so "near-singularity"
is invisible to the user: they see a normal-looking estimate and a p of 0.90.
Why it is the worst case: `coarsecohortcontrols()` adds the control AND its interaction
with every Tz_it dummy (sthlp:1061).  With a cohort-constant control every one of those
terms is collinear with the Tz_it dummy itself, so the whole block drops in the observed
fit; inside a replication the rebuilt cohort pools states with different values, so the
whole block goes live at once and the replication estimates a far richer model.
`coarsecohortcontrols()` is also the option a user reaches for when they have a
cohort-level control, which is precisely the class that is inert.

**(c) One structural argument in FP-1 is stated as fact and was not measured -- I
measured it, and it is right but narrower than claimed.**  FP-1 says the bootstrap's A10
convention "is harmless because re-identification REFINES the partition" and cites the
caution at ado:403-414 as the analogue that permute lacks.  The caution at ado:410 is a
**name-equality tripwire**: `if "\`wdbcv'" == "\`i'"`.  Measured: a bootstrap run with
`controls(cohbase)` (constant within `i`, not NAMED `i`) prints nothing either.  So the
class FP-1 identifies is uncautioned on **both** engines; only the literal `fe(i)` spelling
is cautioned, and only on the bootstrap's nesting leg.


### (d) A 24-dataset SHARP-NULL sweep with RANDOMISED assignment settles the direction
`probe04_multi.do` (6 states, cohorts drawn at random in every dataset, y with no
treatment term at all, M = 90 enumerated, 24 datasets, arms BASE /
`controls(cohbaseY)` where cohbaseY = the cohort's own mean pre-2005 outcome /
`controls(cohidx)` / `fe(early)` / `coarsecohortcontrols(c.cohbaseY)`).
Ratio of the arm's null SD to BASE's null SD, and the p-values:

| arm | SD ratio med [p25,p75] | max | narrower than BASE | p: min / med |
|---|---|---|---|---|
| BASE | -- | -- | -- | .0667 / .6222 |
| `controls(cohbaseY)` | 0.977 [0.889, 1.016] | 1.29 | 16 of 24 | .1333 / .5444 |
| `controls(cohidx)` | 0.981 [0.945, 1.005] | 1.23 | 16 of 24 | .1333 / .5778 |
| `fe(early)` | 0.991 [0.958, 1.039] | 1.29 | 14 of 24 | .1111 / .5889 |
| `coarsecohortcontrols(c.cohbaseY)` | **3.85** [2.14, 7.56] | **229.6** | **0 of 24** | **.4000** / .8222 |

So for `controls()`/`fe()` the distortion is a **few percent, in both directions** -- a
slight median narrowing (2-3%) that is swamped by dataset-to-dataset noise, not the
systematic anti-conservatism FP-1 describes.  For `coarsecohortcontrols()` the null SD is
**never** narrower, is more than doubled on 19 of 24 datasets and more than quintupled on
10 of 24, p is larger than BASE on 20 of 24, and **across all 24 datasets the smallest
p-value the arm could produce was .40**.  The permutation test in that arm cannot reject
anything, silently, with `e(perm_failed) = 0` throughout.

### Severity
Kept at **high**, but on different grounds than FP-1 states: not "anti-conservative
controls()/fe() crossing 5%" (small, and not signed), but "coarsecohortcontrols() with a
cohort-level control silently converts p = .029 into p = .895 with a bitwise-identical
point estimate, zero failed replications and no caution."  That is brief category 3 (an
inference path that does not deliver what it claims), realised silently.
`realistic` = true: `permonsetcohorts` asserts `i` IS the onset grouping, so a
"cohort characteristic" control is the natural thing to write, and
`coarsecohortcontrols()` is the option built for it.

Line number: FP-1 cites `ado:7639` for the varlist-only substitution; on the pinned
build it is **`ado:7641`** (`local repcall \`y' \`repi' \`t''`).  Everything else in
FP-1's code citations checks out.

Repro: `repro_FP1_verify.do` (self-contained, ~6 min).

================================================================================
## VERDICT ON FP-2 -- **PARTIAL**, severity revised DOWN to **low**

Every factual assertion in FP-2 checks out, verified non-quietly on my own panel
(`repro_FP2_cautions.do`):

| run | caution "names the i variable"? |
|---|---|
| A `permute(shuffle) permonsetcohorts fe(coh)` | **no** |
| B `permute(shuffle) permonsetcohorts controls(cohbase)` | no |
| C `bootstrap(standard) cluster(coh)==i att_i fe(coh)` | **YES** |
| D `bootstrap(standard) cluster(coh)==i att_i controls(cohbase)` | **no** |
| E `bootstrap(standard) cluster(state) att_i fe(coh)` | no |

But it is not a separate defect and it is not a documented-vs-actual divergence:
1. It is FP-1's remedy restated.  Everything FP-2 asks for (a caution, or help text) is
   what fixing FP-1 would produce.
2. Row D shows the bootstrap "convention" is a **name-equality tripwire**, not coverage
   of the class.  So the honest statement is "neither engine says anything about a
   covariate constant within `i`; the bootstrap has a narrow tripwire for the literal
   `fe(i)` spelling", not "documented and cautioned for the bootstrap, silent for permute".
3. FP-2's load-bearing inference -- that sthlp:532 "implicitly tells a reader that
   permute() does not have the issue" -- is inference about a reader, not a promise the
   help makes.  Category 5 needs the help to promise X while the code does Y.
4. `known_bydesign.md` §6 already carries "**No display line for the permonsetcohorts
   rebuild of `i`** -- `e(perm_irebuilt)` is the record.  Ruled 2026-09-16, *Not to be
   reopened*."  Adjacent rather than identical (that ruling is about announcing the
   rebuild, not about covariates), but it is the author declining to print on this path.

================================================================================
## NEW FINDING I FOUND MYSELF -- RF-N1 (medium)
### permute()'s fast-mode default silently disables the one diagnostic wooldid owns for
### "this control has no within-cohort-period variation"

`sthlp:1068` documents that `interactivecontrols()` prints a caution naming the control
when every treated cohort-period slope term drops out for want of within-cell variation,
"(that check does not run under `fastnose` or `safety(off)`, where the reduction is
silent)".  Since 2026-09-05 the observed run under `permute()` **takes the fastnose path
by default** unless `permtstat`/`permanalyticse` (ado:5414-5418).  So:

```
wooldid y coh t onset, att cluster(state) interactivecontrols(cohbase)
  -> "Caution: every treated cohort-period slope term for the interactive control coh..."

wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) ///
        permonsetcohorts interactivecontrols(cohbase)
  -> that caution is ABSENT

... same + permanalyticse
  -> the caution is BACK
```
The user never typed `fastnose`.  The caution is suppressed on exactly the path where a
control with no within-cohort variation is most dangerous: under `permonsetcohorts` it
goes live inside every replication (FP-1).  The two documented facts compose into an
undocumented one, and the composition removes the only warning in the package that names
this exact condition.

Why it bites: `interactivecontrols()` is itself an observed-reduced / replication-full
specification under `permonsetcohorts`.  On my panel the observed run collapses to
`controls(i.t#c.x)` (the caution's own words) while every replication estimates the full
interactive model, and the null SD doubles as a result: `.5175` with the rebuild
(`cluster(state)`) against `.2555` without it (`cluster(coh)`), p `.943` against `.667`.
So the caution is gated off on exactly the path where the model it warns about is the one
being rewritten.  (`interactivecontrols()` is NOT in FP-1's inert class -- it does move the
observed point estimate -- which is why this is a separate finding.)

Repro: `repro_NF1_iccaution.do` (self-contained, fast).

================================================================================
## NEGATIVE RESULTS (things I attacked that came back clean)
- **The enumeration is arm-invariant.** Under `permexact` the assignment sequence is
  byte-identical across covariate arms (0 mismatching rows of 211 on the `permsave`
  `onset_*` block), so nothing in FP-1 is an artifact of different draws.
- **`e(perm_irebuilt)`, `e(perm_failed)`, `e(perm_reps)` are truthful** in every arm I
  ran: 1/0/210 under `cluster(state)`, 0/0/6 under `cluster(coh)`.
- **`subgroups()` under `permonsetcohorts` is already refused** with an r(198) that states
  the FP-1 mechanism in full (ado:6022-6034: "a rebuilt cell pools the clusters that drew
  one onset, and those clusters come from different observed cohorts").  The author has
  already reasoned about this exact asymmetry for one option slot.  So it is not a new
  finding -- but it does show the exposure test FP-1 asks for is one the code already
  performs elsewhere.
- **A state-level control is safe.** `fe(state)`, `coarsecohortcontrols(c.a_s)` -- anything
  that varies within the observed cohort -- is live in the observed fit too, and its
  reference distribution is ordinary.  The defect class is exactly "constant within `i`".
