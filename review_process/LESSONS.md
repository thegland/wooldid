# Lessons

What survives from the probe logs, diagnosis threads and campaign scratch that this directory no
longer keeps. Each entry is a fact about the language, the estimator, the harness or the process,
written so it can be acted on without knowing where it came from. Nothing but the last section is
a status claim.

To read an original, find its removing commit and print the file:

```sh
git log --diff-filter=D --oneline -- review_process/
git show <commit>^:review_process/<path>
```

## Stata and Mata

**A `mata:` block nested inside a Stata control structure aborts the do-file.** Inside a `foreach
{ }` loop it stops with `--Break--` r(1); inside a `program define` it stops with "command end is
unrecognized" r(199). Keep every mata block at top level; if a loop needs mata, restructure so the
mata block does the looping.

**A missing Mata function makes Stata re-read the defining ado from the adopath and re-run its
mata block. A program loaded with `do "<path>/wooldid.ado"` has no such association,** so the same
missing function is simply fatal (r(3499)). This splits the evidence base: a failure that
reproduces under `do`-loading may not exist for an adopath-installed user, and vice versa. Every
suite, harness and agent brief here loads by `do`, so that is the arm the evidence describes
unless it says otherwise. Corollaries measured alongside it: functions in a mata library survive
`clear mata` and reload on demand; re-running the ado from inside one of its own programs fails
r(110), because a program cannot be dropped or redefined while it is executing; re-sourcing a
mata-only file from inside a running program does work.

**Never write `local x = e(macroname)`.** The `=` form evaluates its right-hand side as an
expression, and on `e(cmdline)` that aborts the do-file with r(111) ("wooldid not found"). Carry
macro values as strings.

**Do not rely on `e(scalars)`/`e(macros)`/`e(matrices)` to enumerate an estimation surface.** Two
harnesses in this project observed opposite behaviour: one found all three empty for every
estimation command, `regress` included, with Mata `st_dir("e()", ...)` erroring r(3300); another
found them enumerating even `ereturn hidden` results. A dumper that has to be exhaustive builds
its member list explicitly, grepped out of the pinned ado.

**An out-of-range matrix subscript in expression context returns missing with rc 0.** Reading
`r(p)[1,2]` from a 1x1 matrix does not error and nothing ever complains — the value just goes
missing. The shape that makes this bite: `margins` overjoint results differ by treatment type, the
continuous (dydx) `r(p)` being 1x1 with the real p at [1,1] and the binary one 1x2 with the real p
at [1,2]. Read the dimension; never hardcode the subscript.

**Compare in scalars, not locals.** A numeric value stored in a local round-trips through a
decimal representation and loses precision; two values that are bit-identical can compare unequal.

**An empty matrix literal `matrix X = ()` raises a bare "invalid syntax" r(198)** that names
nothing, and an empty `subpop()` raises r(461). Both are what an upstream sample wipe looks like
by the time it surfaces, which is why neither message points at the real cause.

**A `float` holds integers exactly only to 2^24 = 16,777,216, and Stata's default `c(type)` is
float.** A variable that comes back from a collapse, a merge or a bare `gen` can therefore be
float without anyone asking for it, and any comparison of it against a `long` row counter — the
usual way of keeping one chosen row per cell — starts matching the wrong rows once the sample
passes that size. Declare the type on both sides of such a comparison.

**A Stata launched from a bash call that returns is killed with `--Break--`** when the harness
reaps the process group. Keep the calling shell alive until Stata exits.

**Environment, on this machine:** invoke Stata by the direct binary
(`/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp`) because the `stata-mp` PATH symlink is
broken; paths handed to do-files must contain no spaces, since they pass through Stata locals; and
do not start a leg while a full assessment run is live.

## The estimator and its degeneracies

**The whole AME-SE family follows from one FWL identity and two channels.** With cell coefficients
`gamma_gh = s_gh - phi_block` and any reported ame an average `b = sum w_gh gamma_gh`, the true
variance has two parts: within-cell slope noise and baseline noise.

- *Channel A, CRVE cell-score annihilation* (any `cluster()` at or above the cell level): the OLS
  normal equations zero the within-cell residual moments, so the plug-in CRVE score of every cell
  coefficient vanishes and the entire within-cell component is estimated as exactly 0. No
  cohort-level plug-in CRVE variant recovers it — CR2/CR3 divide by a singular `I - H_gg` there.
- *Channel B, plug-in leverage bias* (robust, or clustering strictly finer than the cells): the
  meat sees within-cell residuals, but they are shrunk by fitting two within-cell parameters on a
  handful of observations (six, in the design that measured it), so the within-cell component
  comes in at roughly 0.55-0.65 of truth and SE ratios land at 0.74-0.81. Textbook HC0/HC1
  small-cell bias, made structural by a specification that guarantees tiny saturated cells.

Reported/true SE ratio is approximately `sqrt[(V_base + k*V_slope)/(V_base + V_slope)]` with k = 0
under channel A and 0.55-0.65 under channel B, and `V_slope` falls like 1/(number of cells
averaged) — which is why event-study horizons sit near 0.2 while the pooled ame sits near 0.5.
Before hypothesising a code bug from an unexplained SE ratio on a saturated-cell specification,
check it against this formula. The related standing pathology: V_es-based joint tests are severely
oversized when the cluster VCE is blind to within-cell variation.

**A degenerate variance can present as full rank with plausibly varying entries, so screen on the
design, not on the matrix.** Rank, eigenvalue and condition-number guards are spoofable by exactly
the configurations that most need them — a period-varying baseline block reproduces the
annihilation while destroying the rank signature. The predicate that works is structural (which
cells are contained in which clusters); numeric guards are known-spoofable and are not used.

**Quantify a containment predicate with ANY, never EVERY.** Partial containment reproduces most of
the pathology while defeating a universally quantified guard. Measured: with a cluster variable
that is the cohort for some units and the unit for others (80 clusters), 10 of 15 treated cohorts
had every cell inside one cluster and 5 did not; a drafted "every treated cell inside a single
cluster" predicate did not fire; the result was a full-rank V (eigenvalues .0005-.0245),
horizon-varying SEs at 0.37-0.77 (mean about 0.50) of the unit-clustered SEs on the same data,
zero cautions, and p = 1.2e-08 on one horizon of a true null. Real analogues are ordinary: cluster
variables that cross-cut cohorts (region, industry), or multi-way clustering with one containing
dimension.

**Screen ordering is a correctness property, not a presentation detail.** The unidentified-design
refusal has to run before the containment screen, or the containment screen's "use finer
clustering" advice is given on data where it is wrong. The remedy a caution prints is part of its
contract.

**A screen that reads `e(b)` cannot detect a sample-level drop.** The dropped-treated-dummy screen
iterates over `o.`-marked coefficients *present* in `e(b)`; when rows leave `e(sample)` the
treated levels are absent entirely, so the screen has nothing to inspect and the wreck sails on
into the margins phase. The canonical case: `ccc_absorb(x2)` expands to a separate x2 slope per
treated cohort-by-relative-period cell, and with unit-level `i` every such cell is a single
observation, so reghdfe drops them all as singleton slope groups and the entire treated sample
vanishes. That specification is unestimable by construction — not an estimator bug — and its
symptom downstream is an empty matrix literal or an empty `subpop()`.

**A posted covariance can be exactly rank 1, and what the user then sees is a floating-point
accident.** Under cohort-level `cluster(i)` with `contreat`, every entry of the posted
per-subgroup AME covariance is identical, so a cross-subgroup equality contrast has exactly zero
estimated variance. Depending on how rank detection lands, the same contrast prints `.` (full
drop), or garbage near p = 5e-80 (partial drop), or garbage near p = 1e-109 (no drop). Three
different-looking outputs, one degenerate input.

**When a VCE-dependent rank deficiency looks like a bug, run the iid arm as the discriminator.** A
robust or clustered VCE that zeroes treated residuals manufactures an exact rank deficiency: rank
9 of 10, smallest eigenvalue around 1e-17, a variance ratio of exactly 4 between two pre-period
horizons, and a correlation of 1.0 to machine precision. iid on the same data keeps full rank, a
ratio of 2.22 and a correlation of 0.45. If iid is clean, the deficiency is a property of the VCE
and not of the design — and coarsening `i` to the onset group with unit fixed effects restores
full rank (ratio 2.26, correlation 0.41).

**Monte-Carlo p-values: add one where a statistic is ranked, and nowhere else.** A run that ranks
an observed statistic among recentred null draws needs the add-one; a routine that inverts a
percentile interval by locating 0 among raw draws ranks nothing and has nothing to add it to. The
size cost of getting it wrong is exactly 1/(n+1) on the anti-conservative side: at 20 usable
replications a no-add-one joint p is a 9.5% test at a nominal 5%, shrinking to .001 by n = 1000.

**99, 199 and 999 are the textbook replication counts because `alpha*(R+1)` is an integer there**,
so the add-one and no-add-one conventions have identical rejection sets and identical exact size.
That identity is a knife edge the user does not control: a joint test runs on listwise-complete
replications, and at a request of 99 the usable count came in at 84-89 in 8 of 8 runs, never the
aligned 99. Under add-one the test stays valid whatever the count turns out to be; without it,
validity depends on a number the user cannot see and did not choose. At these replication counts
the `p <= alpha` versus `p < alpha` reading is the larger hazard of the two.

**The bootstrap engine has no rep-0 cell comparator.** It mints its cell accumulator keys from the
rep-0 run and discards that run's cell values; only the permutation engine keeps them. Cell draws
are not retained in `e()` either, so a per-cell recomputation against the rep-0 cell values is not
available on the bootstrap path without a code change.

**`safety(off)` is not uniformly a downgrade, and the two shapes are easy to confuse.** Some
screens have an explicit caution branch and genuinely downgrade — the
`contreatwithin`-with-`byTime` refusal prints a caution instead of an error. The estimability
screens wrap the whole block in `if "safety" != "off"`, so under `safety(off)` nothing prints at
all: the check is skipped, not downgraded. Before documenting a screen as "downgraded to a
caution", check which shape it has.

**A package that aggregates owns the statistical properties of the aggregate.** `regress y i.g,
vce(cluster g)` posts near-zero SEs on cluster-saturated dummies and StataCorp ships it, but there
the pathological numbers sit on raw incidental coefficients the user can see, next to the factor
notation that produced them. A command that aggregates those coefficients into a headline
estimate, prints it in a labelled SE/t/p/CI column and returns it for downstream consumers cannot
fall back on "the regression did it". That is the standing argument behind every refusal and
caution here.

## Probes, harnesses, and verifying a change

**A null result from an instrument that could not have detected the effect is not a null result.**
A differencing harness reported that a guard made things slower; direct timers put the affected
blocks at a stable 5.3% cost. The harness was wrong because the run-to-run variance of a single
bootstrap configuration is about 17%, several times the effect — it could not have detected it at
any sign. Compute the instrument's resolution before believing its null.

**Always build the positive control.** A probe reporting "the check never fired" is
indistinguishable from a dead probe. Pair it with a calibration run that makes the same instrument
fire when it should.

**Write the prediction into the do-file before you run it.** It costs nothing and it is what lets
a probe's result count as evidence rather than as a post-hoc reading. State the expected
magnitude, not just the direction.

**Reproduce from the claim alone.** An independent lens builds its own panel and never reads or
runs the prior agent's script; that is what makes agreement between lenses informative and what
makes a refutation rate mean anything.

**Assert the probe site, don't trust it.** A script that instruments a build by line number should
`assert` the line's content and stop if the build does not match, so a moved target fails loudly
instead of silently instrumenting the wrong line.

**The load-bearing check on a harness fix is a planted decoy, not a green run.** A suite runner
once captured each configuration's saved plots with a directory-wide `svp*` glob while all 42
shards shared one working directory, so one configuration could record another's files as its own
and erase them before the owner dumped. The dangerous direction is the erasure: the victim then
records an empty file list, and a genuine regression compares clean forever with no symptom. The
same hardcoded pattern was simultaneously too broad and too narrow — a configuration whose save
prefix did not match `svp*` was never checked at all, on either side, silently. The fix is a
private working directory per shard plus parsing each configuration's own prefix; the proof is a
decoy file planted in each old shared directory, because a merely-broken capture would also have
gone green.

**Store the recipe, not the dumps.** Suite outputs are deterministic functions of tracked inputs —
the two builds, the suite definition, and the written seeds — so copying 15 MB of artefacts only
stores a snapshot whose provenance decays, while the regeneration command does not. The corollary
is a rule about claims: a "bit-identical" result tied to dumps that have since been overwritten is
a claim owed a re-tie, not a claim proved.

**Old-versus-new suite runs stage git HEAD, never the working tree.** Uncommitted changes are
invisible to the comparison.

**Relative paths in a test are a source of phantom failures.** A test's `local REPO "../.."`,
resolved against whatever directory Stata happened to start in, once produced a spurious 21-of-22
certification result. Resolve from an absolute root.

**Design a probe so the bug can actually fire.** When probing a frame join for a reorder defect,
choose labels whose sort order deliberately crosses the cluster order; with cooperating orders the
defect is invisible no matter how many replications you run.

**A probe corpus is reproducible only if its generator inputs live in the repository.** Two
corpora here failed that test in different ways: in one, 1,089 of 1,168 scripts hardcoded an
`adopath` into a temporary build directory that no longer exists; in another, every script sourced
a variable-generating include that was never copied into the repo at all. Neither runs today. The
single probe fixture that remains regenerable is the one whose generator sources an in-repo
include, sets an explicit seed, and needs nothing external — keep a 4 KB generator over a 400 KB
log every time.

**Test the type, not the scale.** A failure mode that only appears past 17 million rows can be
measured on two hundred: generate counter values just above 2^24, check the storage type the
operation actually returns, and check how the comparison behaves there. Building the
17-million-row panel is the expensive way to learn the same thing, and it is usually the reason
the check never gets run.

**Instrumented builds are recorded as a diff against a pinned blob**, with the run commands, seeds
and data-generating formula written into a README. That pair reconstructs the run; an ad-hoc
do-file pointing at an ephemeral build directory does not.

## How a finding is adjudicated here

**A status line is only true on the day it is written.** "Nothing in it has been decided", "OPEN,
ruled, NOT started", "Nothing below is implemented yet" were all accurate when typed and badly
wrong within days. If a status claim is going into a document that will outlive the week, date it
and say what would make it false.

**A record nobody maintains is worse than no record**, because its stale lines get read as
current. An audit of fourteen process documents turned up 191 stale claims, one of which caused
thirteen already-ruled decisions to be reported to the author as still open.

**Never take a document's own header as evidence of state.** Read the code, or the ledger the code
cites. "All closed" at the top of a file is a status line like any other.

**When records are removed, leave a table naming what went and where the live answer is now**,
plus the `git log --diff-filter=D` / `git show <commit>^:<file>` recipe. Surviving citations to a
removed file are not broken — they name something git still holds — and rewriting forty of them by
hand invites exactly the unverified edits the cleanup was meant to stop.

**Never key a join on identifiers the producers chose.** One campaign's 223 findings carried
agent-assigned ids with nine collisions (the same id used three times in three places), and every
downstream report joined on them had to be regenerated against a unique key.

**A predicate copied from an older remedy menu has to be fixed at its source too**, or the next
draft re-imports the wrong wording.

**Label the strength of your own evidence.** "An arithmetic bound plus observation, not a bitwise
proof" is the right form. And a consistency wart that has never moved a reported number is a wart,
not a defect — say which one you have.

**Keep the reasons refusals were refused, with the measurement attached.** A rejection record with
numbers in it (for instance: a degeneracy guard's threshold sits between a relative dispersion
floor of about 1.6e-15 and the smallest realistic dispersion of about 1e-7) is what stops a future
campaign re-litigating dozens of already-refused tickets.

**The code carries mechanism; the review record carries provenance.** When the help file points a
user at "the design record", something has to be behind the pointer — and when a shipped screen's
scope is marked load-bearing, the reason it is scoped that way (for instance, that the predicate
is true on the ordinary binary configuration, so the mirror-image screen would fire on nearly
every standard run) lives only where it was written down.

**Re-run the implementation audit even when work has been declared done.** The audit shape that
compares every ruled item against the shipped code found twelve items quietly unimplemented after
they had been declared complete twice.

## Open items recorded here because they are recorded nowhere else

Not lessons; unclosed observations that would otherwise be lost with the files that held them.
Each carries the condition that retires it.

1. **`review_process/blue_sky_2026-09-18/findings_raw.json` is the pre-correction copy.** It has
   no unique-id field and repeats nine ids, so a join keyed on `id` reproduces the defect the
   published reports were regenerated to remove. It is kept because it holds the only untruncated
   claim, expected, observed and refuter text. *Retire when the file carries a header saying so,
   or is replaced by the unique-keyed dataset.*

Two further items were found and closed rather than recorded: the help's permutation paragraph
correctly required `permreps()` of at least 20 and then closed by saying neither engine has a
hard minimum, which the ado contradicts at its own refusal - the trailing clause is gone; and the
evidence README asserting that nothing in its directory was tracked by git, when all twenty files
were, went with that directory.
