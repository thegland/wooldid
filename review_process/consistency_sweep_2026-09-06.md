# Internal-consistency sweep: where wooldid's own behaviour switches between features

2026-09-06. Branch `bugfixes`, tree `2d59e70`. Package files read at
`wooldidinstall/wooldid.ado` and `wooldidinstall/wooldid.sthlp` (1,375 lines).

## 0. Purpose and method

The author asked one question: **from a user's perspective, where is wooldid internally inconsistent —
where do default settings or behaviours switch, literally or in spirit, between features, so that someone
who learned how wooldid behaves with feature X is confused when it behaves differently with feature Y, or
with X plus a bonus option?** This report is the answer as a catalogue, not as a verdict. Every item states
what wooldid does in each place side by side, cites the ado line or quotes the help, says why a user
stumbles, and lays out harmonization pathways — typically (a) unify on X, (b) unify on Y, (c) keep both and
document or label — each with its cost and what it would break. Items whose asymmetry the sweep judged
defensible are included on purpose: the author asked for the list, not a ruling, and section 4 collects the
ones the sweep believes are inherent to the methods so they are not raised again.

**Method.** Seven lenses were swept independently — engines (analytic vs `bootstrap()` vs `permute()`),
options (naming, requirement regimes, argument handling), policy (refuse vs caution vs skip, `safety()`,
screens), stored results, display (tables, footers, cautions, plots), feature interplay, and documentation
— each reading the ado and the help end to end and running wooldid against the shipped fixtures where the
display or a return code was the point. Every finding was then handed to a separate adversarial verifier
told to check both halves against the tree, re-running Stata where a claim was about output or a return
code, and to correct or refute rather than confirm.

**Counts.** The sweep produced **132 catalogued items** — 19 engines, 17 options, 17 policy, 20 stored
results, 20 display, 20 interplay, 19 documentation. Of these, **48 carry verifier evidence**: **31
confirmed** and **17 partly confirmed** (a correction to one half, carried into the text below and flagged
on the item's verification line). **0 were refuted outright**; the ~20 sub-claims the verifiers knocked out
of otherwise-standing items are listed in section 5 so nothing is silently dropped. **84 items were not
reached before the shared Stata process cap saturated** and are reported as the finding's author wrote
them, marked **"Not verified (cap)"** — their line numbers were read from source but not re-checked, and no
run was made.

Priority is the finding author's, and where a verifier moved it the item says so
("filed high; verifier: medium"); sections 1 and 3 use the verifier's value. About **thirty-four items are
cross-lens twins** — two lenses found the same fact through different doors. Rather than merge and lose the
count, the later id carries a short pointer body naming the anchor item, the corrections that transfer, and
whatever the second lens added; where both halves were independently verified with different evidence, both
are written in full and cross-referenced.

**Line numbers.** The lenses and their verifiers read the tree at `01a7f7a`/`5b2acf1`, and the citations
below are that tree's. One commit landed afterwards, `b6ef398` (comments and footer text from the F033 and
F018 research records), which shifts every line below 8312 by **+5** and every line below 12169 by a
further **+6**; `wooldid.sthlp` is unchanged. The headline items in section 1 were re-anchored to `2d59e70`.
`b6ef398` also changed one thing this sweep is about, and the affected items say so: the bootstrap footer's
resolution sentence now reads "The rank-based **coefficient** p-values have resolution 2/n ... The joint,
equality and adversarial p-values are add-one shares over their own listwise-complete replication counts
(their smallest value is 1/(n + 1), never 0)" (ado:10113 at `2d59e70`). That closes, on the ado side, the
gap CS-1's verifier found — under `boottstat` nothing printed used to name the convention at all. The help's
half of the same split (CS-121) is still open.

No item here re-opens an author ruling of 2026-09-06; where a ruling is itself one half of an asymmetry
(F030's "no reps floor" against F036's permreps floor, F034's fixed-metric joint tests, F067's on-request
BCa) the item names the ruling and prices the pathway that would reverse it.

## 1. The fifteen most likely to trip a user

| id | theme | feature X does | feature Y does | shortest pathway | priority |
|---|---|---|---|---|---|
| CS-2 | confidence level | analytic tables take `c(level)`, so `set level 90` names the columns `lb_90ci` | every bootstrap/permute rewiring resolves an unspecified `level()` to a literal 95, and the footer then prints "the 95% CIs" under those same 90-named columns | resolve `level` once in the shell before the engine calls, as ado:418 already does for `wooldid_ptp` | high |
| CS-4 | stored results | analytic `e(V)` is `diag(se^2)` and the help says `test`/`lincom` are invalid | under `bootstrap()` `e(V)` is the full bootstrap covariance with `e(df_r)` deliberately aligned so a hand `test` reproduces `bootdisplay(wald)`; under default `permute()` there is no `e(V)` at all | post `e(vsource)`, rewrite the Incompatible Post-Estimation paragraph per engine, and move `e(df_r_analytic)` above the `bootdf()` if/else where it is currently lost | high |
| CS-1 | p-value conventions | bootstrap coefficient rows use plain order-statistic counts and can print exactly 0 | permute's sampled draws and the bootstrap's *own* joint/equality/adversarial tests use the add-one form; `permexact` uses the bootstrap's plain rule | name the convention per block in the help (the footer half shipped in `b6ef398`) | high |
| CS-23 | fast mode | `permute(shuffle) ... adversarial` is legal and reports a permutation adversarial p | adding the redundant word `fastnose` to that same call is r(198) "fastnose is not compatible with adversarial" | move the `fastnose` refusal block after permute's own fast-mode assignment, or exempt permute runs from those clauses | high |
| CS-38 | pre-trend diagnostics | `pretrendpower()` without its prerequisites is refused at parse with a named remedy | bare `adversarial` (no `esprelength`/`espostlength`) clears the option at ado:1042 with no message at all and runs to completion | screen adversarial's prerequisites in `wooldid_validate` as `wooldid_ptpcheck` does, or give ado:1042 a caution | high |
| CS-98 | safety() | under default safety `makecfxplots` with a non-`att` estimand prints "wooldid will not produce these plots unless you specify safety(off)" and produces none | adding `contreatwithin` rewrites `safety` to `off` *after* that caution is printed, so the run prints the sentence and then draws the plots | gate the caution on the post-rewrite value, or move the rewrite above it | high |
| CS-54 | e() naming | `bootdisplay(bc)`/`(bca)`/`(wald)` name their hidden tables `e(boot_bc_mainresults)` etc., so the name derives from `e(boot_display)` | the default `percentile` posts `e(boot_display)` = "percentile" against `e(boot_pctl_mainresults)`, so the derivation fails on the default and only on the default | post `e(boot_family_tag)` as the component to interpolate | high |
| CS-57 | e() naming | `e(perm_ci)` and every `e(boot_ci_*)` use rows `ll`/`ul` | `e(perm_lb_att)`/`e(perm_ub_att)` are read from the kernel's *null band* rows, so `lb` names two different quantities in one run and neither is `e(perm_ci)`'s `ll` | rename the legacy scalars `perm_bandlo_`/`perm_bandhi_`, keeping the old names as aliases | high (filed high; verifier: medium-high) |
| CS-58 | e() shapes | `e(bootdist)` is the full coefficient-draw matrix, stripe-aligned to `e(b)` | `e(permdist)` is *not* its twin — one column per requested estimand, on the `e(perm_stat)` scale (the t-ratio under `permtstat`); the twin is `e(permdist_b)` | deprecate `e(permdist)` in favour of `e(permdist_b)`, or rename it `e(permdist_est)` | high |
| CS-60 | e() types | `e(joint_horiz)` and both adversarial p's are `ereturn scalar` | the per-subgroup joint tests `e(joint_horiz_s#)` are `ereturn local`, so one loop over subgroups gets numbers for one test family and strings for the other | post the subgroup joint p as a scalar at all three posting sites, and fix the display block that dereferences it as a macro | high |
| CS-74 | column labels | `permute()` rewrites `se`→`se_analytic`, `t`→`t_perm`, `p`→`p_perm` | `bootstrap()` renames only the statistic (`t_wald`/`z_wald`/`t_pt`); the column literally named `p` holds an order-statistic rank and `se` holds a bootstrap SD | mark the bootstrap `p` and `se` columns, or publish one engine→colname stripe table in the help | high |
| CS-75 | absent results | an unavailable joint or equality test prints nothing at all | an unavailable adversarial p is posted as missing and the display gate lets it through, so a labelled test line prints with `.` where the p should be — under a caution that says the test is absent | add the `!= "."` guard at ado:9806/9948, or never post the analytic adversarial p under an engine | high |
| CS-20 | replication counts | `permreps()` is refused below 20 and cautioned below 39 | `reps()` has no floor and no resolution caution, so `reps(20)` prints percentile p's of resolution .10 in silence; and `reps(999)` under `permute()` errors "require(s) the bootstrap() option" | accept each spelling as the other's synonym and write one "how many replications" paragraph covering both | high |
| CS-114 | documentation | the `bootstrap()` entry says the joint tests are bootstrap quantities in both forms, with `e(boot_jt_family)` recording which | Example 13 still says "without it the joint tests are absent under `bootstrap()`", so a user pays the full studentized cost for a test the default path now delivers | rewrite Example 13's parenthetical; the option entry is canonical | high |
| CS-25 / CS-40 | safety() | in its advertised sense `safety(off)` downgrades a refusal to a caution and the run is unchanged | it also bypasses whole checks silently, silences a caution while the behaviour stands (`espostlength(0)` is still raised to 1), changes which numbers are computed (the collapse fallback), enables withheld output, is not honoured at all by the collapse refusals, and is set for the user by `contreatwithin` | classify every safety-gated site in the help as refusal-lifted / check-bypassed / caution-silenced / different-number / output-enabled / not-honoured, and print the disabled list at run time | high |

## 2. The catalogue

Grouped by the lens that found each item, which is also the theme grouping the author asked for: engines,
options, policy, stored results, display, interplay, documentation. Within each group items run high →
medium → low priority. Each item carries the two behaviours side by side with the verifier's corrections
already folded in, why a user stumbles, the pathways with their costs, whether the sweep judged the
asymmetry defensible, and a verification line recording what the verifier changed.

### 2.1 Engines

#### CS-1 — three p-value formulas ship in one package, and two of them appear in one printed table

**X (analytic, and the bootstrap's coefficient rows).** Analytic runs print a tail probability of the t/F/chi2
reference distribution. `bootstrap()` coefficient rows — percentile and BC alike — print
`2*min(#draws<=0, #draws>=0)/completed` with plain order-statistic counts and no add-one, so a p can print
exactly 0: ado:12791-12792, `m0 = sum(v :<= 0); m1 = sum(v :>= 0); cipp[1,kcol] = min((1, (2*min((m0,m1)))/nn))`.
The rationale is in the code at 12769-12784 and in the help at sthlp:504 — no add-one "keeps each one dual
to its own interval", measured at 6 boundary disagreements in 35,646 constructed cases against 36 for the
add-one alternative.

**Y (permute, and the bootstrap's own joint family).** `permute()` is not uniformly add-one, which is itself
a "the default flips when a bonus option is added" case: ado:13040-13050 reads
`if (exact) { plo = m0/n; phi = m1/n } else { plo = (1+m0)/(n+1); phi = (1+m1)/(n+1) }`, so under
`permexact` the formula is *literally the bootstrap's* plain `2*min/n` and permute's floor there comes from
self-membership (the observed assignment is in its own reference set), not from add-one; sthlp:407-410 says
so, and `permabsp` switches the same way at 13058. Meanwhile the bootstrap's **own** joint, equality and
adversarial tests — printed a few lines below its coefficient table — use add-one (6519, 6578, 12218,
12343, the last with the explicit comment "unlike the percentile p on raw coefficient draws"). So the count
is three formulas, not four: an analytic tail probability; plain `2*min/n` (bootstrap percentile/BC, and
`permexact`); and the add-one form (sampled permute, and the bootstrap's joint family, which borrows
permute's rule rather than adding a fourth).

**Why a user stumbles.** A user who learned that a wooldid rank p can be 0 (bootstrap) is surprised when
`permute()` floors at 2/(reps+1); a user who learned permute's floor is surprised that the bootstrap table's
0.0000 is literal. Inside one bootstrap run the coefficient p and the joint-test p below it are computed by
different rules, and until 2026-09-06 nothing in the output distinguished them: a grep of printed lines for
"add-one" returned exactly three, two of them permute-only and the third gated on
`e(boot_jt_family) == "fixedmetric"`, so under `boottstat` — where add-one bootstrapped test statistics sit
under a no-add-one coefficient table — no printed line named a convention at all.

**Current tree.** `b6ef398` closes the ado half: the footer at ado:10113 (`2d59e70`) now says "The
rank-based **coefficient** p-values have resolution 2/n ... The joint, equality and adversarial p-values are
add-one shares over their own listwise-complete replication counts (their smallest value is 1/(n + 1), never
0)", and it prints under every rank display and under `boottstat`. What remains is the intra-run split
itself, and the help, which still describes the 2/n resolution without the joint-test qualification
(sthlp:504) — see CS-121.

**Pathways.** (a) Unify on add-one everywhere: costs the bootstrap coefficient p its exact duality with the
printed percentile interval — the measurement at 12778 says add-one triples the boundary disagreements — and
moves every pinned bootstrap p in the battery. (b) Unify on no-add-one everywhere: the joint tests lose the
conservative-never-anticonservative guarantee the F033 record established for listwise-incomplete cases, and
sampled permute gains a reachable 0 with no floor statement (`permexact` cannot reach 0 under any
convention, so this half of the cost is sampled-mode only). (c) Keep all three and make the convention
visible per block: carry the ado's new footer qualification into the help, and add a "p convention" line to
the stored-results section. Display and text only; breaks nothing.

**Defensible?** Yes — each convention is defended on its own terms and the duality argument for the
coefficient p is measured. The inconsistency is that two of them print in one table.

**Verification.** Partly confirmed. Corrections applied: behaviour Y as filed said permute is uniformly
add-one — it is not, `permexact` uses the bootstrap's own formula; and "four rules" is three formulas. The
verifier also sharpened the finding in its own favour on the printed-disclosure count, which `b6ef398` has
since answered.

#### CS-2 — `set level` is honoured by the analytic tables and hard-coded to 95 by every resampling rewiring

**X (analytic).** The worker takes the level in force: ado:525-529, `if `level' != 0 & "`wdrepmode'" == ""
{ set level `level' }` then `local wdlvl = c(level)`; `e(level)` is posted from it (4554) and the interval
columns are named for it (3849, `lb_90ci` under `set level 90`). `wooldid_ptp` is passed
`level(`=cond(`level' == 0, c(level), `level')')` at 418, so the Roth slope interval follows `c(level)` too.
The help promises exactly this: sthlp:489, "the default is `c(level)`, normally 95", and "The analytic
intervals honour level() since 2026-09-06".

**Y (bootstrap and permute).** All six resampling subprograms resolve an unspecified `level()` to a literal
95, never `c(level)`: ado:6369, 6454, 7417, 7681, 8033, 8639 — the first with the comment "the shell's
level() defaults to 0 (meaning c(level)); resolve it here exactly as the coefficient-row rewiring does".
`e(perm_level)` (7642) and `e(boot_level)` (9155) are posted from that 95. The mislabelling mechanism is
that the rewired matrix inherits the **analytic** colnames — `matrix `wm' = e(wdidmainresults)` at 8658 —
while the interval values are recomputed at level 95 (8709), and the rename that would fix the header is
both guarded by `if `level' != 95` and written as a `subinstr ... "lb_95ci"`, so it neither fires nor could
match `lb_90ci` if it did (8772-8779; same shape at 7502-7510, 7997-7999, 8942-8943).

**Why a user stumbles.** `set level 90` is ordinary Stata practice and the help says wooldid follows it. The
user sees columns named `lb_90ci` and reads them as 90% intervals; under `bootstrap()` or `permute()` they
are 95% intervals. It is worse than silent: the footer prints the level literally from `e(boot_level)`
(ado:10100/10103), so the run prints "the ... 95% CIs" directly beneath a table headed `lb_90ci` — a visible
contradiction on one screen. Explicit `level()` is unaffected, because the shell `set level`s it and the two
paths then agree; the defect is confined to an unspecified `level()` under a non-95 session level, which is
why it has not surfaced in the battery.

**Pathways.** (a) Resolve `level == 0` to `c(level)` in the six subprograms, one line per site: breaks
nothing (explicit `level()` behaves as today), makes the rename branch fire correctly, and needs a test
pinning `set level 90` + `bootstrap()`/`permute()` — `test_level.do` (a969126) is the natural home. (b) Have
the shell resolve `level` once before any call. Cheaper-looking but a bigger blast radius than it appears:
the shell only executes `set level` when `level != 0` **and** `wdrepmode == ""` (525), and the worker's
restore logic at 68-70 keys on the `WOOLDID_LEVEL0` global, so making `level != 0` always true changes which
branch of that set/restore pair runs in repmode. (c) Document the exception and refuse a non-95 session level
with the engines unless `level()` is given — loud, preserves current numbers, and contradicts sthlp:489,
which shipped 2026-09-06 and already promises the opposite.

**Defensible?** No. It is a help/code contradiction introduced by F054 (a969126), which taught the analytic
path to honour `c(level)` and left the rewiring subprograms on their pre-existing default.

**Verification.** Confirmed (code-verified; the Stata box was at the process cap, so no run). The verifier
added the colname-inheritance mechanism at 8658 that the finding only implied, the footer contradiction, and
the correction that pathway (b) is more than "slightly larger". See also CS-62 (the `e()` half, and a third
`level` default inside `wooldid_ptp`), CS-36, CS-116, CS-117.

#### CS-3 — `fastnose` is refused outright under `bootstrap()` and is the silent default under `permute()`

**X.** `bootstrap()` refuses the pair with an error: ado:6183-6186, "bootstrap() may not be combined with
fastnose; fastnose requests point estimates only, leaving nothing for the bootstrap to attach inference to
(the replications use it internally)", and sthlp:485 lists it among the incompatibilities. Every bootstrap
replication and the observed run go through the full inference path. (Internally the bootstrap *does* use
fast mode for its own reps — ado:297 appends `fastnose` to the rep option list unless `boottstat`.)

**Y.** `permute()` without `permtstat` **or** `permanalyticse` sets it for the user: ado:5404,
`if "`permute'" != "" & "`permtstat'" == "" & "`permanalyticse'" == "" local fastnose fastnose`. The
consequence is real and was reproduced: on `dta_binary` with `att cluster(i) espre(2) espost(2)
permute(shuffle) permreps(20)`, `e(properties)` is `b`, `confirm matrix e(V)` returns rc 111, `confirm
matrix e(b_es_att)` returns rc 111, and `test _b[att]=0` returns rc 321 — so `e(V)`, the event-study
matrices, the honestdid/pretrends handoff, `lincom`, `test` and `esttab` standard errors are all absent from
a default permute run, and `permanalyticse` is the opt-in that restores them.

**Why a user stumbles.** The same three words — fast mode — name a forbidden combination under one engine
and an unrequested default under its sibling. A user moving a working `bootstrap()` script to `permute()`
loses `e(V)` and the honestdid handoff without asking for anything. The state itself is not hidden: the
permute shell blanks `e(wdid_fastnose)` (7657) and rewrites `e(wdid_vcespec)` (7655) so the header reads
"... analytic standard errors not computed (permanalyticse restores them)", and the footer repeats it at
10030. What is undisclosed at run time is the *consequence* — which `e()` matrices went away — stated only
in the help (sthlp:417).

**Pathways.** (a) Compute analytic SEs by default under `permute()` and add a `permfast` opt-out: costs
every permute run the margins variance work the inference does not consume, which is the reason the default
exists, and is materially slower on poisson event studies. (b) Let `bootstrap()` accept `fastnose` the way
permute does, since bootstrap inference also does not consume the observed analytic SE unless `boottstat`:
costs `e(V_analytic)`/`e(df_r_analytic)` and the observed-run comparator's SEs, and `boottstat` would still
have to refuse it. (c) Keep both defaults and unify the vocabulary and disclosure: name the dropped `e()`
matrices in the header note that already exists, and say in the `fastnose` help entry that `permute()` turns
it on by default rather than burying it mid-paragraph.

**Defensible?** Yes on the behaviour — permute genuinely does not consume analytic SEs and the bootstrap's
observed run does. Not on the vocabulary, or on the silence about the lost `e()` surface.

**Verification.** Partly confirmed, with a live run for the consequence half. Corrections applied: the gate
is ado:5404 and exempts `permanalyticse` as well as `permtstat`; and the state is named in the header, not
only the footer, so pathway (c)'s original "put it in the header" half is already shipped — what remains is
the help text and the `e()` inventory. The verifier judged "high" overstated on that basis. See CS-23 (the
ordering trap), CS-94 (the help's incompatibility list).

#### CS-4 — `e(V)` carries four different provenances and one help paragraph describes all of them

**X.** On an analytic run `e(V)` is exactly `diag(se^2)` — ado:4264, `matrix `ev' = diag(hadamard(...))` —
and the help's Incompatible Post-Estimation Commands paragraph (sthlp:1186-1190) says `e(b)`/`e(V)` are "for
display purposes only" with "all off-diagonal elements set to 0", so `lincom` and `test` are invalid. That
paragraph's only carve-out is `e(b_es_type)`/`e(V_es_type)`; it never mentions `bootstrap()` or `permute()`.

**Y.** There are four provenances, and the blanket warning is wrong for one and materially incomplete for
two more. (1) Analytic: `diag(se^2)`, warning true. (2) `permute()` with `permtstat` or `permanalyticse`:
the *same* `diag(se^2)`, so the warning is literally true — but the real hazard is that the printed p-values
and intervals are permutation quantities while `e(V)` is analytic, which the footer says (10028-10029) and
the paragraph does not. (3) Default `permute()` fast mode: no `e(V)` at all, warning moot. (4) `bootstrap()`:
the full bootstrap covariance is reposted (9047-9067) with the analytic one preserved as hidden
`e(V_analytic)`, and `e(df_r)` is deliberately aligned to the reference distribution actually reported
(9173-9189) *so that* a hand `test`/`lincom` reproduces the `bootdisplay(wald)` family — the footer promises
this at 10111-10113. There the warning is simply wrong, in the user's disfavour.

**The `df_r_analytic` gap.** sthlp:503 and sthlp:1061 both say the analytic residual df "is preserved as
hidden `e(df_r_analytic)`". A grep returns exactly one posting site, ado:9187, and it sits inside the `else`
branch — `bootdf(normal)`. In the `if `wdbdf' < .` branch, which is `bootdf(clusters)`, the default whenever
`cluster()` is given, `ereturn scalar df_r = `wdbdf'` overwrites the analytic value and nothing preserves it.
The help's promise fails on the more common branch.

**Why a user stumbles.** The user reads one sentence — never run `test` after wooldid — and applies it
everywhere. Under `bootstrap()` that advice is wrong in their disfavour (the test is valid and reproduces a
stored family); under `permute()` it is wrong in the other direction (the test runs and returns an analytic
answer while every printed number is permutation-based). One object name carries four provenances with no
marker in `e()` beyond each engine's own scalars.

**Pathways.** (a) Unify on the bootstrap behaviour by posting a genuine covariance on the analytic path too
— large: it needs an analytic cross-estimand covariance that does not exist today (`e(V_es_*)` covers the
event-study block only). (b) Unify on the analytic behaviour, keeping the bootstrap covariance in
`e(V_boot)` only — breaks the deliberate `df_r` alignment and the footer's promise. (c) Keep all four and
post one macro, `e(vsource)` in {`analytic_diag`, `analytic_diag_under_permutation`, `none`,
`bootstrap_full`}; rewrite the Incompatible Post-Estimation paragraph per engine; and move the
`ereturn hidden scalar df_r_analytic` above the `bootdf()` if/else so the help's claim holds. Cheap, breaks
nothing.

**Defensible?** Yes on the behaviour — the bootstrap repost and the `df_r` alignment are deliberate and
well-reasoned (the comment at 9178-9186 cites a measured 12-cluster case). The inconsistency is the single
un-engine-aware help paragraph plus the absent provenance marker.

**Verification.** Confirmed. Correction applied: the title as filed said "true for only one of them" — the
accurate count is four provenances, the warning true for two, incomplete for those two, and wrong for
`bootstrap()`. The `df_r_analytic` branch gap is the verifier's sharpest addition and is a one-line fix. See
CS-59 (the same object from the stored-results lens, adding the missing permute flag) and CS-77 (the footer
registers).

#### CS-5 — the adversarial line test runs under `permute()` but is absent under `bootstrap()` without `boottstat`, though the substrate now exists

**X.** When `adversarial` is requested, `permute()` runs `wooldid_permadv` — `wooldid_advW` on the
permutation-covariance metric of the family's horizon columns — at ado:7874-7937, posts
`adversarial_pval_<fam>[_s<sg>]` plus `e(perm_advW)`/`e(perm_advWobs)`, and the footer states it at 10039
beside the joint and equality tests. Two qualifications: it is computed only under `esfixedbaseperiod`
(`esfflag`, 7891) — without it any worker-posted analytic p is blanked (7918) and only the V-free display
line survives — and it tests the sharp null of no effect at any horizon, not the composite on-a-line null
(the footer says so). It is mode-independent, needing no per-replication covariance.

**Y.** `bootstrap()` computes the adversarial p only under `boottstat`. Line 396 passes `awacc()`/`aobs()`
only when `boottstat` is non-empty; `wooldid_bootjnt` sets
`wdb_advok = ("`bts'" != "" & "`awacc'" != "" & "`aobs'" != "")` at 6648, so without it every requested
adversarial p is set missing, `e(boot_adv_absent)` is posted, and a caution reads "bootstrapping them means
bootstrapping the test statistic, which needs boottstat" (6691). But since `bd8ee3c` (F034, 2026-09-06) the
bootstrap **does** run its joint and equality tests without `boottstat`, on the fixed-metric rank
construction the permutation engine uses (6512, `e(boot_jt_family) = fixedmetric` at 6570; the equality
chains post `e(boot_sgW)` at 6813). The adversarial test was not carried along, and the bootstrap footer at
10114 names "the joint pre-trend tests (and the across-subgroup equality tests...)" and no adversarial.

**Why a user stumbles.** `adversarial esfixedbaseperiod jointtests` under `permute()` prints three test
families; the same command under `bootstrap()` prints two, with the third silently missing and explained
only by a caution that does not survive replay (see CS-42). The caution's stated reason — "needs boottstat"
— is now true of the studentized form only, not of the machinery in the tree.

**Pathways.** (a) Extend the fixed-metric rank test to the adversarial statistic under `bootstrap()` without
`boottstat`, mirroring the joint-test block: the substrate is `wooldid_permadv` at 13377-13425 and the
bootstrap's fixed-metric quadratic-form twin at 13325ff, so this is one more accumulator column and a new
cert pin. The catch the verifier raises is that it is not purely additive — the permutation construction
tests a *different null* (sharp no-effect-at-any-horizon) than the analytic on-a-line test the same option
name delivers elsewhere, so importing it would need its own footer sentence. (b) Withdraw the fixed-metric
joint tests too, so bootstrap without `boottstat` reports no quadratic-form test at all — internally
consistent, a regression against F034. (c) Keep as is and reword the caution to say the joint and equality
tests have a non-studentized form and the adversarial one does not, so the user is not left inferring a
general rule from a specific gap.

**Defensible?** Not as it stands — the caution's reason no longer matches the tree. But (a) is less
automatic than it looks, because of the null change.

**Verification.** Partly confirmed. Corrections applied: the evidence cite ado:12218 is the `boottstat`
studentized construction, not the fixed-metric substrate (that is 13377-13425); permute's adversarial p is
conditional on `esfixedbaseperiod`; and it tests a different null from the analytic test it replaces.

#### CS-6 — `boottstat` and `permtstat` studentize different sets of objects, and permute's contrast column is labelled `t_perm` while being unconditionally blank

**X.** `boottstat` studentizes everything: the coefficient rows, the `gteffects` and `e(cfx*)` cells, the
pairwise subgroup contrasts (each by its own replication pairwise standard error, ado:6939-7018, displayed
as `t_pt` with `se`→`se_analytic`), and the joint, equality and adversarial tests, which become bootstrapped
test statistics (sthlp:499; caution at 9169).

**Y.** `permtstat` studentizes only the coefficient rows and the cfx cells. The pairwise subgroup contrasts
stay coefficient differences in both modes, and their statistic column is set missing *unconditionally* —
`matrix `pw'[`rr', 3] = .` at ado:7977 sits outside every conditional — yet the column is then renamed
`t_perm` at 7993-7994, with the same rename text the main table uses at 7505-7506 where the column *is*
filled under `permtstat`.

**Why a user stumbles.** The parallel names promise parallel scope. A user who adds `permtstat` expecting
the whole surface to studentize finds the contrast table's `t_perm` column empty in every row. Two texts
that describe `t_perm` generically contradict what that table shows: sthlp:414 ("the statistic column
t_perm is the permuted t-ratio under permtstat and blank otherwise") and the footer at ado:10025 (`t_perm` =
estimate/se_analytic, no exception).

**Pathways.** (a) Extend `permtstat` to the contrasts, studentizing each by the replication's own pairwise
analytic SE as `boottstat` does — the permute engine would have to carry a pairwise-SE accumulator it does
not have, and the same degeneracy caveat that limits `boottstat`'s contrasts under cluster ≥ i would apply
(CS-41). (b) Restrict `boottstat`'s reach to match — a large regression, not worth it. (c) Stop labelling a
permanently blank column `t_perm`: drop it from the contrast table, or leave it named `t`, and carve the
contrast table out of sthlp:414 and the ado:10025 footer sentence.

**Defensible?** Yes on the scope difference; not on a column header that can never be filled.

**Verification.** Partly confirmed. Correction applied: the finding said the exception is undocumented and
uncautioned — sthlp:454 and sthlp:413 both state it, and the footer at ado:10062 half-signals it (though
only when more than one subgroup exists, and it never names the blank column). The residual finding is the
labelling plus the two generic statements. Filed medium; verifier: low-to-medium.

#### CS-7 — a bootstrap table's column literally named `p` is a percentile rank

Anchor for the p/se half of the column-naming cluster; see **CS-26** (options lens) and **CS-74** (display
lens), which were verified independently with different citations, and **CS-67** (the cell tables, which no
engine relabels).

**X.** `permute()` rewrites three of the six columns — `se`→`se_analytic`, `t`/`z`→`t_perm`, `p`→`p_perm` —
at ado:7504-7507 (main), 7543-7546 (subgroup) and 7992-7995 (pairwise), so every column announces that it is
not the analytic quantity of the same position.

**Y.** `bootstrap()` renames the statistic (`t_wald`/`z_wald` at 8853-8858, `t_pt` at 8911-8913) and, under
`boottstat` only, `se`→`se_analytic` (8913). It never renames `p`: the header reads `p` whether the number
is a Wald tail probability, a percentile rank, a BC rank, a BCa rank or a studentized rank. `se` likewise
keeps its plain name while holding the bootstrap standard deviation (the 2026-09-04 ruling recorded at
8797-8804 and in sthlp:504).

**Why a user stumbles.** Users copy tables into papers. A permute table's headers say the numbers are
permutation quantities; a bootstrap table reads `estimate se t_wald p lb_95ci ub_95ci` and a reader takes
`p` for the tail probability of `t_wald` — which sthlp:504 explicitly says it is not. The one column whose
provenance most needs flagging is the one left unflagged.

**Pathways.** (a) Rename per family (`p_pctl`/`p_bc`/`p_bca`/`p_wald`/`p_pt`, `se_boot`): breaks every
`colnumb(..., "p")` read and every pinned cert; `wooldid_pfloor` falls back to positional column 4 when it
cannot find "p" (9641-9642), so it degrades rather than errors. (b) Stop renaming under `permute()` too and
put all provenance in the footer — cheaper, loses information already shipped, reverses F053. (c) Keep both
and state the rule in the help's output section: name the column whenever the *quantity* changes (permute's
rule) and additionally whenever the *convention* changes (which catches the bootstrap's `p`).

**Defensible?** No — the half-applied rule is worse than no rule, because `t_wald` signals that renaming
happens here.

**Verification.** Confirmed. The verifier added that provenance *is* recoverable from `e(boot_display)` and
`e(boot_stat)`, just never from the stripe, and that `boottstat`'s `se_analytic` rename means the bootstrap's
rule is "rename when the denominator changes" against permute's "rename whenever the engine changed the
number".

#### CS-8 — the "p = 0 means p < 2/n" disclosure is bootstrap-only, and `permrawp` reintroduced a printable 0 the same day

**X.** Under `bootstrap()`, any table containing an exact 0 gets a line beneath it plus a footer resolution
sentence — `wooldid_pfloor` at ado:9632-9646, called for the main, subgroup and pairwise tables (9767, 9916,
9961).

**Y.** `wooldid_pfloor` exits immediately unless `e(boot_completed)` is set (9635), and the comment above it
justifies the gate: "Bootstrap only: the permutation add-one p never prints 0." That was true of the default
and of `permabsp` — but not of `permrawp`, shipped in the same commit (`dc3f345`, 2026-09-06), whose p is
`ma/n` with no add-one (13056-13059) and is exactly 0 in sampled mode when no replication reaches the
observed absolute statistic. There is no under-table line and no footer resolution sentence for that case;
the footer's raw-convention paragraph (10013-10016) describes `permrawp` as "the number ritest and Stata's
permute prefix report" and stops. Under `permexact` the observed assignment is in its own reference set, so
the hole is sampled mode only.

**A second, sharper wrinkle.** The F036 low-rep caution keys on `absp` alone —
`local wdp_floor = cond("`absp'" != "", 1/(`permreps' + 1), 2/(`permreps' + 1))` at ado:5816 — and `permrawp`
arrives in a separate `rawp()` argument. So a `permrawp` run with `permreps()` between 20 and 39 is told
"the smallest attainable permutation p-value is [2/(reps+1)] ... so no result can reject at the 5% level",
a floor `permrawp` does not have and which the same run may contradict by printing 0. `permside(left|right)`
is missed the same way (one-sided floor 1/(n+1), caution reports 2/(n+1)).

**Why a user stumbles.** A user who learned from bootstrap tables that wooldid explains its zeros adds
`permrawp` to cross-check against `ritest` and gets 0.0000 with no explanation — and, because the footer
emphasises that the default convention has a floor of 2/(reps+1), may read the 0 as impossible rather than
as the resolution limit of a convention they just switched on.

**Pathways.** (a) Generalise `wooldid_pfloor` to key on "a rank-based p is in force" rather than on
`e(boot_completed)`, computing 2/completed for the bootstrap and 1/completed for `permrawp`; and thread
`rawp`/`side` into the cond() at 5816 so the low-rep caution names the right floor. A few lines, breaks
nothing, and it future-proofs the gate. (b) Refuse `permrawp` in sampled mode or floor it at 1/n — defeats
the option, which exists for ritest comparability. (c) Add one sentence to the footer's raw-convention
branch saying a printed 0 means no replication reached the observed statistic, resolution 1/completed.

**Defensible?** No — it is a gate whose stated rationale went stale the day it was written.

**Verification.** Partly confirmed. Correction applied: "no disclosure at all" is too strong — sthlp:481
states the floors including "0 under permrawp", so the gap is on-screen disclosure, not documentation. The
low-rep caution's `absp`-only arithmetic is the verifier's addition and turns an inferred stumble into an
on-screen contradiction. See CS-80, CS-122.

#### CS-9 — `permute()` enforces a replication floor and a completion floor; `bootstrap()` enforces neither, and their failure cautions fire at opposite thresholds

Same asymmetry as **CS-20** (options lens, confirmed there) seen from the engines side, plus the failure
half. `permreps()` must be at least 20 (ado:5807-5810), the run is refused if fewer than max(20, R/2)
replications succeed (7364-7367), a caution fires below 39 (5813-5817), and the failure caution fires only
above 10% failures (7372). `reps()` has no floor by author ruling (6053-6055), no completion floor, and any
single failed replication draws a caution (8614); its resolution disclosure is footer-only (10108).
**Adds:** the failure cautions invert — one failed bootstrap replication in 1000 is loud, up to 10% failed
permutation draws are silent. Pathways: (a) give the bootstrap a floor and a resolution caution — contradicts
the standing F030 ruling; (b) drop permute's floor to a caution — the engine refuses inference below 20
completions downstream anyway, so this only moves the failure later; (c) one "how many replications"
paragraph covering both engines' floors, thresholds and resolutions. Defensible in part (the two floors have
different statistical status). *Not verified (cap).* See CS-32, CS-45, CS-89, CS-123.

#### CS-10 — the few-cluster cautions live on the bootstrap engine only, including the one that recommends `permute()`

Same fact as **CS-39** (policy lens, partly confirmed there); its corrections apply — permute is not wholly
silent, but its two cautions (`permtstat` at `cluster() == i`, ado:5902; the sampled-draw p-floor, 5811-5817)
are gated too narrowly to reach the user who arrives by following the bootstrap's own referral.
**Adds:** the exchangeability limitation is stated in the help's permute entry with a measured 6-of-30
rejection rate at nominal 5% and never fires as a run-time trigger on the default statistic. *Not verified
(cap)* for this framing; the anchor was verified.

#### CS-11 — "every family is computed and stored whichever is displayed" is a bootstrap promise permute does not make, and BCa is the bootstrap's own exception

Overlaps **CS-76** (footer, confirmed) and **CS-115** (help self-contradiction, confirmed); both apply here.
**Adds the permute half:** permute's p conventions are mutually exclusive flags (ado:5741-5743) and only the
selected one lands in `e(perm_p)` — there are no per-convention twin tables and no `e(perm_p_absolute)`. The
one-sided shares `e(perm_p_lo)`/`e(perm_p_hi)` are always posted, so the default equal-tailed p is
reconstructable under any convention, but the absolute-value count is not, so switching to `permabsp` after
the fact means recomputing from `e(permdist_b)`. Pathways: (a) post the absolute-value count as one extra
kernel row — the kernel already computes `m0`, `m1` and `ma` per element, so every convention becomes
reconstructable at negligible cost; (b) drop the bootstrap's all-families storage to match — loses the
documented hidden twins; (c) state the rule per engine in the stored-results help, with the BCa exception
beside it so it is not read as a contradiction of the neighbouring promise. Defensible. *Not verified (cap).*

#### CS-12 — the one-sided p, its Monte-Carlo SE and its Clopper-Pearson interval are a first-class option under permute and a hidden matrix under bootstrap

Same fact as **CS-56** (stored-results lens, confirmed there); its corrections apply, including that the
hidden `e(boot_ptails)` tracks the percentile family regardless of what `bootdisplay()` actually displayed,
so under `bootdisplay(bc)`/`(bca)` there is no MC-precision record for the displayed p at all.
**Adds:** the asymmetry was a deliberate same-day scoping decision (F064/F052, `dc3f345`), so pathway (a)
— promote `e(boot_ptails)`, add a footer sentence, optionally add `bootside(left|right)` — is a resumption
rather than a reversal. *Not verified (cap)* for this framing; the anchor was verified.

#### CS-13 — the `gteffects` and continuous-treatment cells are rewired, blanked or left analytic by engine, and the F053 relabelling never reached them

**X.** Analytic runs fill the `gteffects` columns `estimate/se/p/lb/ub` (ado:10848) and the `e(cfx*)`
columns `b/se/t/pvalue/ll/ul` from margins. `bootstrap()` rewires both to the displayed family
(`wooldid_bootcellsx`, 6346-6435), recomputing the cfx `t` column as b/se under each family (6420).

**Y.** `permute()` blanks the `gteffects` inference columns entirely — the cells are keyed by cohort onset,
which the replications reassign, so they have no permutation distribution (`e(perm_gt_blank)`, footer 10064)
— and permutes the cfx cells only for assignment-invariant lattices (8040-8085). Under `bootdisplay(bca)`
the bootstrap cells themselves fall back to the bias-corrected family with a caution (6357-6358). In every
one of these paths the cell matrices keep their original column names — `matrix colnames `cf' = `:
colnames e(`nm')'` at 8058, twin at 6417 — so the cfx statistic column is still `t` and the p column still
`pvalue`, though the main tables were relabelled by F053 and by the permute rewiring.

**Why a user stumbles.** `gteffects` is the same option under all three engines and returns a table of the
same shape; whether its inference columns are analytic, bootstrap or empty depends on an engine choice made
elsewhere on the command line, and only permute prints a footer line about it. The unrenamed cell columns
then let a user carry a `t` column out of a percentile-bootstrap run believing it is the test statistic
behind the `pvalue` beside it — the confusion F053 shipped to prevent in the main tables.

**Pathways.** (a) Extend the F053 relabelling to the cell matrices (`t`→`t_wald`/`z_wald`/`t_pt`/`t_perm`,
`pvalue`→`p_perm` under permute) — two lines per site, but it breaks user code reading those columns by name
and every cert that pins the stripe. (b) Leave the stripes and add a footer sentence per engine saying what
the cell columns hold, matching the permute footer's existing `gt_blank`/`cfx_blank` lines. (c) Blank the cfx
statistic column under the rank displays the way permute already blanks it without `permtstat` (8056), so no
unlabelled Wald object escapes — costs a column some users read. Defensible in part. *Not verified (cap).*
See CS-67, CS-79, CS-96, CS-102.

#### CS-14 — `e(V_es_*)`, the honestdid handoff, silently changes engine and only the bootstrap discloses it

**X.** On an analytic run `e(b_es_type)`/`e(V_es_type)` are the analytic event-study coefficients and
covariance, documented as the valid inputs for `lincom`, `test` and a positional honestdid handoff.

**Y.** `bootstrap()` replaces `e(V_es_*)` with the bootstrap event-study covariance (ado:9070 onwards; the
analytic one kept as hidden `e(V_es_analytic_*)`), so a subsequent honestdid call performs Rambachan-Roth
robustness on a bootstrap covariance without the user asking — the bootstrap help paragraph mentions this in
the `pretrendpower()` context, the `e(b_es)` stored-results entry does not. Under default `permute()` the
matrices are not posted at all, so the handoff disappears; under `permanalyticse` or `permtstat` they come
back **analytic**, beneath a table whose every printed number is permutation-based (footer 10028-10029).

**Why a user stumbles.** honestdid is a positional handoff: the user passes `e(b_es_att)`, `e(V_es_att)` and
`e(numpre_att)` and gets a robustness region. Nothing in that call reveals which covariance was inside. The
same script yields analytic robustness, bootstrap robustness, or an error about a missing matrix, depending
on an inference option added at the far end of the wooldid command line.

**Pathways.** (a) Post a provenance macro, `e(es_vsource_att)` in {analytic, bootstrap,
analytic_under_permute}, and name it in the stored-results entry — trivial, breaks nothing. (b) Caution the
handoff geometry check when the covariance is not analytic — noisy, and the bootstrap covariance is arguably
the better input under `bootstrap()`. (c) Make `permute()` post a permutation-covariance `e(V_es_*)` so the
surface exists under all three engines — it would reuse the metric the joint tests already build, but calling
it a sampling covariance for honestdid's purposes is a substantive claim needing its own research.
Defensible. *Not verified (cap).*

#### CS-15 — the analytic joint-test refusal targets exactly the clustering `permute()` requires, and deliberately does not gate the adversarial test

Same phenomenon as **CS-41** (policy lens, partly confirmed there) and **CS-99**; those corrections apply.
**Adds two facts.** First, `permute()` itself requires `cluster()` to nest `i` (ado:5890-5894), so every
`permute() jointtests` run sits inside the configuration the analytic screen refuses — by construction, not
by accident. Second, within the analytic engine the same screen deliberately skips the adversarial test
("adversarial is deliberately not gated here", ado:1492), so the adversarial line test runs on the same
rank-deficient clustered covariance without complaint. Pathways: (a) extend the screen to `adversarial`
analytically — internally consistent, but the adversarial statistic's sensitivity has not been measured the
way the joint tests' was, so it would be a refusal without evidence; (b) caution on the carve-out branch when
a run under an engine sits in a configuration the analytic screen would have refused; (c) document that
permute's own cluster requirement puts every permute run inside the refused analytic configuration, so the
two are not alternatives on the same footing. Defensible. *Not verified (cap).*

#### CS-16 — the joint-test lines print identical wording whichever engine produced the number, while the adversarial line branches — on permute only

**X.** The adversarial line does branch: `if "`e(perm_method)'" == ""` prints the analytic composite-null
sentence, else the permutation sharp-null sentence (ado:9806-9808, and 9948-9950 for subgroups).

**Y.** The joint pre-trend lines (9789-9800), the per-subgroup joint lines (9928-9945) and the
across-subgroup equality lines (9967-9977) print one fixed analytic wording under all three engines — "p-value
for joint test that all pre-treatment att horizon estimates are 0" — though under bootstrap they are
fixed-metric rank tests of recentred draws (6520-6577) or bootstrapped statistics under `boottstat`, and
under permute quadratic-form permutation tests. And because the adversarial branch tests only
`e(perm_method)`, a **bootstrap** run falls into the analytic branch, so a `boottstat`-bootstrapped
adversarial p prints with the analytic composite-on-a-line wording.

**Why a user stumbles.** These p-values are what most users screenshot; the line above the number is the only
label it carries into a slide deck, and that label is engine-agnostic. The disclosure lives in footers
further down. For the adversarial line under `boottstat` the wording is actively wrong.

**Pathways.** (a) Append the engine to each test line — "(bootstrap, fixed-metric rank, add-one)" /
"(permutation rank)" / "(Wald)" — display-only, but it touches every joint/equality/adversarial cert pin. (b)
One engine banner above the joint-test block instead of per line — cheaper, fewer pins move. (c) Minimum fix:
make the adversarial branch test for the bootstrap as well, so a `boottstat` adversarial p gets its own
sentence rather than the analytic one. *Not verified (cap).* See CS-82.

#### CS-17 — the two engines use different option-name schemes and store parallel quantities under non-parallel names, including two `e(*_stat)` macros that answer different questions

**X.** `bootstrap()`: `reps()`, `bootalpha()`, `bootdf()`, `bootdisplay(spec)`, `ties`, `bootsave()`,
`boottstat`. Stored: `e(bootdist)`, `e(boot_p_percentile)`, `e(boot_ci_percentile)`, `e(boot_navail)`,
`e(boot_display)`, `e(boot_jt_family)`, and `e(V_analytic)`/`e(V_es_analytic_*)`/`e(df_r_analytic)` as the
archive of what it replaced. **`e(boot_stat)` records the displayed column label** — `t_wald`, `z_wald`,
`t_pt` (9156).

**Y.** `permute()`: `permreps()`, `permexact`, `permstrata()`, `permprobs()`, `permsave()`, `permtstat`,
`permanalyticse`, and the p convention as three mutually exclusive bare flags rather than one keyword.
Stored: `e(permdist_b)` plus the legacy `e(permdist)`, `e(perm_p)`, `e(perm_ci)`, `e(perm_n)`,
`e(perm_pconv)`, no `e(perm_jt_family)`; the archive is a set of **tables**
(`e(perm_analytic_mainresults)` and friends) rather than covariances — and under the fast-mode default those
"analytic" tables contain no analytic standard errors. `e(perm_blank)` is posted while the bootstrap's
identical structural-zero caution (8615) posts no `e(boot_blank)`. **`e(perm_stat)` records the permuted
quantity** — `b` or `t` (192) — not a column label. `bootdf()` and `ties` have no permute analogue;
`permanalyticse` has no `bootanalyticse`. User weights are documented at length for the bootstrap (sthlp:485)
and not at all for permute, which forwards `[wtype=wexp]` to every replication (7075-7078).

**Why a user stumbles.** Anyone scripting both engines learns two naming systems for the same ideas, and the
two `e(*_stat)` macros are the trap: they look like a matched pair and answer different questions, so code
branching on `e(boot_stat)=="t"` and `e(perm_stat)=="t"` reads a column label in one case and a design choice
in the other.

**Pathways.** (a) Alias the divergent names — accept `reps()` under permute, add `bootanalyticse`, add
`permdisplay(...)` mirroring `bootdisplay()`, post `e(perm_jt_family)` and `e(boot_blank)`, and rename one of
the two `*_stat` macros keeping the old as a hidden alias: a deprecation surface and help churn, breaks
nothing if the old names stay accepted. (b) Leave the option names and unify only the stored surface, which
is what downstream code reads — cheaper, and the `*_stat` collision is the highest-value item in it. (c) One
"the two resampling engines side by side" table in the help mapping each option and each `e()` name to its
counterpart, noting explicitly that `e(boot_stat)` is a label and `e(perm_stat)` a choice. Defensible in
part. *Not verified (cap).* See CS-21, CS-33, CS-73, CS-124, CS-125.

#### CS-18 — whether the statistic column says `t` or `z` is decided by a different fact in each engine

**X.** Analytically the label is `t` unless poisson, where it is `z` (ado:3849-3850, twins at 3976-3977,
4019-4020) — the trigger is the estimator. **Y.** Under bootstrap it is `t_wald`/`z_wald` (or `t_pt`),
decided by `local wdblab = cond(`wdbdf' < ., "t", "z")` at 8651, where the df come from `bootdf()`, whose
default is `cond("`cluster'" != "", "clusters", "normal")` (392) — the trigger is whether you clustered.
Under permute it is always `t_perm` (7505-7506), filled or not, whatever the estimator or clustering. So the
same specification flips the letter in both directions: a poisson analytic run says `z` and its clustered
bootstrap says `t_wald`; a robust non-poisson run says `t` and its unclustered bootstrap says `z_wald`.
Pathways: (a) make the bootstrap label follow the estimator and carry the reference distribution in the
footer, where it already appears (10096) — loses the signal that `bootdf()` changed the critical values;
(b) make the analytic label follow the reference distribution as the bootstrap does — arguably more correct,
and it moves the poisson/non-poisson labelling in every analytic cert; (c) say in the help which fact drives
the letter under each engine, and put the reference distribution in the header so the letter is never the
only clue. Defensible — the two engines genuinely have different things to name. *Not verified (cap).*
Twinned with **CS-88** (display lens).

#### CS-19 — `pretrendpower()` is supported under `bootstrap()` and refused under `permute()`, though permute already builds a covariance metric

Same fact as **CS-95** (interplay lens, confirmed there), where the refusal text, the call sites and the
permutation-covariance metric that contradicts the refusal's stated reason are all verified; those citations
supersede this item's. **Adds:** the author's recorded decision — "the author allowed the refusal rather than
a permutation-covariance substitute" (ado:5767-5769) — and F039's "do nothing (documented)" ruling on the
Roth slope under bootstrap, which together mean pathway (c), rewording the refusal to say the diagnostic
needs a *sampling* covariance rather than that randomization inference produces none, is the only one not
already ruled on. Defensible. *Not verified (cap)* for this framing. See also CS-92.

### 2.2 Options

#### CS-20 — `reps()` and `permreps()`: one concept, two names, two requirement regimes, one floor

**X.** `bootstrap()`: `reps(#)` is required (ado:6056-6057, r198 "bootstrap() requires reps(#) with a
positive number of replications", preceded by the comment "reps() is mandatory (author ruling: no floor -
low-rep runs are legitimate)"), has no default, no minimum and no low-replication caution of any kind.
sthlp:487 says so: "No minimum is imposed, so small values can be used to time a specification". `reps(20)`
therefore runs to completion and prints percentile p-values of resolution 2/20 = .10 with nothing said —
and the footer's resolution line is conditional, not routine: `wooldid_pfloor` prints only when a displayed
p is exactly 0, and exits for `bootdisplay(wald)` unless `boottstat`, so a `reps(20)` run whose p's are all
non-zero is told nothing at all about resolution. There is no enumeration alternative.

**Y.** `permute()`: `permreps(#)` is required unless `permexact` (5851-5854), carries a hard floor of 20
enforced with r198 (5807-5810, "permreps() must be at least 20: the engine refuses inference on fewer than
20 completed replications, so a smaller request can only ever burn the run"), and fires a caution whenever
2/(permreps+1) > .05 — below 39 replications — naming the counts needed for .05/.02/.01/.002 (5811-5818,
comment-tagged "author ruling 2026-09-06: caution, floor stays 20"). `permexact` needs no count. The two
floors are the same order of magnitude — 2/20 = .100 for the bootstrap, 2/21 = .095 for permute — but only
one engine says so up front. Under `permute(probabilities)`, `permexact` is refused (5838-5841) and
`permreps()` is unconditionally required (5842-5845), sharpening the asymmetry further.

**Why a user stumbles.** A user who learned `reps(999)` on a bootstrap run types `reps(999)` under
`permute()` and is told "reps() require(s) the bootstrap() option" (6042). Going the other way,
`permreps(25)` is accepted but cautioned (2/26 = .077 > .05) while `reps(25)` is accepted silently with a
strictly worse p-resolution. The engine that protects the user is the one whose option name is harder to
guess.

**Pathways.** (a) Unify on permute: accept `reps()` as a synonym for `permreps()` and vice versa so neither
spelling errors, and port the resolution caution to the bootstrap. Cost: reverses the F030 "do nothing"
ruling, and the caution would fire on legitimate timing runs (mitigable by suppressing it under a small
explicit `reps`). (b) Unify on bootstrap: drop `permreps()`'s floor and caution. Cost: reverses F036 and
re-opens the "p printed at its own floor" problem the caution was built for. (c) Keep both and write one
side-by-side "how many replications" paragraph giving the floor, caution and resolution formula per engine;
add the cross-name synonyms only. Breaks nothing.

**Defensible?** Yes — a permutation p has a hard combinatorial floor and a bootstrap p only a Monte-Carlo
one, so the two floors have different status. The option-name asymmetry is not defended by that.

**Verification.** Confirmed, line by line, including both quoted error strings and the sthlp text. The
verifier added that the bootstrap footer's resolution line is conditional, which makes the silence more
complete than the finding claimed. Anchor for CS-9, CS-44, CS-89, CS-123.

#### CS-21 — the p-value convention is one enum on the bootstrap and four separate flags on permute, and each engine exposes what the other hides

**X.** `bootstrap()`: a single enum, `bootdisplay(percentile|bc|bca|wald)`, choosing p-value and interval
together, validated by `inlist` at ado:6148. No option selects a one-sided p, and none selects the add-one
convention. The one-sided shares are computed and stored, but only as an undocumented hidden matrix —
`ereturn hidden matrix boot_ptails` (9133, rows `p_lo p_hi se_p cp_lo cp_hi`), zero mentions in
`wooldid.sthlp`.

**Y.** `permute()`: the convention is chosen by four separate switches — the default equal-tailed,
`permabsp`, `permside(left|right)`, `permrawp` — with mutual exclusivity enforced by counting flags
(5741-5743: "permside(), permrawp and permabsp are alternative p-value conventions; specify at most one"),
plus `permtci`, which adds a further interval without replacing anything. The one-sided shares are
first-class: `e(perm_p_lo)`/`e(perm_p_hi)` are documented (sthlp:477, 1066) and the convention in force is
recorded in `e(perm_pconv)`. There is no bc/bca/wald family concept.

**Why a user stumbles.** A user who learned from the bootstrap that "the convention lives in one option"
looks for `permdisplay()` and finds nothing; a user who learned `permrawp` gives ritest's c/n convention
looks for `bootrawp` and finds nothing. The mutual-exclusion error names three options that the synopsis
gives no reason to read as one family (they are adjacent at sthlp:72, but nothing says they are exclusive
before the run errors). And the bootstrap's one-sided p exists in `e()` but cannot be found from the help.

**Pathways.** (a) Unify on the enum: add `permdisplay(equaltailed|absolute|left|right|raw)`, keeping the
four flags as undocumented synonyms. Cost: renames a surface shipped 2026-09-06; `permtci` is genuinely
additive and would stay a flag. (b) Unify on flags: split `bootdisplay()` into `bootbc`/`bootbca`/`bootwald`
— worse, the four bootstrap families really are one mutually exclusive display choice. (c) Keep both
grammars, add `bootside(left|right)` to surface `e(boot_ptails)`, document that matrix, and add one help
paragraph mapping the two grammars onto each other.

**Defensible?** No.

**Verification.** Confirmed (one-line offset on the `bootdisplay` inlist: 6148). The verifier added that
the bootstrap is not free of the add-one convention, only unable to select it — the coefficient p's have no
add-one while the joint family does (CS-1) — which strengthens pathway (c), since no relabelling of options
fixes it. It also noted an internal doc split: sthlp:475 lists five values of `e(perm_pconv)` and sthlp:1073
lists two (CS-119).

#### CS-22 — `boottstat` and `permtstat` are parallel names over very different gate surfaces, and `permtci` turns `permtstat` on for you

**X.** `boottstat` carries a wall of gates, all five verified in the current tree: refused with `poisson`
(ado:6139-6142), with `bootstrap(bayes)` (6143-6146), with `bootdisplay(wald)` (6152-6159), with
`bootdisplay(bc)`/`(bca)` (6160-6169, the F003 withdrawal shipped 2026-09-06), and with `fastnose` via
bootstrap's own refusal (6184-6187). With `cluster()` it forces the physical resample-and-re-identify leg and
therefore requires the cluster variable to nest `i` (6223, 6243), which sthlp:485 says roughly doubles the
cost of each replication. It also accepts two companions whose effect is confined to stored results:
`bootdf()` "has no effect on any reported p-value or interval ... and it only aligns `e(df_r)`" (sthlp:503)
and `ties` "touches only the stored coefficient-draw bias-corrected tables, never anything displayed"
(sthlp:495).

**Y.** `permtstat` has no engine gates at all. Its only screen is `fastnose` (5799-5802) plus a caution when
`cluster()` equals `i` (5902). Its side effect runs the opposite way from boottstat's gates: it suppresses
permute's fast-mode default, so the run gets slower for a different reason (5404). And `permtci` implies
`permtstat` (117-118, 4980) with no runtime notice *of the implication* — the help says "Implies permtstat"
(sthlp:457) and the header and footer name `permtstat` (7653, 10025), but nothing says `permtci` caused it
or that the fast path was dropped as a result. There is no bootstrap option that turns `boottstat` on for
you.

**Why a user stumbles.** The names promise parallel objects. A user who learned `boottstat` as a heavily
constrained mode assumes `permtstat` is similarly gated and does not try it; a user who reads `permtci` as
"one more interval" does not realise it re-runs the entire job on the slow inference path.

**Pathways.** (a) Make the side effects explicit: extend the mode note that already prints under
`permtstat` (10025) to say that `permtci` implied it and that the fast-mode default was dropped, with its
cost; optionally require `permtstat` to be typed rather than implied, which breaks existing `permtci`-only
calls. (b) Refuse `bootdf()` and `ties` under `boottstat`, or label them as stored-results-only, so no
accepted option is display-inert — consistent with how `bootdisplay(wald)`/`(bc)` are handled; breaks calls
that pass `bootdf()` habitually. (c) Publish a two-column "studentized inference" table in the help listing,
per engine, the gates, the side effects and the cost.

**Defensible?** Yes for the gates themselves. The poisson and bayes refusals are principled rather than an
oversight: bootstrap replications carry duplication or Dirichlet weights that make `se_r` a different
functional of the data than `se_obs`, while permutation replications relabel onset without reweighting rows,
so that failure mode cannot arise under `permtstat`. `fastnose` is the only gate the two share, and it fires
for the same reason in both.

**Verification.** Partly confirmed. Corrections applied: `permtci`'s implication is documented and the run
is labelled, so "silently" overstates it — the gap is that nothing attributes the mode to `permtci` or names
the dropped fast path; and `bootdf()`/`ties` are *display*-inert, not inert, so pathway (b) is reframed.

#### CS-23 — `fastnose` is refused by `bootstrap()`, is `permute()`'s own default, and typing it explicitly can break a `permute()` run that would otherwise be legal

**X.** `bootstrap()` refuses it outright: ado:6184 (see CS-3).

**Y.** `permute()` sets it itself at ado:5404, and that assignment is placed deliberately **after** the
`fastnose` refusal block at 5367-5391, whose guard is a plain non-empty test on the **user's** token; the
comment at 5393-5403 says the ordering is deliberate. The consequence is that the same effective state gets
opposite verdicts depending on who asked for it:

    wooldid y i t ttre, cluster(i) permute(shuffle) permreps(99) esfixedbaseperiod espre(2) espost(2) adversarial

is legal and reports a permutation adversarial p, while adding the redundant word `fastnose` to that call is
r(198), "fastnose is not compatible with adversarial" (5376-5379). The same holds for `jointtests`
(5380-5383).

**Why a user stumbles.** The natural way to say "I do not need analytic standard errors" is `fastnose`.
Under `bootstrap()` it is an error; under `permute()` it is redundant; under permute-plus-adversarial it is
an error for a reason that no longer applies, because the engine supplies a permutation substitute for
exactly those surfaces. The `fastnose` help entry has to be read as two mutually qualifying halves
(sthlp:266-267). The deliberate ordering makes the trap sharper, not softer: what was considered was
permute's own assignment not tripping the refusals, not the user typing both.

**Adjacent, and worse.** `fastnose` typed together with `permanalyticse` hits no gate anywhere — the 5404
assignment does not fire because `permanalyticse` is non-empty, the user's `fastnose` stands, and fast mode
silently overrides an explicit request for analytic standard errors.

**Pathways.** (a) Make the user's `fastnose` a no-op under `permute()`: move the refusal block after the
permute assignment, or exempt permute runs from the `adversarial`/`jointtests` clauses. A
`fastnose`+`permute`+`unconditionalse` combination would then be caught by a different message; low cost.
Extend it to refuse or resolve `fastnose` + `permanalyticse`. (b) Refuse `fastnose` under `permute()` the
way the bootstrap does, on the ground that the engine owns the mode: one more refusal to learn, nothing
statistical lost, and the rule becomes "the engines choose fast mode, you never do". (c) Keep and document
the ordering in the `fastnose` entry with the exact legal/illegal pair.

**Defensible?** No.

**Verification.** Confirmed by code path and help text; the return codes were not run (the process cap never
cleared), but the path is deterministic — the guard at 5367 is a plain non-empty test, the adversarial branch
errors 198, and 5404 is unreachable before it. Corrections applied: the two mutually qualifying help
sentences are at sthlp:266-267, not 259. The verifier also found that the help's refusal list at sthlp:267
omits `pretrendpower()`, which the code refuses at 5384-5387, and the `fastnose`+`permanalyticse` hole above.

#### CS-24 — an option specified without the option it modifies gets one of three treatments, with no visible rule

**Tier 1, hard error (r198).** `lattice()` without `contreat()` (ado:4905-4909); `makecfxplots` without
`contreat()` (5252-5256); `contreatwithin` without `contreat()` (5409-5413); `poisexpresults` without
`poisson` (5433-5437); `altmethod_adversarial` without `adversarial` (5445-5451, whose comment names "the
modifier-without-parent convention"); every `perm*` companion without `permute()` (4965-4974); every
bootstrap-family companion without `bootstrap()` (6029-6045) — a list that includes the **bare-named**
`reps()` and `ties` as well as the `boot*` ones; `bootalpha()` with `bootstrap(standard)` (6092-6095);
`permprobs()` outside `permute(probabilities)` (5840); `permstrata()` outside `permute(shuffle)`
(5959-5962).

**Tier 2, caution and ignore.** `cfxlattice()` without `makecfxplots` — "Caution: cfxlattice() has no effect
without makecfxplots and is being ignored" (5290-5294) — and `esrelativeto()` without `esfixedbaseperiod`
(5355-5359), same shape.

**Tier 3, accepted with no message, effect confined to stored results or to nothing.** `contreatpoly()` and
`contreatelasticitytype()` without `contreat()`: consumed only inside `if "`contreat'" != ""` blocks
(1264-1275, 2150-2152, 2172, 2550) and value-screened only at 5417-5426, so `contreatpoly(3)` alone is rc 0
and silent while `contreatpoly(0)` alone errors — an orphan that errors on its *value* but not on its
orphanhood. Also `latticeignoreweights` without `lattice()` (used only at 1195); `cfxplottypes()`,
`cfxplotnose` and `makecfxplotsbysg` without `makecfxplots` (the parse block is gated at 5174); `replace`
with no `saveplots()`/`bootsave()`/`permsave()` (no screen exists); `ties` under the default
`bootdisplay(percentile)` and under `boottstat`, which still feed the stored bias-corrected tables
(`e(boot_ci_bc)`, `e(boot_z0)`, `e(boot_ties)` at 9148) though nothing displayed moves; and `bootdf()` under
`boottstat`, which still sets `e(df_r)`.

**Why a user stumbles.** A user who meets the `lattice()`-without-`contreat()` error concludes wooldid
catches orphan modifiers, and trusts that `contreatpoly(3)` without `contreat()` would have been caught too.
The three tiers cut across families: within the continuous-treatment block alone, `lattice()` errors,
`cfxlattice()` cautions, and `contreatpoly()`/`latticeignoreweights`/`cfxplotnose` are silent.

**Pathways.** (a) Unify on error for every orphan modifier — one screen block, one message template,
mechanical to write and test. Cost: breaks do-files carrying harmless leftovers; and `ties`/`bootdf()` under
`boottstat` are *documented* as display-inert, so making them errors contradicts published help text rather
than merely breaking scripts. (b) Unify on caution-and-ignore. Cost: loses the loud protection for orphans
that genuinely change the answer (`contreatwithin` also disables estimability checks — comment at
5409-5410). (c) Keep the tiers but publish the rule — "an orphan that would have changed the answer is an
error; an orphan that is merely inert draws a caution" — and move the whole silent set into the caution tier,
so no accepted option is ever a silent no-op.

**Defensible?** No as a set; several sites are individually defended (see CS-78, which quotes the in-code
rationales).

**Verification.** Partly confirmed; every quoted line and message checks out within a few lines.
Corrections applied: two tier-3 examples were wrong or overstated — `emptycellsoverride` is not a modifier of
`verbose` and is not inert (it suppresses `set emptycells drop` unconditionally at 549-556, a real change to
the reghdfe design matrix), and `ties` under `percentile` changes stored results rather than nothing. `reps()`
and `ties` should be named in tier 1 as bare-named members of the bootstrap companion screen. Anchor for
CS-43, CS-78, CS-111.

#### CS-25 — `safety(off)` means seven different things, one of which happens without the user typing it

**Four leave the computation alone.** (1) An error becomes a caution and the run proceeds unchanged:
`jointtests` with analytic clustered SEs at or above `i` (ado:1523/1527); `bootstrap(bayes)` with a
cell-weighted estimand (6217/6220); a subgroup variable that splits a cohort or varies within `i`-by-`t`
(5938/5947). (2) A domain check is bypassed and the value used as given: the `lattice()`/`cfxlattice()`
non-negative-integer checks (sthlp:752, 791), `customweightmultiplier()` negatives (sthlp:926), the `ttre`
within-`i` constancy and grid checks (sthlp:66). (3) The caution is suppressed and the behaviour is
identical: `espostlength` 0 raised to 1 under `subgroups()` — the raise happens either way, `safety(off)`
only silences the three lines (4936-4941). (4) Not honoured at all: the `fastcollapse` ineligibility
refusals — sthlp:328, "safety(off) does **not** override that refusal, because each entry below is a
demonstrated wrong answer rather than an accepted risk".

**Three change what the run does.** (5) A different number is computed: the cell-collapse degenerate-subgroup
screen — under default safety wooldid restores the full data and recomputes from that point; under
`safety(off)` it keeps the collapsed rows and cautions (5696-5721; sthlp:362). (6) Requested output is
withheld unless `safety(off)` is given: `makecfxplots` with `att_it`/`att_i`/`att_itime` — under default
safety the plots are simply not produced, with a caution reading "wooldid will not produce these plots
unless you specify safety(off)" (4913-4921; the production gate is 2655). (7) The user is put into the off
regime without asking: `contreatwithin` rewrites the working safety value to off (5157-5162; sthlp:739).

**The two-macro wrinkle.** The `contreatwithin` rewrite touches only the macro `safety`; `safetyinput`,
captured at 4989, keeps the user's own value, and cases 1, 2 and 5 all gate on `safetyinput` (the collapse
screen is invoked as `safety(`safetyinput')` at 3094/3303/3321/3364/3419). So a `contreatwithin` run reaches
case 6 only — and does so incoherently, because the caution at 4913-4921 runs *before* the rewrite, so the
run prints "wooldid will not produce these plots unless you specify safety(off)" and then draws them
(CS-98). The ado documents the split as deliberate (comments at 1490, 5228, 5273, 5427), which is arguably an
eighth wrinkle: the same word means two different things inside the code depending on which macro a site
reads.

**Why a user stumbles.** `safety(off)` reads as one dial with one meaning: "stop stopping me". It is at
least seven, and two of them affect numbers rather than messages. The `makecfxplots` case is the sharpest —
a user asks for a plot, gets no error and no plot, and the only way to obtain it is an option whose
documented purpose is to disable safety checks. The `contreatwithin` case fires with no user gesture at all.

**Pathways.** (a) Split into named switches with `safety(off)` retained as "all of them": e.g.
`safety(nochecks)`, `safety(nofallback)`, `safety(quiet)`. Cost: a new surface plus a migration path; the
seven cases need naming. (b) Keep one switch and classify every site in the `safety()` help entry as a table
— refusal lifted / check bypassed / caution silenced / different number / output enabled / not honoured /
set for you by `contreatwithin` — with a column for which macro each site gates on. Cheap, breaks nothing,
and is the smallest thing that makes the option honest. (c) Narrow `safety(off)` to "lift refusals and
bypass checks only" and give the two behaviour-changing sites their own options (a collapse-guard switch;
letting `makecfxplots` draw non-`att` plots with a caution). Two behaviour changes, and the most defensible
end state.

**Defensible?** No.

**Verification.** Partly confirmed; six of the seven sites verified verbatim. Correction applied: the
finding said `contreatwithin` lands the user in cases 5 and 6 — because of the `safety`/`safetyinput` split
it reaches case 6 only. See CS-40 (the same option from the policy lens, with a different taxonomy), CS-37,
CS-98, CS-103.

#### CS-26 — the statistic column is marked per engine and the p-value column is not, and `se` carries two different quantities under one label

**X.** `permute()` renames three of the six result columns — `se`→`se_analytic`, `t`/`z`→`t_perm`,
`p`→`p_perm` (ado:7504-7507; subgroup tables 7543-7546; pairwise contrasts 7992-7995). A script reading
`e(wdidmainresults)` colnames can therefore detect permutation inference. The `z`→`t_perm` substitution means
a poisson run's `z` column is relabelled `t_perm`.

**Y.** `bootstrap()` renames the statistic but not the p-value: `t`/`z` becomes `t_wald`/`z_wald` under
percentile/bc/bca (8855, F053) and `t_pt` under `boottstat` (8911-8913), but the p column stays plain `p` in
every family — even though it is an order-statistic rank, not the tail probability of the statistic printed
beside it. The standard-error column stays `se` and holds the bootstrap standard deviation (8808-8836),
except under `boottstat` where it becomes `se_analytic`. So on the same six-column stripe: `se` means the
analytic delta-method SE on a plain run (3849) and the bootstrap standard deviation on a bootstrap run, both
under the label `se`; `p` means a t/z tail probability on a plain run and a rank on a bootstrap run, both
under the label `p`; while `p_perm`, `t_wald` and `t_pt` do mark themselves.

**Why a user stumbles.** Anyone post-processing `e(wdidmainresults)` by column name can detect permutation
inference and studentization from the names but cannot detect percentile-rank inference, and the two
analytic-looking labels are precisely the two carrying non-analytic numbers on a bootstrap run. A user who
learned from `p_perm` that wooldid marks non-analytic p-values will trust a bootstrap table's `p` column.

**Pathways.** (a) Unify on marking: `p_pctl`/`p_bc`/`p_bca`/`p_wald` and `se_boot`, matching
`p_perm`/`se_analytic`. Cost: breaks every script reading bootstrap tables by column name, and the hidden
family twins share the stripe so all of them move together. (b) Unify on not marking: revert
`t_wald`/`t_perm`/`p_perm` to `t`/`p` and record the family only in `e(boot_display)`/`e(perm_pconv)`. Cost:
reverses F053 and removes the in-table warning the author just added. (c) Keep and publish one help table
giving the colname stripe per engine and per display family — cheapest, and it also serves the `e()`
documentation.

**Defensible?** No.

**Verification.** Partly confirmed. Corrections applied: "`se` means three different quantities" is two
under the label `se` (the third is separately marked `se_analytic`); and "cannot detect percentile-rank
inference" holds at the colname level only — `e(boot_display)` and `e(boot_stat)` are posted and documented.
Line numbers on the permute side run one lower than filed; the bootstrap `se` evidence is at 8808-8836, not
the 8797-8803 ruling comment (a stale comment block just above at 8785-8789 still describes the withdrawn
interval-inverted SE and contradicts the code — a small separate item). See CS-7, CS-74.

#### CS-27 — `cluster()` names the analytic variance, the resampling unit, the randomization unit or the IM aggregation unit depending on the feature, and its requirements change with each

**X.** On a plain analytic run `cluster()` is optional and may name several variables (multi-way). Under
`impvalues` it must be a single variable and becomes the aggregation unit for the Ibragimov-Muller average
(ado:4740-4746; `margins ... over(`cluster')` at 2303). Under `bootstrap()` it is still optional — without
it the resampling unit is the observation, and `bootdf` then defaults to normal — but multi-way is refused
(6206, "bootstrap() supports a single (numeric) cluster() variable; multiway clustering is not available
with bootstrap-based standard errors"), a cell-weighted estimand without `cluster()` is refused (6231),
`boottstat` with `cluster()` requires the variable to nest `i` (6242, r459), and `bootdf(clusters)` requires
`cluster()` to exist at all (6106-6108).

**Y.** Under `permute()` `cluster()` is mandatory and stops describing the variance: it is the experimental
unit whose treatment onset is reassigned. It must be exactly one variable (5852-5855 — the same check does
double duty as the mandatory check, so a user who omits `cluster()` entirely gets a message that reads as if
they supplied too many), must nest `i` (5903-5912, r459, with wording entirely different from the
bootstrap's), must carry a constant onset within itself (~5885), and must code never-treated as plain missing
(~5870). On a default `permute()` run `cluster()` does not enter an analytic variance at all, because the
observed run is a fast (no-SE) run — though it is still forwarded into the replication options (127). Under
`permtstat` or `permanalyticse` it resumes its variance meaning *simultaneously* with its randomization
meaning, which is arguably a fifth state.

**Why a user stumbles.** One option name carries four meanings. A user with a two-way clustered analytic
specification who adds `permute()` is told to pick a single variable, with nothing to signal that the
variable now defines the experiment rather than the variance — and that clustering finer than `i`, perfectly
reasonable for a variance, is now a design error. A user who wanted to cluster the variance on state and
resample on region cannot express it at all.

**Pathways.** (a) Add distinct option names — `bootcluster()`/`permunit()` — defaulting to `cluster()`. Cost:
a new surface, but it makes the overload visible and lets one run cluster analytically on one variable while
resampling or randomizing on another, a real capability today's design cannot express. (b) Keep the overload
and add one help paragraph, "what cluster() means under each engine", plus a footer line on every engine run
naming the resampling/randomization unit and its count — cheap, and the count is already computed. (c)
Require the engines to restate it in their own argument, e.g. `permute(shuffle, unit(state))` — most
explicit, most disruptive.

**Defensible?** Yes — each meaning is right for its feature, and the requirement differences follow from the
mechanics. The overload and the divergent error wording are what trip users.

**Verification.** Confirmed; every cited line and quotation verified.

#### CS-28 — `level()`, `dots()` and `rseed()` are engine-neutral in the code and filed under bootstrap in the help

**X.** The options-by-purpose synopsis (sthlp:72) reads "... `permrawp` `bootstrap(mode)` `reps(int)` `ties`
`level(#)` `dots(int)` `rseed(#)` `bootalpha(spec)` ..." — the three shared options sit inside the bootstrap
block, between `bootstrap()` and `bootalpha()`, after the whole `perm*` family has ended.

**Y.** In the code `level()` applies to every run since F054 (ado:526; sthlp:489 lists the analytic tables,
the plots, the `gteffects` cells and the Roth slope alongside the engine intervals), and `dots()`/`rseed()`
are accepted under either engine (6033-6034 exempt them from the stray screen when `permute()` is present).
The only place the interface states the sharing rule is the bootstrap stray error itself (6042: "dots() and
rseed() are also accepted with permute(); level() applies to every run"). The permutation stray error (4972)
says only "require(s) the permute() option" and teaches nothing.

**Why a user stumbles.** A `permute()` user scanning the option list sees no seed, no dots and no level in
the permutation block and concludes randomization inference is not reproducible or is fixed at 95%. The
correct information is reachable only by triggering an error.

**Pathways.** (a) Split the three into their own synopsis line or dlgtab ("shared inference companions"),
referenced from both engine entries — cheap, breaks nothing. (b) Repeat them in both blocks of the synopsis
— duplication, but the synopsis is a lookup table, not prose. (c) Leave the layout, add one sentence to
`permute()`'s entry naming the three, and give the permutation stray message the same explanatory tail the
bootstrap one has. *Not verified (cap).*

#### CS-29 — `bootsave()`/`permsave()` validate their file names before the replication loop; `saveplots()` validates nothing and fails at plot time

**X.** `bootsave()` and `permsave()` reject a directory part and refuse to overwrite an existing `.dta`
without `replace`, checked in the option-screen programs before any estimation work (ado:6070-6086 and
5786-5798) and again in the writers (6314-6322, 8098-8107). The comment at 6072 records why: "the name rules
... are checked here BEFORE the replication loop runs. Until 2026-09-05 only wooldid_bootsave checked them,
after the loop, so a name clash cost the whole run."

**Y.** `saveplots()` has no pre-check at all. The same `replace` token is handed straight to `graph save` and
`save` at plot time (2727, 2771, 2820, 2860, 10494, and 10504 for the `_cfx.dta`), uncapped. A
`saveplots()` run whose `.gph` or `_cfx.dta` already exists and that omits `replace` aborts with r(602) at
the plotting stage — after estimation and, under `bootstrap()`/`permute()`, after every replication has run.
`saveplots()` is declared `name` in the syntax (506), so it cannot carry a directory part either, but
nothing says so.

**Why a user stumbles.** One word, `replace`, guards three save features with opposite failure timing. The
lesson explicitly learned for `bootsave()` five days before this review was not carried across to the older
feature, so the expensive path is the unprotected one.

**Pathways.** (a) Pre-check the `saveplots()` targets the same way; the plot-name list is derivable from the
estimand, subgroup and `cfxplottypes()` sets before estimation — a moderate amount of name-construction code
duplicated from the plotting section. (b) Wrap the graph saves in `cap` and print a caution instead of
aborting, so a name clash costs a file rather than the run — cheap, but a silent non-write becomes possible.
(c) Document the difference in the `saveplots()` and `replace` entries. *Not verified (cap).*

#### CS-30 — 25 of 83 options have no abbreviation, and three prefixes resolve to the option a user is least likely to mean

**X.** The syntax line (ado:506, duplicated byte-identically at 4670) is the abbreviation spec. It gives
`fastnose` the minimum `fast`, while `fastcollapse` and `fullprecisionslow` have no abbreviation at all — so
`wooldid y i t ttre, fast` silently means "compute no standard errors", not "collapse the cells".
`PERMTci` has minimum `permt` and `PERMTStat` has minimum `permts`, so `permt` yields `permtci`, which
implies `permtstat` (4980) *and* adds an inverted-|t| interval; the user who wanted only the studentized
statistic must type `permts`. `PERMSide` has minimum `perms`, next to `permstr(ata)` and `permsa(ve)`.

**Y.** The 25 options with no abbreviation are not a natural class: `reps`, `ties`, `dots`, `rseed`,
`permreps`, `fastcollapse`, `fullprecisionslow`, `contreat`, `lattice`, `cfxlattice`, `customvce`, `safety`,
`regtol`, `replace`, `update`, `oldsyntax`, `emptycellsoverride`, `altmethod_adversarial`,
`makecfxplotsbysg`, `cfxplotnose`, `fe`, `ccc_absorb` and the four `att` spellings — while every other
`boot*`, `perm*`, `contreat*`, `es*` and `cfx*` option abbreviates. Within one family:
`bootstrap`/`bootalpha`/`bootdf`/`bootdisplay`/`boottstat`/`bootsave` all abbreviate, their two bare-named
companions `reps` and `ties` do not. Separately, the synopsis (sthlp:72) prints `permside` and `permrawp`
without abbreviation markers while the code gives them `perms` and `permraw`.

**Why a user stumbles.** All three collisions parse successfully, so nothing warns. `fast` is what a user
reaching for the cell collapse would try first; `permt` is what they would type for `permtstat`. The missing
abbreviations are unpredictable enough that users type everything in full, which is the outcome the scheme
was meant to avoid.

**Pathways.** (a) Give every option an explicit minimum chosen to avoid the collisions — `FASTNOse` so
`fast` becomes illegal and `fastc` can go to `FASTCollapse`; `PERMTCi` so `permt` becomes illegal. A
breaking change for existing do-files, but a loud one: illegal abbreviations error rather than silently
switching meaning. (b) Leave the abbreviations, add a "watch out" note naming `fast`, `permt` and `perms`,
and fix the synopsis markers for `permside`/`permrawp` — cheap, does not remove the trap. (c) Rename the
colliding options themselves — most disruptive, least likely to be worth it. *Not verified (cap):* parsed
from the syntax line's capitalization; the empirical parse test was written and never run.

#### CS-31 — `bootdf()`'s default flips on `cluster()`, and the message that teaches the option says the default is normal

**X.** The resolved default is context-dependent: ado:392, `if bootdf == "" local bootdf = cond(cluster !=
"", "clusters", "normal")`. The help agrees (sthlp:503): "The default depends on whether you clustered:
bootdf(clusters) when cluster() was given, and bootdf(normal) when it was not."

**Y.** The error a user meets on a typo says otherwise (6117): "bootdf() must be bootdf(normal) - standard
normal critical values, **the default** and the convention Stata's own bootstrap prefix uses - or
bootdf(clusters) ...". The block comment above it repeats the same wrong claim (6097-6098). The neighbouring
refusal has it right (6107): "Use bootdf(normal) - the default without cluster()".

**Why a user stumbles.** Two of the three places a user meets this option assert an unconditional default
that is wrong on every clustered run, and the difference is a t(G-1) versus a z reference distribution — it
moves the printed p-value and the Wald interval. A user who mistypes `bootdf()` once carries away the wrong
model of the default.

**Pathways.** (a) Fix the message — one line, breaks nothing. (b) Also print the resolved value once in the
footer beside the reference-distribution sentence (10101), so the default is never inferred. (c) Drop the
context-dependent default and require `bootdf()` explicitly on clustered runs — a newly required option, and
it reverses a deliberate, well-reasoned choice. *Not verified (cap).*

#### CS-32 — the two engines disagree in both directions about when a failed replication deserves a caution and when a run should stop

**X.** `bootstrap()` cautions whenever any replication fails, with no threshold (ado:8614, "Caution: `nfail'
of `reps' bootstrap replications failed to estimate and were excluded") and refuses only when zero completed
(8610-8612).

**Y.** `permute()` cautions only above 10% failures (7372) — though any failure under `permexact` draws a
caution, because the exact-enumeration reading breaks (7371) — and refuses whenever fewer than
max(20, permreps/2) completed (7364).

**Why a user stumbles.** A bootstrap run losing 3 of 1000 replications prints a caution; a permutation run
losing 90 of 1000 does not. A permutation run with 19 of 40 completing is refused outright; a bootstrap run
with 1 of 1000 completing proceeds and prints p-values of resolution 2. The user is warned loudest where the
problem is smallest.

**Pathways.** (a) Unify on permute's rules: add a completed-replication floor and a 10% caution threshold to
the bootstrap. Cost: bootstrap runs that legitimately lose replications (a resample with no comparison group
is normal with few clusters) would now stop, so the floor probably wants to be lower than max(20, reps/2).
(b) Unify on the bootstrap's: caution on any permutation failure and refuse only at zero — removes a
deliberate randomization-inference safeguard. (c) Keep both and state both rules in one help table;
independently, make the bootstrap caution threshold-based so it stops crying wolf on 3-in-1000. *Not verified
(cap).* See CS-45, CS-90.

#### CS-33 — parallel `e()` objects across the engines are named and shaped differently, and one family's option name and `e()` name are different words

**X.** `bootstrap()` posts `e(bootdist)` as the full coefficient-draw matrix (capped at 10,000 completed
replications, ado:9220-9227), with `e(bootdist_se)` and `e(bootdist_pwse)`; the displayed family is in
`e(boot_display)` (9157); the interval families are `e(boot_ci_normal)`, `e(boot_ci_percentile)`,
`e(boot_ci_bc)`, `e(boot_ci_bca)`.

**Y.** `permute()` posts `e(permdist)`, but it is **not** the twin of `e(bootdist)`: it holds only the legacy
pooled estimands (224), while the full coefficient-draw matrix is `e(permdist_b)` with `e(permdist_se)`
(7636-7640). The convention in force is `e(perm_pconv)`, not a `perm_display`. Separately, the family the
option calls `bootdisplay(wald)` is stored as `e(boot_ci_normal)` (9122) and named that way in the footer
(10107), so the option name and the `e()` name for one object are different words.

**Why a user stumbles.** A script looping over both engines by `e()` name gets the wrong matrix under
`permute()` — a narrow legacy matrix where it expected the full draw matrix — silently, with no error and no
dimension mismatch it can rely on. A user grepping `ereturn list` for "wald" after asking for
`bootdisplay(wald)` finds nothing.

**Pathways.** (a) Add aliases: post `e(bootdist_b)` beside `e(bootdist)`, `e(perm_display)` beside
`e(perm_pconv)`, and rename `e(boot_ci_normal)` to `e(boot_ci_wald)` keeping the old name hidden — low cost,
no breakage. (b) Rename the option to `bootdisplay(normal)` — breaks do-files and the help's several
`bootdisplay(wald)` references. (c) Document the mapping in the stored-results section as a two-column table
of parallel names, with an explicit warning that `e(permdist)` is not `e(bootdist)`'s twin. *Not verified
(cap).* See CS-58 (verified), CS-63, CS-125.

#### CS-34 — `espostlength` quietly corrects its own argument and has an undocumented magic value, where every other numeric option refuses out-of-domain values

**X.** `espostlength(0)` with `subgroups()` and `esprelength() > 0` is silently raised to 1 (ado:4936-4941),
with a three-line caution under default safety and no message at all under `safety(off)`. The escape hatch is
`espostlength(-999)`, taught only by the third line of that caution and rejected as a value elsewhere: the
syntax diagram says `espostlength(int>=0)` (sthlp:72), the entry (sthlp:577) documents the raise but not
`-999`, the code errors on any negative except `-999` (4929-4930), and a grep for `-999` in `wooldid.sthlp`
returns nothing.

**Y.** Every other out-of-domain numeric argument is a hard error with the remedy named and no correction:
`permreps(10)` (5807-5810), `contreatpoly(0)` (5423), `esprelength(-1)` (4923), `esrelativeto(0)` under
`esfixedbaseperiod` (5441), `level(5)` on an engine run (5780/6062), `rseed` out of range (6019), `dots()`
negative (6068), `bootdf()` non-integer (6117).

**Why a user stumbles.** A user who has learned that wooldid refuses out-of-range numbers will not expect one
option to quietly change its argument, and cannot find the override anywhere in the help — only by reading a
caution that `safety(off)` suppresses. A `safety(off)` run therefore silently estimates a post-treatment
horizon the user asked not to have.

**Pathways.** (a) Document `-999` in the `espostlength()` entry, or replace it with a named option
(`esnopost` or similar) — cheapest and removes the worst of it. (b) Turn the raise into a refusal that names
the remedy, matching every other numeric screen — breaks calls that rely on the silent raise, which are
exactly the calls the hack was written for. (c) Keep the raise but print the caution regardless of `safety()`,
so the correction is never silent; `safety(off)` should not hide a change to the estimand. *Not verified
(cap).* See CS-40 sense (iii), CS-103.

#### CS-35 — "every family is computed regardless" holds for three of four, and `boottstat` makes three of four undisplayable

Overlaps **CS-76** (footer, confirmed) and **CS-115** (help self-contradiction, confirmed), whose citations
supersede this item's. **Adds the second direction:** `boottstat` makes three of the four `bootdisplay()`
values illegal — `wald` and `bc`/`bca` are refused (ado:6155-6169) and only `percentile` remains — so the
same option is an unconstrained display switch on one path and a single-valued one on another; and under
`bootdisplay(bca)` the `gteffects`/`e(cfx*)` cells silently carry the bias-corrected family instead, with a
caution (6359). Pathways: (a) always compute BCa — honest, and it costs G extra fast fits on every clustered
run, which the author already judged too expensive; (b) qualify the sentence everywhere it appears and post
`e(boot_bca_available)` so the absence is discoverable from `e()`; (c) add a `bootfamilies(all)` switch so
the full promise is opt-in. Defensible. *Not verified (cap)* for this framing. See CS-11, CS-109.

#### CS-36 — `level()` is now universal but is range-checked only on the engine paths

**X.** Both engines validate the range themselves: ado:5780 in `wooldid_permctl` and 6062 in
`wooldid_bootctl`, each with "Error: level() must be between 10 and 99.99."

**Y.** An analytic run reaches neither check. `wooldid_bootctl` exits at its stray-companion screen (6046)
before the level check whenever `bootstrap()` is empty, and `level()` is not in the stray list;
`wooldid_permctl` only runs when `permute()` is non-empty (4981). So on a plain analytic run `dots(-5)` is
accepted silently and `level(5)` is never seen by wooldid at all — it reaches `set level` in the worker (526)
and fails with Stata's own message rather than wooldid's. The range rule is stated twice in the code and
enforced on neither of the two paths F054 newly made it meaningful for. Note also that `level(0)` is the
sentinel meaning "not specified" (the guard at 6062 is `if `level' != 0 & ...`), so `level(0)` is silently
accepted as 95 on every path.

**Why a user stumbles.** Minor, but F054 made `level()` apply to analytic tables, so a user now naturally
types `level(90)` on a plain call — and a typo gives them Stata's bare `set level` error, from a line inside
the worker, with no mention of wooldid or of the option they typed. The same typo under `bootstrap()`
produces the friendly wooldid message.

**Pathways.** (a) Move the `level()` and `dots()` range checks into `wooldid_validate` ahead of both control
programs and delete the duplicates — small and mechanical; `validate` has a stated compiled-size budget, and
these two checks are cheap. (b) Move only the `exit 0` in `bootctl` to after 6068, so the shared options are
validated on every call while the bootstrap-specific ones stay behind the master gate — one line moved, keeps
the file layout and the check-sync test intact. (c) Document that `level()` is validated by Stata's own `set
level` on analytic runs. *Not verified (cap).* See CS-2, CS-46, CS-62.

### 2.3 Policy: refuse, caution, skip

#### CS-37 — `contreatwithin` silently flips the working safety flag to off, and the casualty the help does not predict is the maxvar pre-check

*Filed high; verifier: medium, with the title narrowed — see the verification line.*

**X.** A user who runs `wooldid ..., contreat(d) contreatcontrols(byCohort) contreatwithin` still believes
safety is on, having typed no `safety(off)`. The help's `contreatwithin` entry (sthlp:739) tells them what is
relaxed: "This option disables some of *safety*'s features (the estimability-relaxation checks; the
data-integrity screens, such as onset constancy within units, still run unless *safety(off)* is given
explicitly)."

**Y.** The code sets the working flag wholesale — `local safety off` at ado:5167, unconditionally at the end
of the `contreatwithin` block — and the flipped value reaches the worker through `wooldid_validate`'s
blanket `c_local` return loop (5459). Six worker gates key on the flipped value: 852, 920, 983, 1030 (which
additionally excludes `contreatwithin` in its own condition), 1446 and 1762. Nothing announces the flip.
The screens hardened against it gate on `safetyinput`, captured at 4989 before the rewrite, and keep running:
`customweightmultiplier` domain (654), onset constancy (698), the subgroups screens (1098, 1122), and the
`jointtests` cluster-at-`i` refusal (1521, whose comment at 1489-1491 says in so many words "so it gates on
safetyinput, not safety (contreatwithin overwrites the latter internally)"). `wooldid_fcdeg` and
`wooldid_bootctl` are likewise passed `safetyinput`.

**Where the help's own rule holds, and where it does not.** The multicollinearity scan at 1762 (hard error at
1816-1822) *is* an estimability check on any reading, and its own error text at 1820 prescribes
`contreatwithin` as the remedy — so its going dark under `contreatwithin` is the documented and intended
half, and pathway (a) below must exempt it or create a loop in which the error tells the user to specify the
option that suppresses the error. The genuinely undisclosed casualty is the gate at 1446: the granular-`i`
caution (1462) and the maxvar r(907) preflight (1464-1467). Those are **resource** checks, not
estimability relaxations, so the help's stated rule predicts they still run; both sthlp:~985 and the in-code
messages say only "(safety(off) suppresses this check.)", and `contreatwithin` is never named as a
suppressor.

**Why a user stumbles.** The user reads that `safety(off)` is what turns off the maxvar pre-check, never
types it, and can blow past that pre-check into a raw r(907) mid-run — the failure the check exists to
prevent. And "some of safety's features" reads as "the ones about control groups", given the sentence it sits
in.

**Pathways.** (a) Unify on `safetyinput`: make every screen gate on the user's own value and give
`contreatwithin` a narrow, named exemption list — which **must** include the multicollinearity scan, for the
reason above. Cost: a handful of gate rewrites plus a battery re-bless. (b) Keep the flip and announce it:
one line at parse time naming the disabled screens, and enumerate them in the `contreatwithin` help entry
rather than saying "some". One display line, one help paragraph; breaks nothing, and leaves a silently
weaker run available. (c) Split the flag into `safetyestim` (relaxed by `contreatwithin`) and `safetyident`
(never relaxed), so the distinction the code already draws by hand becomes a named concept — the cleanest end
state, the largest diff.

**The minimal fix.** A one-token change at 1446, `safety` → `safetyinput`, plus naming the disabled screens
in the help entry.

**Defensible?** Partly. The multicollinearity relaxation is intended and documented; the 1446 resource check
is not.

**Verification.** Partly confirmed; every code and help citation checks out and the flip is silent. The
finding as filed over-scoped the defect — sthlp:739 does supply a rule, and it correctly predicts that 1762
is disabled. Title narrowed accordingly, priority dropped to medium, and the fix shrinks to one token. See
CS-98, which verified the *other* half of the same flip with a live run and stands at high.

#### CS-38 — `adversarial` cancels itself silently on the most common default; `pretrendpower()` refuses the same situation at parse time with a named remedy

**X.** `wooldid y i t ttre, att adversarial` — no `esprelength()`, no `espostlength()`, both defaulting to 0 —
runs to completion and prints nothing about the adversarial test at all. Line 1042,
`if (`espre' + `espl'+ 1 < 3 & "`timetrends'" != "") | `espl' == 0 | `espre' == 0 local adversarial`, clears
the option with no `dis` of any kind. Since `espre = min(min_t, esprelength)` (1008) and
`espl = min(max_t, espostlength)` (1009), a short panel cancels it too, again silently. Adversarial has two
further cancel sites that **do** speak (1048; 2445/2925, and the subgroup twins at 3278/3289) and a degrade
site (1055: without `esfixedbaseperiod` the line is fitted and reported without a p-value), so one option has
four cancellation behaviours, one of them mute.

**Y.** `pretrendpower()` screens the identical prerequisites at parse time, in `wooldid_ptpcheck` (called
from validate at 4799), and refuses with r(198) plus a named remedy: 9253-9254, "pretrendpower() requires
esfixedbaseperiod. Under the pooled reference the comparison period is a cohort-specific pool of pre-periods,
so the single omitted reference period that the Roth (2022) linear-trend normalization requires is not
defined. / Specify esfixedbaseperiod (with esrelativeto() if you want a reference other than -1)."; and 9258
for the `esprelength()` half. Note the domain difference too: `ptpcheck` keys on the user's `esprelength`
option, 1042 keys on the realized `espre`, so `pretrendpower()` errors on the option and `adversarial` folds
on the data.

**Why a user stumbles.** The single most likely first attempt at the feature — typing `adversarial` and
nothing else — produces a clean, complete, apparently successful run in which the requested diagnostic is
simply absent. The user hunts for `e(adversarial_pval_post)`, finds nothing, and cannot tell whether the
option is broken, unsupported on their data, or misspelled. The same user typing `pretrendpower(80)` alone is
told exactly what to add. The silence is run-time only, not undocumented: the `adversarial` help entry
(sthlp:~599-600) says the tests "will not be conducted if an insufficient number of relative time
period-specific treatment effects are requested" — but a user who read that still gets no confirmation that
the clause fired on their run.

**Pathways.** (a) Unify on refusal: screen adversarial's prerequisites in `wooldid_validate` as `ptpcheck`
does — "adversarial requires esprelength() and espostlength() of at least 1" — and keep the later
data-driven cancels as cautions. Small; breaks scripts carrying a decorative `adversarial` on runs where it
never did anything, which is arguably the point. (b) Unify on caution-and-degrade: give 1042 a message in the
house style, mirroring the help sentence above, and leave `pretrendpower()`'s parse errors alone, on the
ground that Roth's diagnostic is a headline number and the adversarial line a by-product. One line, and it
removes the silence, which is the actual harm. (c) Post a machine-readable flag,
`e(wdid_adv_cancelled)` with a reason string, alongside whichever message policy is chosen, so the absence is
visible on replay and in scripted pipelines.

**Defensible?** No.

**Verification.** Confirmed, by code and by live run
(`/Users/hegland/.claude/jobs/be9644d0/tmp/consist/skeptic/adv.log`). Run A, `att adversarial` on
`dta_binary`: completes normally, one ordinary `att` table, no message mentioning adversarial anywhere,
`e(adversarial_pval_att)` missing. Run B, `att pretrendpower(80)`: r(198) with the `esfixedbaseperiod`
message and its remedy line. Run C, `esprelength(3) espostlength(3) adversarial` without
`esfixedbaseperiod`: the 1055 caution does print. The verifier noted an adjacent, unverified oddity — in run
C the caution promises "the fitted line is reported without a p-value" and no adversarial line appeared in
the printed tables, so it presumably lands only in `e(adversarial_line_att)`.

#### CS-39 — the bootstrap cautions twice about cluster counts; `permute()` has no cluster-count caution at all

**X.** `bootstrap()` cautions when fewer than 10 clusters are in the estimation sample (ado:5977, inside the
cluster-count sanity block that starts at 5969, with a hard error at 5975 below 2) and again when fewer than
6 **treated** clusters are available for resampling (6279), citing MacKinnon and Webb (2017, 2018) and
closing "Consider permute() - randomization inference over the treatment assignment - or the wild cluster
bootstrap-t (boottest, on the stage-1 regression outside wooldid), which wooldid does not offer." Both are
pre-loop and keyed on counts of clusters; the F050 ruling is recorded at 6260-6274.

**Y.** `permute()` has **no** caution keyed on the number of clusters or treated clusters. Its cluster-count
gates are hard errors for degenerate designs — 5982 ("found `permnclus' cluster(s) across `permnlvl'
treatment-onset group(s); there is nothing to randomize") and 7156 (a design admitting exactly one
assignment). The cautions that exist are keyed on the number of **possible assignments**: `wdp_caut1` (mata
10980), `wdp_caut2` (10985/11086), `wdp_caut0` (11087). A design with 5 treated clusters and many
never-treated ones has an enormous assignment support, so none of them fires — and that is exactly the
setting the bootstrap's own caution at 6279 tells the user to bring to `permute()`.

**Two cautions that do exist, and why they do not reach this user.** Ado:5902 fires the exchangeability
caution, citing MacKinnon and Webb (2020) and closing "the test presumes exchangeable clusters in both
forms" — but only when `permtstat` is on **and** `cluster()` is the `i` variable, so a user who takes the
bootstrap's referral and runs a default permute (raw-coefficient statistic, cluster at any other level) sees
nothing. And the F036 resolution caution (5811-5818) is explicitly gated on `permexact` being empty, so under
`permexact`, where the floor 2/M is a direct function of how few treated clusters there are (sthlp:481 notes
one treated cluster among 12 floors the exact p at .167), no floor caution fires either.

**Why a user stumbles.** The bootstrap tells the user that few treated clusters are a problem and sends them
to `permute()`, which then says nothing about their five treated clusters on the default path. The user
reasonably concludes `permute()` is immune — but randomization inference with a handful of treated clusters
has its own exchangeability and power limits, which the help does discuss (the F049/F051 caveat) without a
run-time trigger that reaches them.

**Pathways.** (a) Add a treated-cluster caution to `permctl`, beside the existing 5902 caution and worded for
RI ("with `n' treated clusters the reference distribution is coarse and rests heavily on exchangeability of
those clusters"). The counting code exists at 6265-6282 and can be lifted. (b) Unify on permute's policy —
drop the bootstrap's count-based cautions and rely on the footer and help. Withdraws a caution the author
shipped on 2026-09-06 after ruling for it; not credible. (c) Cheaper than (a): widen the 5902 caution's gate
— drop the `permtstat` and `cluster()==i` conditions and reword per statistic — and extend the 5817 floor
caution to `permexact` using `wdp_npd`, which is already computed and displayed. That covers both the
exchangeability and the resolution halves with no new counting code. Optionally, extend the bootstrap
caution's last sentence to say what permute's own limitation is, so the referral is not read as a clean bill
of health.

**Defensible?** No — the referral currently points at an engine that is silent on the design it was
recommended for.

**Verification.** Partly confirmed. Corrections applied: the <10-cluster caution is at 5977, not 6260; and
"permute() then says nothing at all" is too strong — the two narrowly-gated cautions above exist, which is
what turned pathway (c) into the cheap option. See CS-10.

#### CS-40 — `safety(off)` means four different things depending on which screen it reaches

**X, the advertised sense.** It downgrades a refusal to a caution and proceeds: `jointtests` at cluster-level
`i` (ado:1522-1523 error under `if "`safetyinput'" != "off"`, caution at 1527); the split-unit `subgroups`
screen (1122-1124 error, 1127 caution); `bootstrap(bayes)` with cell-weighted estimands (6217 error inside
`if "`si'" != "off"`, 6220 caution); `contreatwithin` with `byTime` `contreatcontrols()` (5160 error, 5159
caution). sthlp:977-987 frames this sense, and sthlp:584 and 485 repeat it for the `jointtests` and bayes
screens.

**Y, three other senses.** (i) A caution is suppressed outright so the condition is never mentioned: the
granular-`i` notice inside the gate at 1446 (fires at 1462); the `aweight`+`unconditionalse` caution at
4721-4724, whose own text reads "Set option safety(off) to suppress this caution" and whose behaviour is
unchanged either way — purely verbosity. (ii) A whole error-class check is skipped with nothing printed in
its place: the `customweightmultiplier` non-negativity assert (654); the dropped-dummy multicollinearity scan
and its error inside the gate at 1762; the maxvar preflight whose r(907) error at 1464-1467 sits inside the
same 1446 gate. (iii) Only the **notice** is suppressed while the behaviour change still happens: at
4937-4941, `local espostlength 1` is at 4938, outside every safety test, and only the three explanatory lines
at 4939-4941 are gated — so under `safety(off)` wooldid still silently rewrites the user's `espostlength(0)`
to 1 (the escape is `espostlength(-999)`, which sets `esposttrip` and skips the block).

**Why a user stumbles.** A user types `safety(off)` to get past one specific refusal — say the `jointtests`
cluster screen — and in the same stroke loses the multicollinearity guard, the maxvar guard and the
weight-domain guard, with no list of what was switched off. Sense (iii) is the sharpest: a user reading
"Specify safety(off) to avoid seeing this message" gets exactly that — the message is gone, the rewrite of
their option is not.

**Pathways.** (a) Unify on "downgrade, never silence": every safety-gated screen still runs under
`safety(off)` and prints a caution; nothing disappears. Cost: `safety(off)` runs become noisy, which is
contrary to how some users employ it (the 4724 and 4941 texts sell it as a mute button); the
multicollinearity scan is the expensive one, though less so than it was — since the 2026-09 rework the loop
is guarded by `if strpos("`bvars'", "o.") != 0` at 1766, so on specifications with no omitted term it is a
single `strpos`, not a wordcount-length loop. (b) Unify on "safety(off) means silence **and** skip", and
move the mute cases that are really about verbosity (espostlength, granular-`i`) onto a separate `quiet`-style
control, so `safety()` only ever governs refusals — a new option and a help section; the clean end state.
(c) Keep the four senses and publish a table in the help — screen, default behaviour, behaviour under
`safety(off)` — plus a single line at run time under `safety(off)` naming which screens were disabled for
this call. Cheapest fix to the actual confusion, since the user's problem is not knowing the list.

**Defensible?** No.

**Verification.** Partly confirmed; every ado line number is exact and the `espostlength` case holds
precisely as described. Correction applied: "sthlp 977-987 describes only the first sense" overstates the
documentation gap — that paragraph attributes the multicollinearity, granular-`i` and maxvar screens to
"Default safety", the suppress-the-check sense is documented in roughly nine further passages (sthlp:66 twice,
696, 714-715, 752, 791, 930, 956), and two of the suppressed cautions say so in their own printed text. The
accurate charge is scatter plus the missing run-time list, not silence; sense (iii) is the one genuinely
undocumented item. See CS-25 (the same option under a seven-way taxonomy), CS-34, CS-103.

#### CS-41 — the cluster-at-`i` covariance degeneracy is a refusal, an exemption and a caution-only in three sibling places

**X.** Analytic `jointtests` with clustering at or above `i` is **refused** under default safety
(ado:1522-1523, `error 198` at 1524) with a two-paragraph message and a `safety(off)` escape that downgrades
it to a caution (1527). The measured basis is in the comment at 1485-1488: "measured rejection rates 11-66%
at nominal 5%".

**Y.** The same underlying degeneracy is (i) **exempted** under `bootstrap()` and `permute()` — the gate at
1498 reads `if "`jointtests'" != "" & "`bootstrap'`permute'" == ""`, on the ground that those engines are the
remedy the screen's own text recommends (documented at sthlp:428-429 and 485); (ii) reported as a **caution
only**, with no refusal and no safety gate, for the `boottstat` pairwise contrasts, which studentize by
exactly the degenerate analytic pairwise SE in every replication — 7050-7059, gated only on
`if "`bts'" != "" & "`clus'" != ""`, and its own text names the connection: "the same limitation behind the
refusal of jointtests with analytic clustered SEs"; and (iii) a **caution only** again for `permtstat` at
`cluster()` equal to `i` (5902).

**Why a user stumbles.** A user is refused outright for `jointtests cluster(state)` with `i` = state, learns
that wooldid treats this configuration as unusable, then runs `bootstrap(standard) boottstat cluster(state)
subgroups(g)` and gets studentized contrasts and an across-subgroup equality test built on the very standard
error the refusal was about — flagged by one caution among several. The severity signal inverts between two
calls that differ by one option.

**Pathways.** (a) Unify on refusal: refuse `boottstat` with `cluster()` at or above `i` when subgroup
contrasts or the equality test are requested, with a `safety(off)` downgrade matching 1522's shape. Cost: it
would over-refuse unless the screen actually tests the realized pairwise SEs — the degeneracy is
design-dependent, not universal, and 7050 calls the caution "a disclosure, not a guard". (b) Unify on
caution: downgrade 1522 too, since `safety(off)` already permits it. Cost: contradicts the measured 11-66%
size distortion, and the author has re-ruled in favour of the refusal. (c) Unify the **disclosure** rather
than the policy: give all four sites one shared paragraph in the help under a single heading ("clustering at
or above the cohort level") and have each site point to it, so the user sees one phenomenon with three
treatments and why. Cheaper than it looks — sthlp:501 already carries the boottstat sentence in bold; what is
missing is the cross-link to the `jointtests` paragraph at sthlp:582 and 985-986.

**Defensible?** Site (i), the engine exemption, is defensible on its own terms — the engines genuinely do not
use the analytic covariance for those tests. Sites (ii) and (iii) are the ones that treat the same quantity
differently from the refusal.

**Verification.** Partly confirmed. Corrections applied: (iii) is looser than filed — 5902 fires only when
`cluster()` equals `i` exactly, not at or above, and its stated harm is that `permtstat` is not protected
either, not that the SE is degenerate, so it belongs on the list as a related third treatment with a weaker
claim of identity than (ii), which names the connection itself. See CS-15, CS-99.

#### CS-42 — half the bootstrap's non-results are `e()` flags that survive replay; the permute engine's equivalents are print-once cautions

**X.** The bootstrap posts a flag for every family it could not compute — `e(boot_jt_absent)` (ado:6574,
whose comment reads "a flag, not just a caution, so a REPLAYED run still carries the disclosure"),
`e(boot_adv_absent)` (6689), `e(boot_sgjt_absent)` (6817) — and the footer re-prints them on replay (10115).
Its rank-p resolution is a standing footer line (10112) plus `wooldid_pfloor`'s under-table sentence
(9626-9645).

**Y.** The permutation engine's identical absences — the joint pre-trend test (7770-7771), the
across-subgroup equality test (7858-7859), the adversarial line test (7916-7917) — are accumulated as `wdpc`
locals and printed once by the loop at 8016-8019, with no `e()` flag and no footer sentence. Four
estimation-time cautions are likewise unposted (`wdpmiss`, `wdplow`, `wdpdeg` and the blank `ame` cells,
7660-7662/7664), as are the two pairwise-contrast mapping failures (8009-8015). And `wooldid_pfloor` exits
immediately when `e(boot_completed)` is empty (9635), so the "p = 0 above means p < 2/n" explanation never
prints under `permute()` — including under `permrawp`, the one permute convention where 0 is attainable
(CS-8).

**Replay reaches only footers and `e()`.** The replay path calls `wooldid_display` parts 1/roth/2/perm/boot
(44-50), so nothing outside a footer or an `e()` flag is redisplayed — the analytic path's own cautions (the
Ibragimov-Muller degeneracy notice at 2279, the "SEs not identified for N of M estimates" at 4302) are
equally lost.

**One text detail belongs here.** The bootstrap's stated policy is "never an analytic p under a bootstrap
header" (comments at 6583, 6691, 6819; footer 10115), but the `gteffects` keying-failure caution at 8320
reads "e(gteffects) inference is left analytic and blanked", which asserts both at once — self-contradictory
wording rather than a policy violation, and worth a fix either way.

**Why a user stumbles.** Someone replaying a stored bootstrap run, or handing `e()` to a downstream script,
is told which joint tests are missing. The same person replays a permute run and sees a footer that describes
the tests in general terms while the record of which ones were not computed is gone from the screen and was
never in `e()`.

**Pathways.** (a) Unify on the bootstrap's mechanism: post `e(perm_jt_absent)`, `e(perm_sgjt_absent)`,
`e(perm_adv_absent)` and a `e(perm_degenerate)` list, add the matching footer sentence, and extend
`wooldid_pfloor` to key on `e(perm_completed)` as well (2/(R+1), 1/(R+1) under `permabsp`, c/n under
`permrawp`). Slightly cheaper than it looks: `e(perm_incomplete)` and `e(perm_offcentre)` already establish
the post-plus-footer pattern inside `wooldid_permwire`. (b) Unify downward — drop the `boot_*_absent` flags
and rely on estimation-time cautions everywhere. Loses replay fidelity the author explicitly built and makes
scripted consumption of `e()` less safe. (c) Adopt a rule and write it in the help: "any result wooldid
declines to compute is recorded as an `e()` flag and restated in the footer; cautions are for conditions,
flags are for absences", then audit both engines and the analytic path against it. The audit is the work; the
rule is the payoff, because it also settles the analytic path's replay silence.

**Defensible?** No.

**Verification.** Partly confirmed. Correction applied: permute posts and re-prints **five** footer flags,
not two — `e(perm_gt_blank)` (7573), `e(perm_cfx_blank)` (8084), `e(perm_precaveat)` (7727),
`e(perm_incomplete)` (7644) and `e(perm_offcentre)` (7645), all re-printed at 10061-10067 — and two of them
sit beside condition-type cautions, so the mechanism is not confined to display-blanking cases. The cited
caution range 7660-7665 is too wide (7663 and 7665 are posted). See CS-65, CS-85.

#### CS-43 — an option that does nothing without its master is a hard error nine times and a caution-and-ignore twice

Same taxonomy as **CS-24** (options lens, partly confirmed there) and **CS-78** (display lens, confirmed with
a live run); their corrections and citations apply. **Adds the count and the recorded rationales:** nine
refusal sites (ado:4907, 5416, 5435, 5441, 5448, 6092, 4962-4978, 6030-6045, plus `permprobs`/`permstrata`)
against two cautions (5293, 5357), and the code's own stated rule at 5290-5292 — "inert option: without
makecfxplots the variable is never passed to margins ... caution, do not error" — which is equally true of
`poisexpresults` without `poisson`, recorded at 5433 as an author ruling (SI-9) and refused. Pathways: (a)
unify on error, breaking scripts that carry a superset of options across specifications; (b) unify on
caution, reversing the SI-9, VAL-G19 and NFADV-07 rulings that were made to close silent-accept paths; (c)
write the rule down — "companion options fail loudly except where an established call pattern already carried
them harmlessly" — and list which are which. Defensible per site, incoherent as a set. *Not verified (cap)*
for this framing.

#### CS-44 — `permreps()` has a hard floor and a resolution caution; `reps()` has neither, and its resolution warning arrives only after the run

Same asymmetry as **CS-20** (confirmed there). **Adds the timing contrast:** permute's two protections fire
before a single replication runs (ado:5808, 5817), while the bootstrap's resolution disclosure arrives after
the loop as a footer sentence (10112) plus `wooldid_pfloor`'s under-table line when a p prints as exactly 0
(9632-9645). The two are live author rulings pointing opposite ways — F036 (caution, floor stays 20) and F030
(do nothing) — so pathway (a), a pre-loop bootstrap caution when 2/reps > .05, needs the author to revisit
rather than a build. *Not verified (cap)* for this framing. See CS-9, CS-89, CS-123.

#### CS-45 — the two engines draw the "too many replications failed" line in different places and caution on different fractions

Same fact as **CS-32** (options lens). **Adds the framing:** the two policies are exact mirrors — permute has
a strict refusal (r459 at 7364-7366, below max(20, R/2)) and a lax caution (above 10%, 7372); the bootstrap
has a lax refusal (only at zero, 8610-8612) and a strict caution (any failure, 8614). A user with a fragile
design sees `bootstrap` finish and print inference after most replications failed, then runs the same design
under `permute` and is stopped dead by r(459), with nothing in the table or the help predicting which engine
will refuse. Pathway (c) is mostly a harmonization of one sentence: `e(boot_failed)` already exists (10090)
and the permute footer already prints "`e(perm_completed)' of `e(perm_reps)' replications completed" (10006).
*Not verified (cap).* See CS-90.

#### CS-46 — the three run-wide options are validated inside the bootstrap's private screen

Same defect as **CS-36** (options lens), from the policy side. **Adds the mechanism:** `rseed()` is
validated on every call because `wooldid_bootctl` is invoked unconditionally (ado:4996) and its range check
sits at 6016-6025, deliberately ahead of the stray screen (the comment at 6014-6015 says "rseed() is
validated before the stray screen: it is shared with permute()"). `level()` and `dots()` are checked at 6062
and 6066, **after** the block at 6031-6041 that ends `exit 0` whenever `bootstrap()` is empty — and `level()`
is deliberately omitted from the stray list, whose own error text says "level() applies to every run" (6042).
So the pattern that protects `rseed()` was available and not applied to the two options F054 newly made
meaningful on the analytic path. *Not verified (cap).*

#### CS-47 — two eras of refusal wording sit side by side: the impvalues family states no reason and no remedy, the engine family states mechanism, remedy and alternative

**X.** The older combination refusals are bare declaratives: "Error: Ibragimov and Muller p-values not
compatible with subgroups." (ado:4778), "...not compatible with continuous treatment effects plots." (4793),
"...are not compatible with the jointtests option." (4805), "Error: Adversarial Pre-Trends Testing not
allowed with Ibragimov and Muller p-values." (4729). No mechanism, no remedy, no alternative. One new engine
screen is written the same way: "Error: permute() is not compatible with impvalues." (5733).

**Y.** The 2026 refusals give the reason, the fix and a substitute, often in one long sentence: `boottstat` +
`poisson` (6140), `boottstat` + `bootstrap(bayes)` (6144), `bootdisplay(bc|bca)` + `boottstat` (6169),
`bootdf(clusters)` without `cluster()` (6107), `permute` + `pretrendpower` (5773), and the eight-line
`fastcollapse` refusal at 5117-5126. Where a combination is refused twice from opposite sides the newer side
is the informative one: `fastnose` + `impvalues` (5373) and `pretrendpower` + `impvalues` (9249) explain what
fast mode and the Roth diagnostic need; the impvalues-side twins do not. Across the file, 45 of 217 `dis as
error` lines carry a remedy verb.

**Two smaller inconsistencies in the same theme.** Informational continuation lines print `as error` in three
places ("Note: anticipation(...) shifted every effective onset earlier...", 874, 994, 1035); and the remedy
line under the `unconditionalse` refusal prints `as text` (4717) where the equivalent lines at 711 and 1107
print `as error`.

**Why a user stumbles.** A user who has met wooldid's newer messages expects every refusal to say what to do
next. Hitting "Error: Ibragimov and Muller p-values not compatible with subgroups." they cannot tell whether
the combination is unimplemented, statistically invalid, or merely unverified, nor whether some other option
restores it. The display-class slips read as errors where they are notes and as notes where they are the
actionable half of an error.

**Pathways.** (a) Unify on the new house style: rewrite the impvalues family and 5733 to "Error:
<combination>. <one clause of mechanism>. <what to do>." — about ten message rewrites plus any tests pinning
the text; no behaviour change. (b) Unify downward to terse messages and put all mechanism in the help —
throws away the most-praised part of the recent work; not credible. (c) Adopt a written message contract in
the process checklist (every refusal names the combination, one clause of why, and one remedy; informational
lines print `as text`; remedy lines print in the same class as their error) and apply it opportunistically.
*Not verified (cap).*

#### CS-48 — `customvce()` is refused by the bootstrap for a reason that also holds under `boottstat`, and is silently inert — then quietly honoured — under permute

**X.** `bootstrap()` refuses it: ado:6192, "bootstrap() may not be combined with customvce(); the reported
standard errors would be bootstrap-based, contradicting the custom-VCE request." That rationale is right for
the coefficient-draw families — but under `boottstat` the reported standard error **is** the analytic one
(the footer at 10097 says so), so a custom VCE would play there the same role it plays under `permtstat`.

**Y.** `permute()` has no `customvce` gate at all: `wooldid_permctl` (5722-6005) refuses `unconditionalse`
(5746) and never mentions `customvce`, and the shell passes it into every replication at 128. Under permute's
fast-mode default neither the observed run nor the replications compute standard errors, so `customvce()` is
accepted, carried into the reps, and has no effect whatsoever — no error, no caution, no consequence. Under
`permtstat` it is genuinely honoured and becomes the studentization denominator and the displayed
`se_analytic` column.

**Why a user stumbles.** Three outcomes for one option across two engines — refused, inert-and-silent,
honoured — and the two permute outcomes are selected by whether `permtstat` happens to be present, which the
user thinks of as a p-value convention rather than a VCE switch. A user who specifies `customvce()` under
permute's default and sees no complaint reasonably believes their custom VCE shaped the run.

**Pathways.** (a) Make `permctl` explicit: refuse `customvce()` when the run takes the fast-mode default,
where it can have no effect, and allow it under `permtstat`/`permanalyticse` where it does — one screen,
matching how `permctl` already reasons about `fastnose` at 5799-5801. (b) Relax `bootctl` symmetrically:
allow `customvce()` under `boottstat`, where the analytic SE is the reported quantity, keeping the refusal for
the coefficient-draw families — more work, and it needs verification that the custom VCE survives the
physical-resampling leg, but it removes a refusal whose stated rationale does not apply. (c) Refuse
`customvce()` under both engines unconditionally and say in the help that custom VCEs and replication-based
inference do not mix — withdraws a combination that works today under `permtstat`; simplest to state and
cheapest to hold. *Not verified (cap).*

#### CS-49 — the cell collapse is the only screen family with no `safety(off)` escape, and the only one that skips silently by default

**X.** The collapse's ineligibility list carries two policies at once. With an explicit `fastcollapse` an
ineligible combination is refused with an eight-line r(198) ending "Drop fastcollapse, or drop the
conflicting option. (No safety(off) escape: these are measured wrong answers, not risks.)" (ado:5126).
Without it, the collapse is **silently skipped** and the run proceeds on the full data, with the reason named
only under `verbose` (5055-5057, 2143). The subgroups screen inside the collapse (`wooldid_fcdeg`, 5695-5715)
is a third policy again: under default safety it silently falls back to the full data (message only under
`verbose`), under `safety(off)` it proceeds collapsed with a caution.

**Y.** Every other refusal in the package that rests on a statistical judgment offers `safety(off)`:
`jointtests` at cluster-level `i` (1523), `bootstrap(bayes)` with cell-weighted estimands (6217), the
split-unit subgroups screen (1124), the onset-constancy screen (711), `contreatwithin` with `byTime` (5163).
And nothing else in the package silently declines to do the thing the user asked for.

**Why a user stumbles.** Two learnings clash. A user who has internalized "safety(off) gets me past anything"
is stopped cold by `fastcollapse` with no escape and must read a long message to find that there is none. A
user who has internalized "wooldid tells me when it declines to do something" gets a run that silently
declines the collapse — visible only in `e(wdid_fcskip)` or under `verbose` — and then wonders why an
identically-specified run is ten times slower than a colleague's.

**Pathways.** (a) Keep the refusal, end the silence: print a one-line caution naming the ineligible option
whenever the default collapse is skipped, instead of gating that line on `verbose`. Cost: noise on ordinary
`contreat()` runs — precisely the cost the author avoided when ruling the collapse on by default ("an
ordinary contreat() run must not start erroring", 5056). (b) Keep the silence, end the special-case refusal:
make explicit `fastcollapse` on an ineligible request a caution-and-skip too, so `fastcollapse` reads as
"prefer the collapse" rather than "require it". Destroys the guarantee the refusal provides, and the author
ruled explicitly against it. (c) Make the state a **status** rather than a warning: one line in the run
header — "cell collapse: on" / "off - fullprecisionslow" / "off - ineligible: lattice()" — consistent with
the `e()` flags already posted (F061), and it sidesteps the noise objection because it is not a caution.
Defensible on its merits — a measured wrong answer is genuinely not a risk to be accepted — but it is the one
place where both of the package's general rules fail to hold. *Not verified (cap).*

#### CS-50 — `gteffects` is refused, blanked or fully supported depending on which sibling option accompanies it

**X.** With `contreat()` it is a hard refusal: ado:4786, "Error: gteffects is not available with contreat():
the cohort-by-period interaction coefficients are dose-zero intercepts under a continuous treatment, not
treatment effects. Aggregate e(b_es_ame) or the cfx materials instead."

**Y.** With `permute()` it is accepted, computed, and then blanked — inference columns set to missing
(7565-7573), `e(perm_gt_blank)` posted, a note above the table (9989) and a footer sentence (10064). With
`bootstrap()` it is fully supported: the cells carry their own bootstrap family (9990, 6356-6359), except
that a cell-keying failure produces a caution whose text says the inference is "left analytic and blanked"
(8320) — both of the engines' two stated policies at once, and against the bootstrap's own rule that no
analytic p appears under a bootstrap header (6583, 6691, 6819).

**Why a user stumbles.** One option, three outcomes, chosen by a different option elsewhere on the command
line. A user refused `gteffects` under `contreat()` concludes wooldid refuses what it cannot compute; the
same user under `permute()` gets a full `gteffects` table whose se/p/lb/ub columns are empty and must find the
note to learn why. The 8320 wording then leaves them unable to tell whether the numbers are analytic or
absent.

**Pathways.** (a) Unify on blank-with-disclosure: accept `gteffects` under `contreat()` and blank the cells
with the same note/footer machinery permute uses, since the objection at 4786 is about interpretation rather
than computability. Cost: publishes a table of dose-zero intercepts a reader may mistake for effects — the
reason the refusal exists — so the note would have to be emphatic. (b) Unify on refusal: refuse `gteffects`
under `permute()` too. Removes a table users legitimately want for its point estimates and weight totals,
which permute does compute correctly. (c) Keep all three and give the package one vocabulary for them —
"refused", "reported without inference", "fully supported" — used identically in the messages, footers and
help, and fix 8320 to say which of the three it means (a one-line correction). *Not verified (cap).* See
CS-13, CS-102, CS-128.

#### CS-51 — `makecfxplots` with a non-`att` estimand is a caution that also withholds the requested plots, a fourth policy class the help does not describe

**X.** `makecfxplots` with `att_it`, `att_i` or `att_itime` prints a caution about interpretation and then,
under default safety, adds "As such, wooldid will not produce these plots unless you specify safety(off)."
(ado:4913-4918). The gate that withholds them is at 2655. So the run succeeds, the caution is a caution
rather than an error, and a requested output simply does not appear.

**Y.** Everywhere else the package chooses one of three: refuse the combination (`makecfxplots` without
`contreat()`, 5254), produce the output with a caution attached (the many interpretation cautions, e.g. 4734,
5314), or skip a computation and disclose it in the footer or `e()` (the engines' absent tests).
"Caution and withhold" is a distinct fourth behaviour, and the `safety()` help section (sthlp:977-987)
describes only refuse-versus-caution.

**Why a user stumbles.** The user asked for plots, got no error, got a caution about interpretation — and no
plots. Whether they connect the missing files to the caution depends on reading the second sentence, which
only prints under default safety. A user under `safety(off)` sees the first sentence alone and does get the
plots, so the same caution text means opposite things about what was produced.

**Pathways.** (a) Unify on refusal: make the combination r(198) with the reasoning and the `safety(off)`
escape spelled out, matching 1522's shape — breaks scripts that today request several estimands plus
`makecfxplots` and quietly get the `att` plots only. (b) Unify on produce-with-caution: draw the plots for
every estimand and keep only the interpretation warning — publishes plots the author judges hard to
interpret. (c) Keep the behaviour and make the message unambiguous in both safety modes ("... For that reason
wooldid is producing these plots for att only; safety(off) produces them for the other estimands as well."),
and name the fourth class in the `safety()` help section. *Not verified (cap).* See CS-25 sense (6), CS-98.

#### CS-52 — `subgroups()` screens run over treated units analytically and over every in-sample unit under permute

**X.** The analytic subgroup screens check treated units only, each scoped `if `ifstatement' & `ttre' != .`:
negative values (ado:1071), non-integers (1086), variation within `i`-by-`t` cells (1104), units split across
subgroups (1123).

**Y.** `permctl` re-runs the same four screens over **all** in-sample units and says why: 5926, "requires the
subgroup variable to be non-negative (or missing) on every in-sample row, not only on treated units' rows:
any cluster can be assigned treatment in a replication"; 5931, 5939, 5948 likewise. The `safety(off)` escapes
are preserved (5939, 5948), and 5951/5954 add two permute-specific cautions.

**Why a user stumbles.** The analytic run is clean; adding `permute(shuffle) permreps(999)` — which the user
thinks of as a p-value change — errors at parse on the subgroup variable, on rows of units never treated in
their data. The messages do explain why, so this is a well-handled case; the stumble is expectational rather
than informational, and it is the only place in the package where adding an inference option tightens a
data-integrity screen.

**Pathways.** (a) Unify on the strict domain: run the analytic screens over all in-sample rows too, so the
check is a property of the data rather than of the requested inference — refuses analytic runs that are
correct today (coding a subgroup variable only for treated units is a normal pattern) and would need a
deprecation cycle. (b) Say it once in the `subgroups()` help entry: "under permute() these screens apply to
every unit, because any cluster can be assigned treatment in a replication" — one sentence, and the messages
already carry the reasoning, so it only removes the surprise. (c) Downgrade the permute widening to a caution
where the analytic screen passed — replications would then run on data the engine considers ill-formed,
producing subgroup statistics that mix compositions, which is what the screens exist to prevent. Defensible.
*Not verified (cap).* See CS-106.

#### CS-53 — the no-standard-error state is fatal when the user asks for it and routine when wooldid chooses it

**X.** `fastnose` is refused in combination with five features that need a variance: `unconditionalse`
(ado:5369), `impvalues` (5373), `adversarial` (5377), `jointtests` (5381), `pretrendpower()` (5385). Each
message says fast mode "does not estimate" the covariance those features need.

**Y.** `permute()` without `permtstat` puts the observed run into exactly that state by default (5393-5400),
and `adversarial` and `jointtests` remain fully available, rebuilt as permutation quantities after the
replications (7770-7771, 7916-7917, footer 10056). `pretrendpower()` is the one that stays refused (5773),
and the refusal's rationale — no event-study covariance — sits alongside the permutation-covariance metric
the joint and adversarial tests are built on in the very same engine.

**Why a user stumbles.** "fastnose is not compatible with adversarial" teaches that `adversarial` needs
standard errors. The user then runs `permute()` — which the help says takes the fast path by default — and
`adversarial` works. Nothing in either message distinguishes the two no-SE states, or explains why
`pretrendpower()` alone does not get the substitution the adversarial test gets.

**Pathways.** (a) Reword rather than rebuild: change the `fastnose` messages from "fast mode does not
estimate the covariance" to "fastnose reports point estimates only and provides no substitute; under
permute() these surfaces are rebuilt as permutation quantities instead" — five message rewrites. (b) Extend
the substitution to `pretrendpower()` and lift 5773, as the bootstrap's own ptp refusal was lifted on
2026-09-04 once the ordering constraint was solved. Real research, not a fix: the permutation covariance is a
null-distribution covariance rather than a sampling covariance, so the substitution is not obviously valid.
(c) One help sentence under `fastnose` and `permute()` naming which surfaces are refused, which are
substituted, and which are substituted under one engine but refused under the other. Defensible — the engines
really do supply a substitute that `fastnose` does not. *Not verified (cap).* See CS-3, CS-23, CS-95.

### 2.4 Stored results

#### CS-54 — the default bootstrap family is "percentile" in every macro and "pctl" in every hidden matrix, so the obvious name derivation fails on the default and only on the default

**X.** `bootdisplay(bc)` sets `e(boot_display)` = "bc" and the hidden tables are `e(boot_bc_mainresults)` /
`e(boot_bc_sgresults_#)` / `e(boot_bc_gteffects)` / `e(boot_bc_pwsgs)`; `bootdisplay(wald)` →
`e(boot_wald_mainresults)` (8862, 9006, 7026); `bootdisplay(bca)` → `e(boot_bca_mainresults)` (8865, 9009,
7029). Deriving the table name from the macro works for all three.

**Y.** `bootdisplay(percentile)` — the **default** — sets `e(boot_display)` = "percentile" (9157, from `bdisp`
defaulted at 8791) and `e(boot_gteffects_family)` = "percentile" (6401), but the tables are
`e(boot_pctl_mainresults)` (8863), `e(boot_pctl_sgresults_#)` (9007), `e(boot_pctl_pwsgs)` (7027),
`e(boot_pctl_gteffects)` and `e(boot_pctl_<cfxname>)` (6397, 6428), because of
`local ftag = cond("`fam'" == "percentile", "pctl", "`fam'")` at 6388 and 6414. So
`e(boot_`e(boot_display)'_mainresults)` resolves to `boot_percentile_mainresults`, which does not exist.
Under `boottstat`, `e(boot_display)` still reads "percentile" while the displayed tables are
`e(boot_tpctl_*)` (8916, 9041, 7032). And under `bootdisplay(bca)` the cells fall back to bc (the remap at
6357-6358 is local to `wooldid_bootcellsx`), so one run posts `e(boot_display)` = "bca" and
`e(boot_gteffects_family)` = "bc" — a disagreement the help documents as a cell fallback (sthlp:504) but not
as a macro disagreement.

**The split runs deeper than the hidden tables.** The **same run** posts the long word in the visible-matrix
namespace and the short tag in the hidden-table namespace: `e(boot_ci_percentile)` (9123) and
`e(boot_p_percentile)` (9125) beside `e(boot_pctl_mainresults)`. The wald family splits the other way —
`e(boot_ci_normal)` against `e(boot_wald_mainresults)`, both named in one display line at 10107 (see CS-63).
So no single token derivation works across the visible and hidden halves for either the default or the wald
family.

**Why a user stumbles.** The natural scripted pattern — read `e(boot_display)`, build the twin name, `matlist`
it — works on three of four families and silently fails with r(111) on the one nobody had to type. And a user
reading `e(boot_gteffects_family)` to learn what was displayed gets a different answer from
`e(boot_display)` with no error.

**Pathways.** (a) Unify on the long word: rename the hidden twins `boot_percentile_*` / `boot_tpercentile_*`.
Cost: `boot_percentile_sgresults_#` is 27 characters and safe, but `boot_percentile_` is 16, so the 32-char
guard at 6428 drops the usable cfx-name budget from 22 characters to 16 — more cfx names silently lose their
percentile twin, the failure mode CS-72 describes. Moves 8 files (`grep -rlF boot_pctl_mainresults` over
`test/` returns 6, plus the help and the build copy). (b) Unify on the short tag: make `e(boot_display)` and
`e(boot_gteffects_family)` emit `pctl`/`tpctl`/`bc`/`bca`/`wald`. Cost: the `bootdisplay()` option values
would no longer equal the recorded macro, breaking the other half of the intuition, and battery footer
baselines move. (c) Keep both and post `e(boot_family_tag)`, exactly the twin-name component, documented as
the thing to interpolate. Cheapest; adds a fourth macro to a surface that already has three, and fixes only
the hidden half.

**Defensible?** No.

**Verification.** Confirmed. Two corrections applied: the pathway (a) arithmetic is 22 → 16 usable cfx
characters, not 17; and the asymmetry extends beyond the hidden tables to the visible `ci`/`p` matrices, in
both directions.

#### CS-55 — the per-element inference kernel is a visible matrix under bootstrap, a hidden one under permute, and visible again for permute's own pairwise contrasts

**X.** Under `boottstat` the studentized kernel is posted **visible** as `e(boot_tstud)`, six rows
`n t_obs p_pctl qlo_pctl qhi_pctl live` (ado:9163-9164); the help advertises it (sthlp:501, 1061) and the
caution text points users at row `n` (9169-9170). Its contrast twin `e(boot_tstud_pwsgs)` is also visible
(7034-7036) — unlike its siblings `boot_tpctl_pwsgs`, `boot_wald_pwsgs`, `boot_pctl_pwsgs`,
`boot_bc_pwsgs` at 7025-7032, which are all `ereturn hidden`.

**Y.** Under `permute()` the equivalent object — richer, eighteen rows
`n s_obs p qlo qhi sd lb ub se_shift offcentre bandlo bandhi live p_lo p_hi se_p cp_lo cp_hi` — is posted
**hidden** as `e(perm_kernel)` (7626-7627), so `ereturn list` does not show it. But the identical eighteen-row
object for the pairwise contrasts, `e(perm_pwsgs_kernel)`, is posted **visible** (8004-8006). The row
vocabularies also drift: `t_obs` vs `s_obs`, `p_pctl` vs `p`, `qlo_pctl` vs `qlo`; only `n` and `live` are
shared. The footer at 10057 names `e(perm_pwsgs_kernel)` to users and no display line anywhere names
`e(perm_kernel)`. The default (non-`boottstat`) bootstrap posts no per-element kernel at all.

**Why a user stumbles.** A user who learned to inspect `e(boot_tstud)` after a studentized bootstrap runs
`ereturn list` after `permute()` and concludes no kernel exists, then reconstructs p-values by hand from
`e(permdist_b)`. If they have subgroups they **do** see `e(perm_pwsgs_kernel)` and reasonably infer the
main-coefficient version was omitted rather than hidden.

**Pathways.** (a) Unify on visible: drop `hidden` from 7627. `ereturn list` after permute grows by one wide
matrix; the churn is lower than it looks, since unhiding does not break `matrix X = e(perm_kernel)` readers,
and the 7 files that mention the name are mostly a build copy, MC scripts and review records rather than
`ereturn list` baselines. (b) Unify on hidden: hide `e(perm_pwsgs_kernel)` and `e(boot_tstud)` too — the
boottstat caution (9170) and the permute footer (10057) both point users at these by name and would have to
say "hidden"; loses a documented debugging aid. (c) Keep both, state the rule in the help, and harmonize the
row stripes onto one vocabulary (`n, s_obs, p, qlo, qhi, live`) so a script can read either — row-name
changes break any test pinning `e(boot_tstud)` rownames.

**Defensible?** No.

**Verification.** Confirmed, line by line. One correction that softens the stumble: the help *does* document
`e(perm_kernel)` at sthlp:1068, explicitly as hidden, so it is discoverable — just not enumerable by
`ereturn list`, which is how most users inventory a run.

#### CS-56 — the Monte-Carlo precision record for a p-value: four visible named matrices under permute, one hidden matrix under bootstrap, shipped in the same commit

**X.** `permute()` posts `e(perm_p_lo)`, `e(perm_p_hi)` (one-sided shares), `e(perm_p_se)` (MC standard error
of the share) and `e(perm_p_ci)` (Clopper-Pearson, rows lo/hi) as four **visible** stripe-aligned matrices
(ado:7599-7614), the footer names all four (10036), and the help documents them three times (sthlp:477, 481,
1066).

**Y.** `bootstrap()` computes the same five numbers and posts them as one **hidden** matrix,
`e(boot_ptails)`, rows `p_lo p_hi se_p cp_lo cp_hi` (9126-9133, whose own comment calls it "the bootstrap
analogue of the permutation tail record"). Nothing printed mentions it — the bootstrap footer's only
precision sentence is the 2/n resolution line — and `boot_ptails` appears **nowhere** in `wooldid.sthlp`.
Both shipped as F064 in `dc3f345`, where the ship note calls the bootstrap side "hidden".

**A wider asymmetry than visible-versus-hidden.** `boot_ptails` takes its colnames from
`e(boot_p_percentile)` and its rows are percentile shares about 0. Under `bootdisplay(bc)` or `(bca)` the
**displayed** p is BC or BCa, and there is no MC-precision record for it at all, hidden or otherwise. So
permute's record tracks the convention actually shown — the footer even says "doubled under the equal-tailed
rule" — while the bootstrap's hidden record tracks one family regardless of what was displayed.

**Why a user stumbles.** A user who reported `e(perm_p_se)` in a robustness appendix switches engine and
finds no `e(boot_p_se)`; `ereturn list` shows nothing and the help says nothing, so they conclude the
bootstrap p has no Monte-Carlo standard error and either compute it by hand or drop the column.

**Pathways.** (a) Unify on the permute shape: split `e(boot_ptails)` into visible
`e(boot_p_lo)`/`e(boot_p_hi)`/`e(boot_p_se)`/`e(boot_p_ci)` and add a footer sentence. Costs more than four
renames, because a faithful visible `e(boot_p_se)` would have to follow `bootdisplay()` rather than always
reporting the percentile family. Two source files reference `boot_ptails` today. (b) Unify on the bootstrap
shape: collapse the four permute matrices into a hidden `e(perm_ptails)`. Costs more than filed — 9 files
reference `perm_p_lo`/`perm_p_se`, and the permute footer and help name the four matrices. (c) Keep both and
document the mapping in one help sentence: "the bootstrap analogue of `e(perm_p_*)` is the hidden
`e(boot_ptails)`, rows p_lo p_hi se_p cp_lo cp_hi". Cheapest; leaves the asymmetry.

**Defensible?** No.

**Verification.** Confirmed on code, footer and help. Corrections applied: the pathway file counts (2 rather
than 3 for `boot_ptails`; 9 rather than 6 for `perm_p_*`), and the display-family scope point above, which
the finding did not raise. The verifier also noted permute posts *both* a hidden kernel and the four visible
matrices, where the bootstrap has only the hidden one. Anchor for CS-12, CS-83.

#### CS-57 — a two-row interval is named `ll`/`ul`, `lo`/`hi`, `lb`/`ub` or `lb_95ci`/`ub_95ci` depending on which object you opened, and `e(perm_lb_att)` is not the lower bound of the reported interval

*Filed high; verifier: medium-high — the printed output labels the quantity correctly, so the risk is confined
to users reading scalars programmatically.*

**X.** `e(perm_ci)` (the reported shift interval) and every `e(boot_ci_*)` use rows `ll ul` (ado:7584;
9116-9118 and 9138); the analytic result tables use **columns** `lb_95ci`/`ub_95ci`, renamed to the level in
force (3849-3850).

**Y.** In the same permute run, `e(perm_nullband)` (7597) and `e(perm_p_ci)` (7613) use rows `lo hi`, while
`e(perm_kernel)` names the very rows that become `e(perm_ci)` `lb`/`ub` and the very rows that become
`e(perm_nullband)` `bandlo`/`bandhi` (7628). And the legacy pooled scalars `e(perm_lb_<est>)` /
`e(perm_ub_<est>)` are read from kernel rows **11-12** — the null band — not rows 7-8:
`ereturn scalar perm_lb_`wdpe' = `wdpk'[11, `wdpc']` (205-206). So `lb` in a scalar name and `lb` in the
kernel stripe denote two different quantities in one run, and neither is the `ll` of `e(perm_ci)`.

**Why a user stumbles.** `matrix li e(perm_ci)` then `e(perm_ci)["lo",...]` fails. And `di e(perm_lb_att)`
returns a number that looks like a CI bound and is the 2.5th percentile of the null draws. The printed output
is not to blame: the footer heads the block "Pooled estimand summaries (p, `level'% null band and SD of the
permutation distribution)" and prints "null band [lb, ub]" (10068-10071), and the help calls them the null
band too (sthlp:415, 1074). The misleading surface is the **scalar name**, read out of a saved-results
do-file by a user who has learned that `ll`/`ul`/`lb`/`ub` mean interval bounds everywhere else.

**Pathways.** (a) Unify on `ll`/`ul` wherever a two-row interval is posted (`perm_nullband`, `perm_p_ci`,
kernel rows 7-8 and 11-12 relabelled `ll`/`ul` and `band_ll`/`band_ul`). Cost: 8 files reference
`perm_nullband`, 4 reference `perm_p_ci`, and the kernel stripe is pinned in the permute cert (Q19-Q23). (b)
Unify on `lo`/`hi` and change `e(perm_ci)`/`e(boot_ci_*)` — much larger blast radius (13 files touch
`boot_ci_percentile`) and it diverges from `r(table)`'s `ll`/`ul`, which the cfx matrices already inherit.
(c) Rename the legacy scalars `e(perm_bandlo_<est>)`/`e(perm_bandhi_<est>)` with the old names kept as
aliases for one release. This is the half that actually misleads, and it costs the help's legacy paragraph
plus any user do-file reading `perm_lb_att`.

**Defensible?** No for the scalar names; the four stripes themselves are partly historical and partly
inherited from `r(table)`.

**Verification.** Confirmed; every code citation verbatim. Correction applied: the finding said the value
"prints beside the estimate" as if it were a CI bound — the footer and help both label it as the null band, so
the risk is programmatic, not on-screen. Priority softened accordingly.

#### CS-58 — `e(permdist)` and `e(permdist_b)` are different matrices on different scales, and the bootstrap has one draw matrix and hides its comparator

**X.** `bootstrap()` posts one draw matrix, `e(bootdist)` (colnames `obn`, stripe-aligned to `e(b)`, capped at
10,000 completed replications, ado:9220/9227), plus `e(bootdist_se)` (9222) and `e(bootdist_pwse)` (7041)
under `boottstat`. Its rep-0 comparator, `e(boot_comparator)`, is **hidden** (9211), and no per-cell draws
are stored at all.

**Y.** `permute()` posts **two**: `e(permdist_b)` (7640), colnames `obn`, always on the coefficient scale;
and the retained legacy `e(permdist)` (224), one column per requested **estimand** and on the `e(perm_stat)`
scale — i.e. the t-ratio scale under `permtstat`, since the source cell is filled
`_b[main_`ee']/_se[main_`ee']` at 7320 against `_b[main_`ee']` at 7312, and sthlp:1075 says so. The
comparator doubles too: visible `e(perm_comparator_b)` (7617) and visible `e(perm_comparator)` (227).
Permute additionally stores per-cell draws as hidden `e(permdist_<cfx>)`/`e(permcmp_<cfx>)` and their `se`
twins (8067-8079). The 10,000-cap applies to all three draw matrices alike, so that is not an asymmetry.

**Why a user stumbles.** `svmat e(bootdist)` after a bootstrap gives K coefficient columns; `svmat
e(permdist)` after a permute gives one column per estimand, and under `permtstat` on a completely different
scale — a plausible-looking but wrong "distribution of the ATT". The `_b` suffix reads like "the b-part of
permdist" rather than "the real one".

**Mitigations worth recording.** `e(permdist)` carries estimand colnames (7389), so `matrix list` labels the
columns and the trap is the **scale**, not the shape; `svmat` is where the labels are lost. The scale is
surfaced at display time for the objects `permdist` backs — the footer appends ", on the t-ratio scale" to
the pooled-estimand summary header when `e(perm_stat)` is `t` (10068). And both engines expose the rep-0
comparator in their save files (`bootsave` sthlp:506, `permsave` sthlp:462, each writing a `_rep = 0` row),
so "hidden" is an `e()`-surface fact only.

**Pathways.** (a) Deprecate `e(permdist)`/`e(perm_comparator)` in favour of the `_b` versions (aliases for one
release, or drop). Cost: the live consumers are `test/suite/newfeat/ri_crosscheck.do` (20 hits),
`test/test_permute.do`, `test/test_permute2.do`, the help and the standing-issues record; the ~11 and ~29
raw file counts are mostly stale logs and documentation. (b) Give the bootstrap the same legacy twins — adds
objects nobody asked for. (c) Rename for clarity: `e(permdist_est)`/`e(perm_comparator_est)` for the
estimand-keyed pair, freeing `e(permdist)` to alias `e(permdist_b)` so the two engines' names match. One help
paragraph plus the `permsave` cross-reference; breaks any script reading the current `e(permdist)`.

**Defensible?** No — the name is the trap.

**Verification.** Confirmed. Corrections applied: the pathway file counts above; the colname mitigation; the
footer's t-scale disclosure; the save-file comparator; and the fact that the 10,000-cap is symmetric and
should not be read as another asymmetry. See CS-33.

#### CS-59 — `e(V)` has three regimes across the engines and only the bootstrap aligns `e(df_r)` with what it printed

Companion to **CS-4** (engines lens, confirmed), which carries the four-provenance framing and the
`df_r_analytic` branch gap. This item is the stored-results view and adds what `e()` does **not** carry.

**X.** Analytic: `e(V)` posted, diagonal-only. Bootstrap: `e(V)` reposted with the bootstrap covariance
(ado:9047-9067), the analytic one preserved as hidden `e(V_analytic)`, a visible `e(V_boot)` also posted,
`e(boot_v_posted)`/`e(boot_N_V)` recording whether the repost succeeded and on how many replications,
`e(V_es_*)` replaced by bootstrap subblocks with `e(V_es_analytic_*)` hidden and `e(boot_ves_failed)` on
failure (9077-9098), and `e(df_r)` set to the reference distribution actually reported (missing under
`bootdf(normal)`) with the analytic value kept as hidden `e(df_r_analytic)` (9174-9187).

**Y.** Permute in the default fast mode: **no** `e(V)` and no `e(V_es_*)` — and, the verifier adds, no
`e(b_es_*)` either, because the post at 4348 sits inside the same `f2ok` gate that all-missing variances
fail. The footer says so in prose (10029-10030) and the help says so in the `permute()` entry (sthlp:417), but
there is no flag: a grep finds no `e(perm_v_posted)`, no `e(perm_v_source)`, no `e(V_perm)` and no `perm_df*`
anywhere in the ado, and the Permutation-inference stored-results block (sthlp:1063) never mentions `e(V)` or
`e(df_r)` at all. Under `permanalyticse` or `permtstat`, `e(V)` exists but is the **analytic** covariance
under a permutation header, and `e(df_r)` is the analytic residual df beside permutation p-values, with no
alignment and no `df_r_analytic` twin.

**Why a user stumbles.** `esttab` and `estout` read `e(b)` and `e(V)`; after analytic and bootstrap runs they
work, after a default permute run there is nothing to read and the user gets a coefficient-only table or an
error, with no `e()` flag to test for. A user who ran `test` after a bootstrap and got the documented
`bootdisplay(wald)` answer runs it after `permute() permanalyticse` and gets an analytic Wald test that
contradicts the printed `p_perm`, with only prose warning them.

**Pathways.** (a) Always post something: give permute a zero-off-diagonal `e(V)` built from
`e(perm_se_shift)` or `e(perm_sd)` so `esttab` works, plus `e(perm_v_source)`. Invites exactly the misuse the
fast-mode design avoids — a shift-implied SE is not a sampling SD, and the author already rejected the
analogous inverted-interval SE for the bootstrap (CS-71). (b) Make the bootstrap not repost `e(V)` either,
leaving only `e(V_boot)` — breaks the documented `esttab`/`test` story the bootstrap help spends a paragraph
on. (c) Post machine-readable flags on the permute side — `e(perm_v_posted)` 0/1 and `e(perm_v_source)` in
{none, analytic} — mirroring `e(boot_v_posted)`, and align or blank `e(df_r)` under permute the way the
bootstrap does (or post `e(perm_df_note)`); optionally flag the withheld `e(b_es_*)` too. Modest, and it makes
the regime testable instead of prose-only.

**Defensible?** Yes on the three regimes; not on the missing flag.

**Verification.** Partly confirmed (code and help; no run — the box was at the cap, but the no-`e(V)` chain is
traced end to end: missing diagonal → `ereturn post` failure at 4275 → b-only fallback at 4306, with the
explanatory caution deliberately suppressed under `fastnose` at 4308). Corrections applied: `e(b_es_*)` is
withheld alongside `e(V_es_*)`; and the absence is documented at sthlp:417, not footer-prose only.

#### CS-60 — the same joint pre-trend p-value is an `e()` scalar without `subgroups()` and an `e()` local macro with it

**X.** `e(joint_horiz)`, `e(joint_horiz_i)`, `e(joint_Choriz)`, `e(joint_Choriz_i)` are posted with `ereturn
scalar` (ado:4444-4447), so `di e(joint_horiz) < 0.05` and `scalar p = e(joint_horiz)` work. The adversarial
p-values are scalars at **both** levels — `e(adversarial_pval_att)` (4459) and
`e(adversarial_pval_att_s#)` (4483) — and `e(joint_att_allsgssame)` is a scalar (4494-4501).

**Y.** The subgroup versions of the joint horizon tests — `e(joint_horiz_s#)`, `e(joint_horiz_i_s#)`,
`e(joint_Choriz_s#)`, `e(joint_Choriz_i_s#)` — are posted with `ereturn local` (4470-4473). The help codifies
it: sthlp:1111-1112 says the overall tests are "presented in the scalars", sthlp:1115 says the subgroup ones
are "presented by subgroup, when used, in the local macros", while sthlp:1119 and 1125 say "in the scalars"
for the adversarial and `allsgssame` families. So one statistic changes `e()` type at the subgroup boundary,
and only this one family does.

**Both engines replicate the split deliberately.** `wooldid_bootjnt` at 6541-6543 and `wooldid_permjnt` at
7752-7754 each do
`local en = cond("`sgtag'" == "", "`ename'", "`ename'_s`sgtag'")` then
`if "`sgtag'" == "" ereturn scalar `en' = `jp'` / `else ereturn local `en' = `jp'`, with the comment at 6611
naming it "(ereturn LOCALS - the analytic macro convention)".

**Why a user stumbles.** `if e(joint_horiz_s1) < .05 { }` compiles and silently compares a string to a
number; a loop building a table of joint-test p-values across subgroups gets string cells or a type error,
while the identical loop over adversarial p-values works.

**Pathways.** (a) Unify on scalar. Note the real scope: **three** posting sites, not one — 4470-4473,
~6541-6543 and ~7752-7754 (plus the `boottstat` twin near 6591-6620) — and the display code at 9928-9942
dereferences the subgroup values as macros (`%-6.4f `e(joint_horiz_s`sgcase')'`), so it must be edited too.
Existence tests of the form `"`e(joint_horiz_s1)'" != ""` mostly survive, since macro-quoting an `e()` scalar
yields its string form; code doing exact string comparison breaks. Sizing:
`grep -rlF 'e(joint_horiz_s'` hits `test/test_refmargin.do`, `test/test_fastcollapse.do`, the build copy and
several `review_process/` records. (b) Unify on macro — wrong direction; scalars are the Stata norm for
p-values. (c) Keep and document loudly in the help synopt with a worked `real("`e(joint_horiz_s1)'")`
example. Cheapest, leaves a silent type trap.

**Defensible?** No.

**Verification.** Confirmed (type read from the posting verbs and the help, not a live `ereturn list`; the
verbs are unambiguous). Corrections applied: the pathway's claim that the `allsgssame` existence idioms at
6733/7862 need rewriting is wrong — `allsgssame` is already a scalar and those idioms are unaffected; the
site that *does* need editing is the display block at 9928-9942, which the finding omitted. See CS-101.

#### CS-61 — `e(perm_pconv)` takes five values, the stored-results block documents two, and `e(perm_absp)` reads 0 for three of them

*Filed high; verifier: medium — the five values are enumerated in the option entry, so the defect is an
intra-help lag plus an undescribed scalar.*

**X.** The bootstrap's convention macro is complete in the code: `e(boot_display)` is exactly the
`bootdisplay()` value and `e(boot_jt_family)` is `studentized`|`fixedmetric`. (Its help does not enumerate the
values either — sthlp:1061 says only "records which family was displayed" — so the bootstrap side is complete
in the code, not in the help.)

**Y.** `e(perm_pconv)` is set to one of five strings by a `cond()` chain over the five-valued `wdpabs` code
(ado:7648; `wdpabs` set at 7421, 7685, 8037: rawp = 2, side left = 3, right = 4, absp = 1, else 0). The
companion scalar is `ereturn scalar perm_absp = (`wdpabs' == 1)` (7649), so it reads 0 under `permrawp`,
`permside(left)` and `permside(right)` — indistinguishable from the equal-tailed default. The help is split
against itself: the `permabsp` option entry lists all five ("equaltailed or absolute; raw, left or right
under the two options below", sthlp:475), while the stored-results block says "`e(perm_pconv)` (equaltailed
or absolute)" (sthlp:1073), lists `e(perm_absp)` with no description at all, and describes `e(perm_p)` as
"the equal-tailed permutation p, or the absolute-value p under permabsp" (sthlp:1065) — two of the five
conventions its content can carry.

**Why a user stumbles.** A user auditing saved permute models writes `if e(perm_absp) == 1` to find the
absolute-value runs, and a `permrawp` run — which the help's own account calls an absolute-value share — is
classified with the equal-tailed ones. Reading only the stored-results block, they see a binary macro beside
a binary-looking scalar and conclude the run was equal-tailed when it was `permside(right)`.

**Pathways.** (a) Fix the help: list all five values in the `e(perm_pconv)` and `e(perm_p)` entries, and
describe `e(perm_absp)` as "1 only under permabsp; the convention in force is e(perm_pconv)". Two help lines
plus a clause. Widening `e(perm_absp)` to a coded scalar would change its meaning for existing scripts, so
deprecating in place is safer. (b) Post a second boolean per convention (`e(perm_rawp)`, `e(perm_side)`) —
symmetric, but it multiplies redundant state for a problem the existing macro already solves. (c) Fix only the
help and leave `e(perm_absp)` undescribed — minimum viable; the trap remains.

**Defensible?** The scalar's name is defensible on its own terms (it means literally "permabsp was given");
the two-valued stored-results enumeration is not.

**Verification.** Partly confirmed. Correction applied: "reading the help, they do not learn that
e(perm_pconv) can say 'raw' at all" is refuted — sthlp:475 states all five. The defect is that the option
entry is current and the stored-results section is not, which is the same lag as CS-119. Priority dropped to
medium.

#### CS-62 — `e(level)` follows `set level` and `e(boot_level)`/`e(perm_level)` default to 95, so the two disagree

The `e()` half of **CS-2**, verified independently by a second lens; both are confirmed and their citations
agree. This item adds three things.

**First, the shell passes the raw option at eight call sites**, not six: ado:195, 229, 236, 394, 396, 438 and
also the two controller calls at 4982 (`wooldid_permctl`) and 4996 (`wooldid_bootctl`). The correct
resolution already exists twice in the tree as a model — 418 uses
`level(`=cond(`level' == 0, c(level), `level')')`, and 4595-4596 pass the already-resolved `level(`wdlvl')`.

**Second, there is a third `level` default convention.** `wooldid_ptp` declares `level(real 95)` in its own
syntax line (9279), so the same option has three defaults in one file: 0 in the shell, 95 in the engines, 95
in `wooldid_ptp` — which is why 418 had to hand-resolve. A single `local level = cond(`level' == 0,
c(level), `level')` in the shell before any engine call removes all three.

**Third, the syntax line is `Level(real 0)`, not Stata's `cilevel` keyword** (506, and the byte-identical
twin at 4670), and `grep -n "c_local level"` returns nothing, so `wooldid_validate` does not resolve it
either — an omitted `level()` really is 0 all the way down. The engine calls at 195/229/236/394/396/438 sit
in the shell **before** the worker's `set level` at 525-528 ever runs, so `level` there is unambiguously the
raw 0. And the subprogram default appears at **six** sites (6369, 6454, 7417, 7681, 8033, 8639), not two.

**Pathways.** As CS-2, with (a) — resolve once in the shell — now the clear winner, since it also removes the
`wooldid_ptp` third default. Pathway (c), "document that `level()` must be given explicitly under the
engines", is off the table: sthlp:489 already promises the opposite, in a sentence shipped 2026-09-06.

**Defensible?** No.

**Verification.** Confirmed (code only; the reproducer at
`/Users/hegland/.claude/jobs/be9644d0/tmp/consist/skeptic_level/lvl.do` was written and never run, the box
being at the cap — the path has no branch that could rescue it). See CS-2, CS-36, CS-46, CS-116, CS-117.

#### CS-63 — one bootstrap family answers to "wald", "normal" and "z/t", "normal" separately names the reference distribution, and it is the only family with no stripe-aligned p matrix

**X.** The percentile, bias-corrected and BCa families each get a matching pair of stripe-aligned matrices:
`e(boot_ci_percentile)`+`e(boot_p_percentile)`, `e(boot_ci_bc)`+`e(boot_p_bc)`,
`e(boot_ci_bca)`+`e(boot_p_bca)`.

**Y.** The normal/t family is selected by `bootdisplay(wald)`, stored as `e(boot_wald_mainresults)` (ado:8862),
and its interval matrix is `e(boot_ci_normal)` (9122) — three tokens for one family. There is no
`e(boot_p_wald)` or `e(boot_p_normal)` at all, so its p-values exist only inside column 4 of the hidden
table. Meanwhile `bootdf(normal)` uses the same word for a different concept entirely (z versus t critical
values), and `e(boot_stat)` reports z/t for that setting; sthlp:503 and 504 use "normal" for the two
different things in adjacent paragraphs.

**Why a user stumbles.** A user who found `e(boot_ci_bc)` and `e(boot_p_bc)` reaches for `e(boot_ci_wald)`
(absent) and then `e(boot_p_normal)` (absent). Someone reading `bootdf(normal)` in a colleague's do-file may
believe the Wald display was requested when only the reference distribution was set.

**Pathways.** (a) Unify on "wald": rename `e(boot_ci_normal)` → `e(boot_ci_wald)` and add `e(boot_p_wald)`,
which is already computed as column 4 of `e(boot_wald_mainresults)`. Cost: `boot_ci_normal` appears in the
footer string at 10108 and in the help's stored-results synopt; 9 files reference it, plus the battery footer
baselines. (b) Unify on "normal": rename `bootdisplay(wald)` → `bootdisplay(normal)` and
`boot_wald_mainresults` → `boot_normal_mainresults` — collides head-on with `bootdf(normal)` and with the
F053 `t_wald`/`z_wald` column labels shipped 2026-09-06; worse. (c) Keep the names, add `e(boot_p_wald)` for
symmetry plus one help sentence saying wald = normal = the `ci_normal` matrix, and rename the `bootdf()`
values to `bootdf(z)`/`bootdf(clusters)`/`#` with an alias for back-compatibility. *Not verified (cap).* See
CS-54.

#### CS-64 — permutation gives you `e(perm_p_att)` as a scalar; the bootstrap gives you nothing per-estimand

**X.** After `permute()`, the pooled estimands have direct scalars: `e(perm_p_att)`, `e(perm_lb_att)`,
`e(perm_ub_att)`, `e(perm_sd_att)`, and under `permtci` `e(perm_tlb_att)`/`e(perm_tub_att)` — one set per
requested estimand, posted by the `foreach` over `wdp_est` at ado:203-219, so `di e(perm_p_att)` just works.

**Y.** After `bootstrap()` there is no `e(boot_p_att)`, no `e(boot_lb_att)`, nothing keyed by estimand. The
same number requires `matrix P = e(boot_p_percentile)` then `di P[1, colnumb(P, "main_att")]`, and the user
has to know that the display family and the matrix family are chosen separately. The bootstrap's
per-coefficient surfaces are all 1xK matrices (9122-9143, 9207).

**Why a user stumbles.** A user who wrote `local p = e(perm_p_att)` in a table-building do-file and then
switched the engine to `bootstrap()` gets an empty local, silently, because Stata returns missing for an
unset `e()` scalar rather than erroring.

**Pathways.** (a) Unify upward: post `e(boot_p_<est>)`, `e(boot_lb_<est>)`, `e(boot_ub_<est>)` from the
**displayed** family plus a macro saying which family they came from — three new scalars per estimand and a
decision about which family they carry when `bootdisplay()` is not the default, which is exactly the
ambiguity that kept them out. (b) Unify downward: withdraw the permute legacy scalars — they predate the
matrix surface, the help calls them retained, and user scripts break. (c) Add a one-line help note under the
bootstrap synopt showing the `colnumb()` idiom, so the asymmetry is at least discoverable. *Not verified
(cap).* Note that CS-57 is about these same legacy scalars from the other direction.

#### CS-65 — only the bootstrap records "this test is absent" in `e()`, and the two engines use opposite flag polarity

Companion to **CS-42** (policy lens, partly confirmed), whose correction — permute posts five footer flags,
not two — applies here. **Adds two facts this lens found.** *Polarity:* permute posts the **positive**
`e(perm_adv) = 1` when the adversarial test ran (ado:7937); the bootstrap posts the **negative**
`e(boot_adv_absent) = 1` when it did not (6689). A script cannot test the two the same way. *Observed-statistic
twins:* `e(perm_jointWobs)`/`e(perm_sgWobs)`/`e(perm_advWobs)` always accompany their replication matrices
(7776-7777, 7866-7867, 7934-7935), while `e(boot_jointWobs)`/`e(boot_sgWobs)` exist only in the fixed-metric
branch (6568-6569) — the studentized branch posts `e(boot_jointW)` alone (6638) — and there is no
`e(boot_advWobs)` at all (6698). And a grep confirms no `e(perm_jt_absent)`, `e(perm_sgjt_absent)`,
`e(perm_adv_absent)` or `e(perm_jt_family)` anywhere in the file. Pathways as CS-42, plus: choose one
polarity and post the missing `Wobs` twins. *Not verified (cap)* for this framing.

#### CS-66 — three column vocabularies for the same six numbers inside one analytic run

**X.** `e(wdidmainresults)`, `e(wdidsgresults_#)` and `e(pwsgs)` carry columns
`estimate se t p lb_95ci ub_95ci relyr` (ado:3849-3850, 3976-3977, 4019-4020), with `z` for `t` under
poisson. `e(gteffects)` carries `i onset caltime relyr estimate se p lb_<lvl>ci ub_<lvl>ci inwindow nobs
w_att w_att_it w_att_i w_att_itime` (10848) — the same words for estimate/se/p/bounds, but no statistic
column at all.

**Y.** Every `e(cfx*_*)` matrix carries `b se t pvalue ll ul` — margins' `r(table)` row names, transposed
(2720, 2757; the column map is stated in the comment at 6404) — so the estimate is `b` not `estimate`, the
p-value is `pvalue` not `p`, and the bounds are `ll`/`ul` not `lb_95ci`/`ub_95ci`, and they do not move with
`level()`.

**Why a user stumbles.** A user who wrote `svmat e(wdidmainresults), names(col)` then `keep estimate lb_95ci
ub_95ci` repeats it on a cfx matrix and gets "variable estimate not found". The bound names in the main
tables change with `level()` while the cfx bound names never do, so a level-aware script has to handle three
schemes.

**Pathways.** (a) Relabel the cfx matrices to the analytic vocabulary at the point they are transposed. The
cfx plot code reads columns positionally in most places but by name in the plotting rename block
(10556-10680), and both engines' cell rewiring addresses columns 1/2/4/5/6 positionally, so the rename is
mostly cosmetic — but the cfx colnames are pinned in the cfx cert and in `saveplots` datasets users may
already have. (b) Unify on `r(table)`'s `b se t pvalue ll ul` everywhere — much larger, since
`lb_95ci`/`ub_95ci` appear in the `saveplots` dataset contract, the plot code and dozens of battery
baselines. (c) Keep all three and add a "column dictionary" to the help's Working With Output section
mapping the three vocabularies — cheapest, and probably right, since (b) would break saved plot data.
Defensible in part: the `r(table)` names are a deliberate "this is margins output" signal. *Not verified
(cap).* See CS-79, CS-84.

#### CS-67 — the engines relabel the main-table columns to disclose what they hold, but only some columns, differently per engine, and never in the cell tables

Overlaps **CS-26** and **CS-74** (both confirmed), whose corrections apply. **Adds the cell tables, which are
this item's real content:** under **both** engines `e(gteffects)` and every `e(cfx*)` keep their analytic
column names verbatim — `matrix colnames `gt`fam'' = `: colnames e(gteffects)'` at ado:6398 and
`= `: colnames e(`nm')'` at 6423-6424, with permute's cell path at 8055-8056 doing the same — so a
bootstrapped `gteffects` table still says `se` and `p` for a bootstrap SD and a rank. Pathways: (a) extend the
marking to `p`/`se` on the main tables **and** to the cell tables under both engines — every bootstrap battery
baseline and boot cert re-pins; (b) unify on the bootstrap discipline (rename only the statistic) — loses
permute's `se_analytic` disclosure, which the footer at 10029-10030 leans on; (c) at minimum extend the
rename to the cell tables so a matrix never carries analytic names over engine numbers, and say in the help
which columns are renamed per engine. The cell tables are the sharpest gap, because they are what plots and
appendices are built from. *Not verified (cap)* for this framing. See CS-13.

#### CS-68 — the cluster count is `e(Nclus)`, not Stata's `e(N_clust)`, and wooldid does the renaming itself

**X.** wooldid reads the underlying regression's cluster count as `e(N_clust)` — Stata's universal
convention, which `esttab`/`estout`/`outreg2` all look for by that name — at ado:1887.

**Y.** It re-posts it as `e(Nclus)` (4437, `foreach scalar in N Nclus ll ll_0 `r2list'`), and that is what the
help documents and the display header prints (9699, 9722-9723). Neither engine posts the number of
**resampling** units either: `bootdf(clusters)` uses G-1 and the F050 caution counts treated clusters, but G
reaches `e()` only implicitly through `e(boot_df)`; permute posts `e(perm_npossible)` but not the cluster
count.

**Why a user stumbles.** `esttab, stats(N N_clust)` prints a blank cluster row after wooldid and a filled one
after `reghdfe`/`areg` in the same table; the user assumes wooldid did not cluster.

**Pathways.** (a) Post `e(N_clust)` as well, keeping `e(Nclus)` as an alias — near-zero risk, it is an
addition; ~20 files reference `e(Nclus)` and none break. (b) Rename outright and drop `e(Nclus)` — breaks
those files and any user script; not worth it when (a) is free. (c) Do nothing and note the difference beside
the `esttab` compatibility paragraph. Separately, add `e(boot_nclus)`/`e(perm_nclus)` if the author wants G
recoverable — today it is only inferable from `e(boot_df)`+1 under `bootdf(clusters)`. *Not verified (cap).*
Overlaps practice-sweep PS-1, which verified the same fact against the comparator packages.

#### CS-69 — three different ways to name an estimand family inside `e()`, and the newest feature joined the oldest scheme

**X.** The modern surfaces spell the family out: `e(b_es_att)`, `e(V_es_ame_it)`, `e(numpre_att)`,
`e(es_ndegen_ame)`, `e(roth_power_att_it)`, `e(adversarial_pval_ame)`, `e(perm_p_att)`,
`e(joint_att_allsgssame)` (ado:4348-4360, with `f2ky` in {att, ame, att_it, ame_it}).

**Y.** The joint horizon tests use a positional legacy scheme in which `att` carries **no** token, `att_it` is
`_i`, and the continuous families are marked by an infix `C`: `e(joint_horiz)`, `e(joint_horiz_i)`,
`e(joint_Choriz)`, `e(joint_Choriz_i)` (7731-7734,
`local ename = cond(strpos("`fam'", "ame"), "joint_Choriz", "joint_horiz")`). The F065 step-down, shipped
2026-09-06, adopted that legacy scheme rather than the modern one: `e(perm_stepdown)`, `e(perm_stepdown_i)`,
`e(perm_Cstepdown)`, `e(perm_Cstepdown_i)` (7792).

**Why a user stumbles.** A user who found `e(roth_power_ame_it)` reaches for `e(joint_horiz_ame_it)` or
`e(perm_stepdown_ame_it)` and finds nothing; the `C` infix is not guessable and `_i` reads as "unit-level"
rather than "att_it". Looping `foreach f in att ame att_it ame_it` over `e()` names works for half the
surface.

**Pathways.** (a) Unify on the explicit scheme: add `e(joint_horiz_att)`/`e(joint_horiz_ame)`/... and
`e(perm_stepdown_att)`/... as the primary names, keeping the legacy names as aliases. The legacy joint names
are documented, long-standing and built by concatenation in the display code (6733, 7841, 7862), so those
sites need a lookup, and ~40 files mention them. **The step-down is only 7 files and 2 help lines — cheap to
move now, before it acquires users.** (b) Unify on the C-prefix — would mean renaming `b_es_ame` to `b_es_C`,
which collides with the honestdid handoff convention. (c) Keep both and put a family-key table in the help
("att = no token / _i = att_it / C = ame / C_i = ame_it"). *Not verified (cap).* See CS-100.

#### CS-70 — there is no stated rule for what is hidden, and the help contradicts itself about the `wdid_*` namespace

**X.** The help's Replay paragraph says the display state "travels in hidden `e(wdid_*)` macros, which are
implementation details and not part of the documented `e()` surface" (sthlp:1183).

**Y.** Two paragraphs elsewhere document `e(wdid_fcol)` (sthlp:281) and `e(wdid_fcguard)` (sthlp:364, 318) as
things the user should read to learn what the cell collapse did — and F061 (`de0d5bb`, 2026-09-06) shipped
`e(wdid_fcncell)` and `e(wdid_fcskip)` as user-facing collapse diagnostics that appear **zero** times in the
help. Conversely, plainly user-facing objects are hidden with no stated reason: `e(boot_ptails)` (9133),
`e(perm_kernel)` (7627), `e(boot_comparator)` (9211), `e(boot_ase_obs)` (9167), every
`boot_<family>_<table>` twin (8862-8865, 9006-9009, 7026-7029) and every `perm_analytic_<table>` twin (7513,
7554, 8002, 8061). And a matched pair splits: `e(perm_comparator_se)` — the observed analytic SE the
permutation studentization divides by — is **visible** (7622), while its exact bootstrap counterpart
`e(boot_ase_obs)` is **hidden** (9167).

**Why a user stumbles.** A user told the `wdid_*` namespace is off-limits will not look there for the collapse
diagnostics the author just shipped; a user reading `ereturn list` to inventory a run misses two-thirds of the
bootstrap and permutation surface, because the tables backing the printed output are all hidden.

**Pathways.** (a) Adopt and state a rule, then enforce it — "hidden = replaced-table twins, analytic remnants
and replay state; visible = anything the help names". Under it, promote `e(boot_ptails)`, `e(perm_kernel)` and
`e(boot_ase_obs)` to visible (all three are named in the help or a footer) and demote nothing; three `hidden`
keywords removed, and battery baselines that snapshot `ereturn list` move. (b) An alternative rule — "visible
= the printed numbers and their direct inputs; everything reconstructable is hidden" — would hide
`e(perm_p_lo..ci)` and `e(perm_pwsgs_kernel)`, a bigger move. (c) **Required regardless of which pathway is
chosen:** fix the help — qualify sthlp:1183 to "the display-state `e(wdid_*)` macros", and add the two missing
F061 names to the stored-results section. *Not verified (cap).* See CS-55, CS-56, CS-120.

#### CS-71 — three parallel dispersion and count objects that differ in name, in existence, or in a convention the author explicitly rejected on the other engine

**X.** `permute()` posts `e(perm_n)` (replications carrying each element), `e(perm_sd)` (SD of the null draws,
a named 1xK vector) and `e(perm_se_shift)` (the standard error **implied by inverting the reported
interval**), ado:7586-7594; the footer advertises `e(perm_se_shift)` at 10034.

**Y.** `bootstrap()` posts the same count under a different name, `e(boot_navail)` (9205-9207); it posts no
`e(boot_sd)` at all — the bootstrap standard deviation is reachable only through `diag(e(V_boot))` or by
computing it from `e(bootdist)`; and it deliberately does **not** publish an interval-implied standard error,
because the author ruled on 2026-09-04 that the inverted-interval construction "is not what the literature or
Stata's own bootstrap prefix does" (the ruling comment at 8798-8803). So the same construction is rejected on
one engine and published as a named matrix on the other.

**Why a user stumbles.** A user comparing engines on one estimate finds an SE-like number under permute with
no bootstrap counterpart, and no bootstrap `sd` vector to line up against `e(perm_sd)`, so they compare
`e(perm_se_shift)` to `sqrt(e(V_boot)[k,k])` — two constructions the author's own ruling says are not the
same object. The `n`/`navail` split defeats a simple `foreach x in perm boot` loop.

**Pathways.** (a) Unify names and add the missing one: rename `e(boot_navail)` → `e(boot_n)` (or the reverse)
and post `e(boot_sd)` as the named diagonal of the bootstrap covariance — 5 files reference `boot_navail`, 4
reference `e(perm_n)`, and adding `e(boot_sd)` is free. (b) On `perm_se_shift`, either withdraw it for
consistency with the bootstrap ruling, or state in one help sentence why the shift interval's implied SE is
publishable where the percentile interval's is not — the shift interval is symmetric by construction under
`permabsp`, so the inversion is exact there; that is the defensible distinction and it is currently unstated.
6 files reference `perm_se_shift` and the footer names it. (c) Document the mapping (`perm_n` ↔
`boot_navail`; `perm_sd` ↔ `diag(V_boot)`; `perm_se_shift` ↔ nothing, by ruling). Defensible if (b)'s
sentence is written. *Not verified (cap).* See CS-59 pathway (a).

#### CS-72 — the 32-character `e()` name cap is guarded in five places with three different behaviours and no flag saying a twin was dropped

**X.** Some twins are guarded and skipped whole: `e(perm_analytic_sgresults_#)` for a 9+ digit subgroup value
(ado:7554), `e(perm_analytic_<cfx>)` (8061), `e(boot_<ftag>_<cfx>)` (6428),
`e(boot_tpctl_sgresults_#)` (9041).

**Y.** The permute cell-draw twins are guarded per **prefix** — 9 characters for `permdist_`/`permcmp_` and 11
for `permdistse_`/`permcmpse_` (8067, 8074) — so a 22-23 character cfx name keeps the coefficient draws and
silently loses only the `permtstat` se draws, leaving a partial object; the comment there says so, and the
help documents the split. Names with no user-supplied component (`boot_wald_mainresults`, `perm_kernel`) are
not guarded at all. In every case the skip is silent: no `e()` macro records that a twin was dropped, so a
script that finds `e(boot_pctl_<longname>)` missing cannot distinguish "name too long" from "family not
computed".

**Why a user stumbles.** A user with a long continuous-treatment variable name gets a complete set of twins for
one run and a partial set for the next (after adding `permtstat`), with no message.

**Pathways.** (a) Unify on all-or-nothing per object family: guard every twin group on its longest member, so
either all four permute cell twins post or none do — users with 22-23 character cfx names lose draws they
currently get, and the help sentence documenting the split is replaced. (b) Unify on longest-prefix truncation
with a collision-safe index, as the ado already does for internal names (comment at 12030), so nothing is ever
silently dropped — real work, and truncated names are their own usability problem. (c) Keep the guards but
post one macro, `e(wdid_longnames)`, listing every twin skipped for length — turns a silent gap into a
testable one; one new macro and five one-line additions. Defensible. *Not verified (cap).* Note pathway (a)
of CS-54 tightens this budget from 22 characters to 16.

#### CS-73 — `e(boot_stat)` and `e(perm_stat)` look like the same macro and answer different questions; seed and rng-state presence differ within each engine

**X.** `e(boot_stat)` records the printed **column label** — z, t, z_wald, t_wald, t_pt (ado:9156). What was
actually bootstrapped is recorded separately, as the scalar `e(boot_tstat) = 1` (9158-9159), absent
otherwise.

**Y.** `e(perm_stat)` records **what was permuted** — `b` or `t` (192) — not what the columns say, which is
`t_perm` either way. There is no permute analogue of `e(boot_stat)` and no bootstrap analogue of
`e(perm_stat)`'s question except the differently-shaped `e(boot_tstat)`. Separately, within both engines the
reproducibility pair is asymmetric in presence: `e(boot_seed)`/`e(perm_seed)` are posted **only** when
`rseed()` was given (184-185, 9192, 407), while `e(boot_rngstate)`/`e(perm_rngstate)` are posted on every
run.

**Why a user stumbles.** A script logging `e(*_stat)` across engines to record "what statistic was used" gets a
column label from one and a design choice from the other; comparing "t_wald" with "b" is meaningless. On
seeds, a user auditing reproducibility across a batch finds `e(boot_seed)` empty on the runs where they relied
on the session seed, and has to know that `e(boot_rngstate)` is the field that is always there.

**Pathways.** (a) Unify the question: rename `e(perm_stat)` → `e(perm_permuted)` (values b|t) and add
`e(perm_stat)` as the column label, matching `e(boot_stat)`; or symmetrically add `e(boot_permuted)`. Cost:
`e(perm_stat)` is read by the display code (10027, 10032) and by 4+ test files and the help names it, so a
rename needs the footer branches updated. (b) Unify the shape of the "what was bootstrapped/permuted" flag:
replace the scalar `e(boot_tstat)` with a macro `e(boot_stat_kind)` = b|t so both engines answer with a macro.
(c) Keep both, document the difference in one sentence in each synopt, and post `e(boot_seed)`/`e(perm_seed)`
unconditionally (the literal "not set" when none was given) so the reproducibility pair is always complete —
the seed part is a one-line change with no consumers to break. Defensible in part. *Not verified (cap).* See
CS-17.

### 2.5 Display

#### CS-74 — `permute()` renames the p and se columns; `bootstrap()` renames neither

**X.** Under `permute()`, every results table's columns are rewritten: `se` → `se_analytic`, `t`/`z` →
`t_perm`, `p` → `p_perm` — ado:7503-7507 (main), 7542-7546 (subgroup), 7991-7995 (pairwise). All four
substitutions are **unconditional**, not gated on `permtstat`, so a plain permute run with a blank statistic
column still prints the header `t_perm`; sthlp:414 confirms ("the statistic column t_perm is the permuted
t-ratio under permtstat and blank otherwise"), and sthlp:1063 documents the stripe.

**Y.** Under `bootstrap()`, the p column keeps the bare name `p` and the se column keeps the bare name `se`
in every display family; only the statistic column is relabelled (8854-8855 for the `_wald` suffix on the
pctl/bc/bca twins, 8910-8913 for `t_pt` and `se`→`se_analytic` under `boottstat`, with the subgroup twins at
8996 and 9037). A grep for `subinstr local ... "p" ` returns exactly three hits — 7507, 7546, 7995 — all
permute. Analytically the stripe is `estimate se t p lb_<lvl>ci ub_<lvl>ci relyr` (3849-3850, 3976-3977,
4019-4020), so the bootstrap's headers coincide with the analytic ones on `se` and `p` exactly.

**Why a user stumbles.** A user who learns from a permute run that wooldid stamps the engine onto the column
names then reads a bootstrap table headed `p` and cannot tell from the table whether that p is a rank count
over 30 replications or an analytic two-sided t probability. Symmetrically, a user piping wooldid tables into
a script keyed on the column name `p` finds the script works analytically and under bootstrap and silently
misses the column under permute.

**Pathways.** (a) Unify on the permute scheme: rename the bootstrap's columns `p_boot` (or
`p_pctl`/`p_bc`) and `se_boot`. Cost: `wooldid_pfloor`'s `colnumb(pm,"p")` lookup at 9641 falls back to
positional column 4 (9642), so it degrades rather than errors; every pinned bootstrap cert and user script
keyed on the name moves, and the hidden per-family twins share the stripe. (b) Unify on the bootstrap scheme:
drop `_perm` from the permute columns and carry the engine only in the header and footer — loses the
at-a-glance marker permute users have, breaks permute certs, and re-creates the ambiguity F6 deliberately
removed. (c) Keep both and add one engine → column-name table to the help's Working With Output section so
the mapping is stated once.

**Defensible?** No.

**Verification.** Confirmed (code-read; the claim is about matrix colnames, set literally in source). Line
numbers run one lower than filed on the permute side; cite 8854-8855 and 8910-8913 on the bootstrap side, and
9641 for the `pfloor` lookup. The verifier added that the `t_perm` rename fires even without `permtstat`,
when the column is blank, and that the parallel subgroup relabels at 8996 and 9037 must be touched by pathway
(a). See CS-7, CS-26, CS-67.

#### CS-75 — an unavailable adversarial p prints as a test line with `.`; an unavailable joint test prints nothing at all

**X.** When the adversarial p cannot be formed, both engines post the scalar as missing — bootstrap
`ereturn scalar `ename' = .` (ado:6681, inside the family loop at 6660-6684, reached whenever `boottstat` is
absent), permute the same at 7919 and again at 7924. The display gate is
`if "`e(adversarial_pval_`typetext')'" != ""` (9806, and the subgroup copy at 9948) — a missing `e()` scalar
macro-expands to ".", which is non-empty, so the labelled line prints with `%-6.4f` of missing, i.e. a dot.

**Y.** When a joint pre-trend or equality test cannot be formed, **nothing is posted** (bootstrap 6538-6563
posts only under `if `jp' < .`; permute the same at 7752), so the line simply does not appear. Nothing stale
survives, because the analytic joint tests are gated off at the source under either engine (2509, 2989, 3297,
3608). Meanwhile the estimation-time caution for the adversarial case says the opposite of what the table
shows — "They are absent rather than shown as analytic tests, which would not be bootstrap quantities"
(6690-6691) — and the bootstrap footer repeats that claim for the joint tests (10115).

**That a missing `e()` scalar expands to "."** is settled by the codebase's own defensive idiom eleven lines
earlier in the same program: `if "`Nclus'" != "." & "`Nclus'" != ""` (9722). The adversarial gate lacks that
second test.

**Why a user stumbles.** Two tests fail for the same reason in the same run; one vanishes and one appears as a
labelled p-value of ".". A reader who takes "." for a computed-but-tiny or computed-but-degenerate number,
rather than for "not computed", draws a conclusion from a line that the caution above it says does not exist.

**Pathways.** (a) Unify on silence: add the `!= "."` guard at 9806/9948 so an unavailable adversarial test
disappears like its siblings, leaving the caution and `e(boot_adv_absent)` as the record. Cheapest, matches
the stated intent, and the scalar can still be posted as missing for programmatic callers. (b) Unify on
visibility: print a line for every requested-but-unavailable test with an explicit "(not computed)" marker
instead of a dot — more lines, nothing silently missing. (c) Keep the dot and add one footer sentence saying
a dotted test line means the test could not be formed. **(d), the verifier's addition and arguably the
cheapest:** the analytic adversarial p **is** still computed and posted under the engines (4459 is inside
`if "`adversarial'" != ""` with no engine gate, unlike the joint-test gates), so the bootstrap block
overwrites a live value with missing rather than declining to create one. Suppress the analytic adversarial p
under `bootstrap()`/`permute()` the way the joint tests already are, and the existing `== ""` gate makes the
line disappear on its own with no display code changed — making the two features structurally, not just
cosmetically, the same.

**Defensible?** No — a dot is readable as a computed number.

**Verification.** Confirmed (code read; the process cap was full, no run). The verifier also found that
`e(boot_adv_absent)` (6689) is posted but **never read** by any display or footer code — the footer at 10115
fires on `e(boot_jt_absent)`/`e(boot_sgjt_absent)` only — so on a replay the dotted adversarial line carries
no disclosure at all. See CS-42, CS-65.

#### CS-76 — "every family is computed whichever is displayed" is printed on every bootstrap run, and BCa is an exception to it

**X.** The bootstrap footer prints unconditionally: "  Every family is computed whichever is displayed:
intervals in e(boot_ci_normal)/e(boot_ci_percentile)/e(boot_ci_bc), p-values in
e(boot_p_percentile)/e(boot_p_bc), and the complete tables in the hidden
e(boot_wald_mainresults)/e(boot_pctl_mainresults)/e(boot_bc_mainresults)." (ado:10107 at the verified tree; 10112 at `2d59e70`) —
a universal claim whose enumeration silently omits BCa. The `gteffects` note makes the same claim for the
cells: "every family is stored in e(boot_*_gteffects)" (9990).

**Y.** BCa is computed only when asked for: line 365 passes
`jkacc(`=cond("`bootdisplay'" == "bca", "`wdbjk'", "")')`, so `e(boot_accel)`, `e(boot_ci_bca)` and
`e(boot_p_bca)` exist only under `bootdisplay(bca)`; footer line 10101 says so only on a bca run. The cells
never get it at all — `wooldid_bootcellsx` downgrades `bdisp` bca → bc and cautions (6355-6359) — so
`e(boot_bca_gteffects)` never exists. The help contains both halves in one paragraph (sthlp:504): "Every
family is computed and stored in e() on every bootstrap run regardless" and "because of that jackknife cost
the BCa family is computed only when it is requested, unlike the other three."

**A fourth instance the finding missed.** The user-facing `bootdisplay()` validation error at ado:6149 also
ends "...every family is computed and stored in e() whichever is displayed", and calls four legal values
"all three". So the stale claim appears in four user-visible places, not three — plus the source comments at
6123 and 8783.

**Why a user stumbles.** A user reads the footer on a default percentile run, concludes every family is in
`e()`, writes a script reaching for `e(boot_ci_bca)` as a robustness check, and gets matrix-not-found on runs
that never asked for bca. On a bca run the same user reads "every family is stored in e(boot_*_gteffects)"
under the `gteffects` table while a caution three screens up says those cells carry bc instead.

**Pathways.** (a) Fix the wording in all four places: "Three families are computed on every run (normal,
percentile, bias-corrected); BCa is computed only when bootdisplay(bca) asks for it, and never for the
gteffects/cfx cells." Pure text, no behaviour change, and it also repairs the help's self-contradiction and
the error text's miscount. (b) Make the claim true by computing BCa unconditionally — G extra fits on every
clustered bootstrap run, the reason it is on-request. (c) Say "every displayable family" and enumerate.

**Defensible?** No, as text. The on-request BCa **behaviour** is defensible and deliberate (F067).

**Verification.** Confirmed, all quotations verbatim. Corrections applied: the fourth instance at ado:6149;
and two mitigations worth recording — the always-printed footer's enumeration never names a bca matrix, so
the false part is the universal quantifier rather than a specific pointer; and on a bca run the footer at
10104 does say the cells carry the bias-corrected family (never on a non-bca run, where 9990 still promises
otherwise). See CS-11, CS-35, CS-109, CS-115.

#### CS-77 — what `e(V)` is, and what the run tells you about `test`/`lincom`, differs by engine — and the analytic run says nothing

**X.** Analytically `e(V)` is diagonal by construction (ado:4264) and the run prints nothing about it. The
only analytic `dis` lines touching `e(V)` are failure cautions (4302, 4308-4309); the statement that it is
diagonal lives solely in the help's Incompatible Post-Estimation Commands paragraph (sthlp:1186-1189).

**Y.** Under bootstrap the footer actively instructs the reader: "  e(V) holds the bootstrap covariance
(complete-replication listwise ...)" (10109) and "  e(V) and e(df_r) carry the standard-deviation family ...,
so test and lincom run by hand reproduce the bootdisplay(wald) numbers rather than the columns printed
above." (10113) — presuming `test`/`lincom` are meaningful, which they are, because `e(V)` is a genuine full
covariance. Under permute the footer says "  The standard-error column is the ANALYTIC standard error ...;
e(V) and e(V_es_*) hold the analytic covariance, so test or lincom run by hand on them are analytic, not
permutation, quantities." (10028-10029) — **without mentioning that this analytic `e(V)` is diagonal**; and
in the default fast mode there is no `e(V)` at all (10030).

**The bootstrap's sentences are themselves conditional**, so it is really four regimes: 10109 fires only when
`e(boot_v_posted) == 1` (the else branch at 10110 says the covariance is in `e(V_boot)` and none was posted
to `e(V)`), and 10113 — the explicit `test`/`lincom` instruction — fires only when the display is not `wald`
or under `boottstat`, so a plain `bootdisplay(wald)` run gets the first sentence and not the second.

**Why a user stumbles.** A user who reads the bootstrap footer's invitation to run `test`/`lincom`, then
reruns analytically or under permute and follows the same advice, gets numbers from a diagonal matrix — the
one case the help calls invalid — with no on-screen warning. Three regimes for one stored object are
announced in three registers: prose invitation, prose caveat that omits the key fact, and silence.

**Pathways.** (a) State the `e(V)` regime in one sentence on every run, analytic included ("e(V) is diagonal
on this run; use e(V_es_*) for joint work") — costs an analytic footer where there is none today, and every
analytic cert re-pins. (b) Make the permute footer's `e(V)` sentence say "analytic **and diagonal**" — a
two-word text fix that removes the worst of the three gaps immediately. (c) Leave the runs alone and
consolidate the regimes into one help paragraph the footers point at.

**Defensible?** No, for the omission in (b); the register differences follow from the engines having
different things to say.

**Verification.** Confirmed, all quoted text verbatim. Corrections applied: the bootstrap sentences are
display-conditional, making four regimes rather than three; and the help is not wholly silent on the permute
case — sthlp:418 carries the same sentence as the footer, with the same omission, so the gap is consistent
between help and footer rather than a footer-only slip. The generic diagonality statement at sthlp:1186-1189
is never cross-referenced from either engine's material, which supports pathway (c). See CS-4, CS-59.

#### CS-78 — the same "this option has no effect without that one" situation is an error in three places and a caution in two

**X, errors.** "Error: lattice() has no effect without contreat(); specify a continuous treatment variable or
remove the lattice() option." (ado:4907); every `perm*` companion without `permute()` → "Error: ...
require(s) the permute() option." (4959-4972, one error line at 4972 built from a `permstray` accumulator);
every bootstrap companion without `bootstrap()` → "Error: ... require(s) the bootstrap() option ..."
(6031-6042, likewise at 6042).

**Y, cautions, in the identical grammatical form.** "Caution: cfxlattice() has no effect without makecfxplots
and is being ignored." (5293) and "Caution: esrelativeto() has no effect without esfixedbaseperiod and is
being ignored." (5357).

**Why a user stumbles.** The user cannot predict which mistyped-companion mistake stops the run and which
costs them a long run whose output silently ignores the option they cared about. `esrelativeto()` is the
worst case: a user who wrote `esrelativeto(-3)` and forgot `esfixedbaseperiod` gets a full set of tables
computed against wooldid's **default varying base period** instead of the reference they asked for, with one
caution line scrolled far above them.

**Pathways.** (a) Unify on refusal: make the two cautions errors. Cost: breaks scripts passing a superset of
options across specifications — a real pattern for `esrelativeto()` and `cfxlattice()` in loops — and
reverses the tolerance the current behaviour encodes deliberately. (b) Unify on cautioning: downgrade
`lattice()`-without-`contreat()` to a caution. Cost: silently ignoring `lattice()` changes the estimand's
presentation, which is why it errors. (c) Keep the split and make it legible: one fixed phrasing for
"refused" and another for "ignored", and list the ignored-option set in one place in the help so the boundary
is discoverable rather than per-option folklore.

**Defensible?** Per site, yes; as a set, no. The code carries explicit in-line rationales: 5288-5291 says
`cfxlattice` "is an inert option: without makecfxplots the variable is never passed to margins, and such
calls ran fine before these screens existed - caution, do not error" — i.e. cautioning was chosen precisely to
avoid the cost pathway (a) would impose — while 4905-4906 says `lattice()` errors because "the
lattice-construction code all sits inside the contreat block, so the option would otherwise be silently
ignored".

**Verification.** Confirmed, and exercised in Stata on `dta_binary`: `lattice(3)` alone → r(198) with the
quoted message; `permreps(20)` alone → r(198) "require(s) the permute() option"; `reps(20)` alone → r(198)
"require(s) the bootstrap() option"; `esprelength(2) espostlength(2) esrelativeto(-2)` → the quoted caution,
rc 0, full run completing, with the caution printed roughly 40 lines of log before any results table (the
"scrolled far above" claim holds). A grep of `wooldid.sthlp` for "no effect without" and "being ignored"
returns nothing, so pathway (c)'s premise — the ignored-option set is undocumented — is also correct.
Corrections applied: the in-code rationales above (so "defensible" is per-site rather than false), and the
`esrelativeto` stumble is a fallback to the default reference, not a third reference. Anchor with CS-24; see
CS-43, CS-111.

#### CS-79 — three tables printed by one run use three different column vocabularies

Same fact as **CS-66** (stored-results lens) from the display side; see there for the pathways.
**Adds the user-facing contract:** the help's Tips paragraph tells users to "use svmat ... names(col)" to
remake graphs from these matrices, which makes the column names a published interface rather than an
internal detail. And under `level(90)` the main and `gteffects` tables read `lb_90ci`/`ub_90ci` while the cfx
matrices still read `ll`/`ul` — the only pair that gives no clue what level it is at. Cheapest partial fix:
add only the level suffix to the cfx bounds, so `level()` is visible everywhere it applies. *Not verified
(cap).*

#### CS-80 — the "p = 0 means p < 2/n" line exists only for the bootstrap, though `permrawp`'s floor is exactly 0

Same defect as **CS-8** (engines lens, partly confirmed there); its corrections apply — notably that
sthlp:481 does state the floors, so the gap is on-screen rather than documentary, and that the F036 low-rep
caution's `absp`-only arithmetic tells `permrawp` and `permside` users a floor they do not have.
**Adds the routine's own words:** `wooldid_pfloor`'s header comment (9629-9631) states the rationale as
"Bootstrap only: the permutation add-one p never prints 0", which was true of the default and `permabsp` and
never of `permrawp`, shipped the same day. Pathway (a) — generalise `pfloor` to key on `e(perm_pconv)` and
find the column as `p_perm` as well as `p` — is a few lines, since the routine already takes the matrix as an
argument. *Not verified (cap)* for this framing. See CS-122.

#### CS-81 — resolution disclosure fires at a different time, on a different trigger, and with incomplete arithmetic per engine

**X.** `permute()` warns **before** the run and only when the floor exceeds .05 (ado:5815-5817). It is gated
on `permexact` being empty, so an enumerated design never gets it even though the exact floor 2/M can be .167
(sthlp:481 states that number); and its arithmetic covers `permabsp` only
(`cond("`absp'" != "", 1, 2)`) although `permside()` and `permrawp` are passed into the same program (4982)
and have floors of 1/(R+1) and 0.

**Y.** `bootstrap()` gives no pre-run warning at all — `reps()` explicitly has no floor by ruling (6053-6058)
— but its footer states the resolution unconditionally **after** the run (10106).

**Why a user stumbles.** A user cautioned about resolution before a permute run expects the same courtesy
before a bootstrap run and gets nothing until the footer; a user reading the permute caution under
`permside(right)` is told the floor is 2/(R+1) when it is 1/(R+1); a user running `permexact` with one
treated cluster among twelve gets no caution at all for a hard floor of .167.

**Pathways.** (a) Give every engine both a pre-run resolution caution and a footer resolution sentence, with
the convention-correct floor. Note the `permexact` leg cannot simply be added at 5815 — M is not known there,
it is enumerated later in `wooldid_permx` (~7140), so the caution has to move after enumeration or key on the
enumerated assignment count. (b) Unify on the bootstrap practice: drop the permute caution and state
resolution in both footers — reverses F036, and the pre-run warning is the one that can still save the run's
runtime. (c) Keep both placements but fix the permute caution's arithmetic to cover `permside`/`permrawp` and
extend it to `permexact`. *Not verified (cap).* See CS-8, CS-20.

#### CS-82 — the adversarial line is the only test line whose wording changes per engine, and it changes for permute only

Same fact as **CS-16** (engines lens). **Adds the precise consequence:** because the branch at ado:9806-9808
(and 9948-9950) keys on `e(perm_method)` alone, a **bootstrap** run falls into the analytic branch, so a
`boottstat`-bootstrapped adversarial p prints with the analytic composite-on-a-line wording over a number
that is either a bootstrapped rank p (sthlp:499) or missing (CS-75). The joint lines (9789-9800), per-subgroup
joint lines (9928-9945) and equality lines (9967-9977) print one fixed analytic wording under all three
engines. Minimum fix: make the adversarial branch test the engine, not `e(perm_method)`. *Not verified (cap).*

#### CS-83 — the permutation Monte-Carlo tail record is advertised in the footer and the help; the bootstrap twin is hidden and undocumented

Same fact as **CS-56** (stored-results lens, confirmed there), whose corrections apply — including that the
hidden bootstrap record tracks the percentile family regardless of the displayed one, so pathway (a) costs
more than four renames. **Adds the display asymmetry in one line:** the permutation footer devotes a sentence
to it (ado:10032) and the help documents all four matrices (sthlp:481, 1066), while `boot_ptails` appears in
no footer and zero times in the help — so a user wanting the same reassurance for a 200-replication bootstrap
finds nothing in the footer, nothing in the help, and (because it is hidden) nothing in `ereturn list`, and
concludes the engine does not compute it. *Not verified (cap)* for this framing.

#### CS-84 — the header field labelled "Standard Errors:" carries three kinds of content, and the analytic run has no footer at all

**X.** One header line serves all three engines: "Wooldid Estimation for Outcome: y; Standard Errors:
<e(wdid_vcespec)>" (ado:9718). Analytically the field holds an option echo, `cluster(i)` (1532-1534). Under
bootstrap it holds compact metadata, "cluster bootstrap, cluster(i), reps(30)" (9150). Under permute it holds
a sentence with an option hint: "randomization inference (permute(shuffle), 30 replications); analytic
standard errors not computed (permanalyticse restores them)" (7653-7655) — a field labelled Standard Errors
announcing that none were computed.

**Y.** Below the tables, the analytic run prints **no footer at all**, the bootstrap footer runs roughly nine
to twelve lines (10086-10117) and the permutation footer roughly twelve to eighteen (10005-10080, growing by
one line per pooled estimand and one per step-down family). So the amount of inference documentation a reader
gets for the same specification is 0, ~10 or ~15 lines depending only on which engine was named. The analytic
run states its confidence level nowhere in prose — only inside the column names — while both engine footers
state it in words.

**Why a user stumbles.** A user who has read a permute footer and learned what wooldid tells them about p
conventions, interval construction and what `e(V)` means reruns the same model analytically and is told
nothing — not the reference distribution, not the level, not that `e(V)` is diagonal. The header field also
reads as a category error under permute, where "Standard Errors: randomization inference ...; analytic
standard errors not computed" answers a question the label did not ask.

**Pathways.** (a) Give the analytic path a short footer (level, reference distribution, what `e(V)` is and is
not) so all three engines close the same way — new output on the most-run path, and every analytic cert
re-pins. (b) Relabel the header field to something engine-neutral ("Inference:") and keep the per-engine
content — cheap; changes one string and every pinned header. (c) Move the permute header's prose into the
footer, leaving the field to hold a short token comparable to the bootstrap's. *Not verified (cap):* the line
counts are read from the code's branch structure, not from a run. See CS-77.

#### CS-85 — which standing disclosures survive replay depends on the engine

**X.** The permutation engine deliberately routes its standing disclosures into the footer — the source
comment at ado:7659 reads "cautions (estimation-time; the footer carries the standing disclosures)" — so
typing `wooldid` alone reprints the p convention, the shift-interval caveat, the `se_analytic` explanation,
the Monte-Carlo precision line and the pre-treatment caveat.

**Y.** The bootstrap engine leaves several standing disclosures as estimation-time cautions that replay never
repeats: the `boottstat` studentization paragraph that explains the `t_pt` and `se_analytic` columns (9170),
the bayes cell-sampling-variance caution (6220), the fewer-than-ten-clusters and fewer-than-six-treated
cautions (6260, 6279), the "cells carry bc under bootdisplay(bca)" caution (6359) and the whole
`wooldid_bootjnt` aggregator (7068). A replayed `boottstat` store therefore prints columns headed `t_pt` and
`se_analytic` with only the shorter footer branch (10093) to explain them.

**Two wrinkles in the same area.** The Roth block deliberately replays some cautions and not others (the
rank-deficiency and interrupted-pass ones are reconstructed at 9843-9852 and 9836-9840; the per-type skip
cautions at 9360-9427 are live-only, and the source says so). And the replay branch emits a blank line before
`part(perm)`/`part(boot)` (shell 47, 49) that the live path does not — so replay is not in fact
"byte-for-byte as originally shown", which is what the help's Replay paragraph promises while also saying
replay "does not repeat estimation-time Caution messages".

**Why a user stumbles.** A user who archives a run and replays it later gets a self-explaining screen under
permute and a partially unexplained one under bootstrap — the very columns whose names are least
self-evident lose their explanation.

**Pathways.** (a) Move the bootstrap's standing disclosures into `part(boot)` so they replay, keeping only
genuinely run-specific counts as live cautions — this is what the permute engine already did, and it needs
`e()` flags for the few that are data-dependent, some of which exist (`e(boot_adv_absent)`,
`e(boot_jt_absent)`). (b) Move the permute footer's standing text back into cautions — reverses the
2026-09-05 chunk-1 design. (c) Keep both and soften the help's "byte-for-byte" claim to name what replay
omits; also drop the two stray blank lines at shell 47/49. *Not verified (cap).* See CS-42.

#### CS-86 — continuous-treatment plots are `marginsplot` analytically and `twoway` under either engine, and the plot-data file appears only under an engine

**X.** Analytically the cfx plots come from `marginsplot` on the live margins results (ado:2725, 2769, both
gated `if "`wdrepmode'" == "" & "`bootstrap'`permute'" == ""`), and `saveplots()` writes only the `.gph`
(2726, 2770). The event-study plot data is written on every engine as `<name>_ES.dta` with the columns
`postperiod relyr estimate lb_<level>ci ub_<level>ci rownamesvar` (10578, 10653-10659).

**Y.** Under `bootstrap()` or `permute()` the cfx plots are redrawn by `wooldid_cfxplots` as `twoway`
rarea/connected renderings of the rewired matrices (10395-10510), and only that path writes
`<name>_cfx.dta` — with a different column vocabulary again: `matname lattice b se p ll ul` (10462-10476,
saved at 10505). So one run can emit two plot-data files naming the same quantities differently
(`estimate` vs `b`, `lb_95ci` vs `ll`), one carrying `se`/`p` and the other not, and the cfx one appears only
when an engine is named. The help documents the substitution in the `gteffects` paragraph (sthlp:967), while
the Working With Output section's `saveplots` sentence mentions only "a dataset containing the information
used to produce event study plots".

**Why a user stumbles.** A user who builds a downstream pipeline on `<name>_cfx.dta` from a bootstrap run
finds the file simply absent when they drop `bootstrap()` for a quick analytic pass; a user writing one
reader for both saved datasets has to map `estimate`↔`b` and `lb_95ci`↔`ll`; and the two figures for one
specification differ cosmetically because they came from different graph commands.

**Pathways.** (a) Always draw cfx plots through `wooldid_cfxplots` and always write `_cfx.dta`, dropping
`marginsplot` — one code path instead of two, and the analytic run already has the `e(cfx*)` matrices it
would need; analytic cfx figures change appearance for existing users and lose `marginsplot`'s conventions.
(b) Give the analytic path a `_cfx.dta` writer while leaving `marginsplot` in place, so the file's existence
no longer depends on the engine — cheapest useful step. (c) Harmonize the two saved datasets' column names on
the main-table scheme regardless of which path wrote them. *Not verified (cap).* See CS-113, CS-66.

#### CS-87 — `makeplots` without event-study horizons is skipped analytically and entered under either engine

**X.** The worker guards the event-study plot call on horizons actually existing:
`if `espre' + `espl' > 0 & "`makeplots'" == "makeplots" & "`bootstrap'`permute'" == ""` (ado:4614) — so an
analytic run with `makeplots` but no `esprelength()`/`espostlength()` draws nothing and writes no
`<name>_ES.dta`.

**Y.** Both engine blocks in the shell defer the same call and guard it on `makeplots` alone —
`if "`makeplots'" == "makeplots" { ... wooldid_esplots ... }` (permute at 258-270, bootstrap at 454-466).
`wooldid_esplots` then drops every row with `relyr == .` (10573) before deriving its graph list, so on such a
run it proceeds with an empty plot frame; and the shells treat any non-zero rc from `esplots` as an error of
the whole call, propagating it and setting the replay-staleness marker (263-269, 459-465).

**Why a user stumbles.** A user who habitually passes `makeplots` in a shared option string learns
analytically that it is harmless without event-study horizons, then adds `bootstrap()` and hits whichever
behaviour the empty frame produces — and if it is an rc, the entire bootstrap run's output is lost after the
replications have already been paid for.

**Pathways.** (a) Move the worker's `espre + espl > 0` condition into both shell guards so all three paths
agree — one line each, no behaviour change on runs that have horizons. (b) Put the guard inside
`wooldid_esplots` itself (exit quietly when no `relyr` rows survive), which fixes the hand-call surface too.
(c) **Verify first:** run analytic and bootstrap with `makeplots` and no horizons and record what each does
before choosing — if the engines merely draw nothing this is cosmetic; if either errors it is a run-losing
bug. *Not verified (cap)* — and this is the one item in the sweep where the verification would change the
severity, not just the wording.

#### CS-88 — whether the statistic column says `t` or `z` is decided by the estimator analytically and by the reference distribution under bootstrap

Twin of **CS-18** (engines lens); see there for the citations and pathways. **Adds the concrete pair:**
`bootdf()` defaults to `clusters` when `cluster()` is given (sthlp:503), so the same poisson specification
prints `z` analytically and `t_wald` under `bootstrap(standard) cluster(i)` — and a poisson user who has
learned "wooldid prints z for poisson" reads the change of letter as a change of estimator or of the poisson
handling rather than as a change of which quantity the label describes. Defensible. *Not verified (cap).*

#### CS-89 — `permreps()` below 20 is a hard error; `reps()` has no minimum at all

Same asymmetry as **CS-20** (confirmed there). **Adds the source comment** that makes it an explicit ruling
rather than an oversight: ado:6053-6058 reads "reps() is mandatory (author ruling: no floor - low-rep runs
are legitimate; see the help file on how many replications inference needs)", so `bootstrap(standard)
reps(3)` runs and prints a table while `permreps(10)` is refused with a message about burning the run.
Pathway (c) — one sentence in each option's help entry ("unlike permreps(), reps() has no minimum") — is the
only one that does not reverse a ruling. *Not verified (cap).* See CS-9, CS-44, CS-123.

#### CS-90 — failed replications are reported in the bootstrap footer and never reported under permute

**X.** Bootstrap: "  N of M replications completed" and, when any failed, " (K failed and excluded)"
(ado:10088-10090), reading `e(boot_failed)`.

**Y.** Permute: "Randomization Inference (permute(shuffle)): N of M replications completed; K possible
assignments." (10006) — `e(perm_failed)` is posted at 189 and **never displayed anywhere**, even though the
permute help says failures are common by design ("Reassigned onsets can change which periods are estimable in
a given draw ... Replications that fail to estimate are excluded and counted (e(perm_failed))", sthlp:460).

**Why a user stumbles.** The difference between completed and requested is the user's only on-screen clue
that replications failed. Under bootstrap the shortfall is named and explained on the same line; under
permute the reader has to notice N < M themselves and then guess whether the missing draws failed or were
never requested — in the engine where failure is most expected.

**Pathways.** (a) Mirror the bootstrap wording in the permute header line: "(K failed and excluded)" — one
line, no behaviour change. (b) Drop the bootstrap parenthetical and rely on N < M in both — removes a useful
signal. (c) Document that `e(perm_failed)` carries it — weakest, since the number is already computed and one
token from being printed. *Not verified (cap).* See CS-32, CS-45.

#### CS-91 — whole footer blocks exist under one engine only, and the bootstrap's contrasts sentence prints even with no subgroups

**X, permute-only.** The pooled-estimand recap, which reprints estimate / p / null band / SD for each estimand
already in the table above (ado:10068-10077), and the max-|studentized| step-down block over the post
horizons (10046-10061). Neither has a bootstrap analogue — the bootstrap surface has no multiplicity
adjustment at all and no per-estimand recap.

**Y, bootstrap-only.** The percentile-asymmetry paragraph ("A percentile interval is not symmetric about the
estimate ... this is the ordinary reading of a bootstrap table, not a wooldid peculiarity", 10099) and the BCa
jackknife line (10101).

**The one part that is not inherent.** Where permute gates its contrasts sentence on subgroups actually being
present — `if wordcount("`e(wdid_sglist)'") > 1` at 10057 — the bootstrap prints its contrasts sentence
unconditionally, hedged only by the words "when requested": "  Subgroup contrasts, when requested, are
bootstrapped as coefficient differences ..." (10116-10117). A bootstrap run with no `subgroups()` therefore
closes with a sentence about a feature it did not use.

**Pathways.** For the inherent blocks: (a) build the missing twins (a bootstrap step-down over the same
horizons; a bootstrap pooled-estimand recap) — real work, real value, new statistics to validate; (b) leave
them and note in the help which footer blocks are engine-specific. For the contrasts sentence: (c) gate
10116-10117 on `wordcount(e(wdid_sglist)) > 1` exactly as 10057 is — a one-line fix with no downside.
Defensible for the first two; not for the third. *Not verified (cap).* See CS-108.

#### CS-92 — `pretrendpower()` prints a full diagnostic block under `bootstrap()` and is refused outright under `permute()`

Same fact as **CS-95** (interplay lens, confirmed there) and **CS-19**; those citations supersede this
item's. **Adds the display-side observation:** under bootstrap the Roth block prints between `part(1)` and
`part(2)`, exactly where the analytic path puts it (ado:417-427, with the ordering comment at 411-416), so
the user sees an unbroken layout; under permute they are refused at parse and must split into two calls — the
switch the bootstrap's own fewer-than-6-treated-clusters caution (6279) recommends. Cheapest pathway: have
that caution mention that the `permute()` route it recommends does not carry `pretrendpower()`. One clause.
Defensible. *Not verified (cap).*

#### CS-93 — the same annotation class is printed in different SMCL styles, and one footer line breaks its own block's indentation

**X.** "Note:" lines print `as result` in some places — "Note: summary stats above are for untransformed Y
..." (ado:9758), "Note: all Estimates Shown Apply ... Weights" (9765, 9914) — and `as text` in others (the
three `gteffects` notes, 9988-9991). The Roth block uses a third form, lowercase and indented four spaces:
"    note: mdt_wald powers the joint pre-trend test above ..." (9876, 9892, 9896, 9898).

**Y.** Inside the permutation footer every line is `dis as result` and begins with two spaces, except the
`permtci` hedge, which is `dis as text` and starts flush left (10077). Separately, the replay branch inserts a
blank line before `part(perm)` and `part(boot)` (shell 47, 49) that the live path does not, so the two
renderings of one store differ despite the help's "byte-for-byte as originally shown".

**Why a user stumbles.** In a colour terminal the same class of annotation appears in two colours, which reads
as two severities; the one flush-left line in the permutation footer reads as belonging to the section below
rather than to the footer. Individually trivial, collectively the output looks assembled rather than designed.

**Pathways.** (a) Pick one convention per class — cautions `as text` with the "Caution:" prefix, notes
`as text` with "Note:", footer body `as result` with two-space indent — and sweep. Cost: a mechanical pass
plus cert re-pins wherever output is compared byte-for-byte, which is most of the suite. (b) Fix only the two
outliers (10077's style and indent; the two stray blank lines at shell 47/49). (c) Leave as is and record the
convention in the process notes so new lines stop diverging further. *Not verified (cap).* See CS-85.

### 2.6 Feature interplay

#### CS-94 — fast mode is refused outright by `bootstrap()` and turned on by `permute()`

Third view of the `fastnose` asymmetry, after **CS-3** (the lost `e()` surface) and **CS-23** (the ordering
trap); both were verified and their corrections apply. This item adds the **documentation** half, which is
where its own verification landed.

**X.** `bootstrap()` + `fastnose` is a hard r(198) (ado:6184), and sthlp:485 lists `fastnose` among the
incompatibilities.

**Y.** `permute()` without `permtstat`/`permanalyticse` sets it for the user (5404), and the observed run
then posts no `e(V)`, no `e(b_es_*)` and no `e(V_es_*)`. Plain `permute() fastnose` is also accepted — only
`permtstat`/`permtci` refuse it (5799-5801).

**The gap this lens found.** The `fastnose` help entry's own incompatibility list (sthlp:267) names only
`unconditionalse`, `adversarial`, `jointtests`, `impvalues`, `permtstat` and `permtci`. It omits
`bootstrap()`, refused at 6184, **and** `pretrendpower()`, refused at 5385. So a user reading the `fastnose`
entry alone learns that `permute()` turns fast mode on and never learns that `bootstrap()` forbids it — the
refusal is documented only under `bootstrap()`. And the two exemptions are stated asymmetrically: the gate at
5404 exempts `permanalyticse` as well as `permtstat`, while the sthlp:416 paragraph opens "Without permtstat"
and mentions `permanalyticse` only in the following clause.

**Pathways.** (a) Compute analytic SEs by default under permute with a `permfast` opt-out — costs the
measured margins saving that motivated the ruling. (b) Let `bootstrap()` accept `fastnose` for the observed
run, since the replications already run fast — the observed run then has no analytic `V_es` to hand
honestdid, and `boottstat` would still need it. (c) **Name the asymmetry in one place:** add `bootstrap()`
and `pretrendpower()` to the `fastnose` entry's list, state the `permanalyticse` exemption symmetrically, and
add one display line naming the dropped `e()` matrices. A two-word list fix plus one optional line; cost near
zero, breaks nothing.

**Defensible?** Yes on the behaviour, for the reason in CS-3.

**Verification.** Confirmed; both code sites and both help passages verbatim. Corrections applied:
"silently" should read "with a differently-worded header line", since ado:7657 blanks `e(wdid_fastnose)` so
the ordinary fast-mode notice at 9720 does not print and 7655 substitutes its own; the lost `e()` matrices
remain undisclosed at run time.

#### CS-95 — `pretrendpower()` is bootstrapped under `bootstrap()` and refused outright under `permute()`

**X.** `bootstrap()` supports it. ado:4595-4596 splits the call — without bootstrap the full `wooldid_ptp`
runs there; with bootstrap only the data-dependent half runs (`hbaronly`), and the shell re-enters
`wooldid_ptp` at 418 after the rewiring, passing the two `hbar` scalars back through
`e(wdid_ptphbar_att)`/`e(wdid_ptphbar_it)`, with the comment "running the diagnostic here feeds Roth's own
formula a bootstrap covariance rather than an analytic one" (the 2026-09-04 ruling lifting the plan-A13
refusal, also noted at 6296). sthlp:485: "`pretrendpower()` **is now supported**: the Roth diagnostics are
computed after the bootstrap rewiring, so Roth's own formulas are fed the bootstrap event-study covariance
rather than the analytic one, and the reported power, detectable trend and slope figures are bootstrap-based."

**Y.** `permute()` refuses it with r(198): ado:5773, "permute() is not compatible with pretrendpower(): the
Roth (2022) diagnostics are computed on an event-study variance-covariance matrix, and randomization
inference produces none", with the comment at 5767-5769 recording that "the author allowed the refusal rather
than a permutation-covariance substitute". sthlp:442 lists it among permute's incompatibilities.

**The contradiction.** `wooldid_permjnt` (prog def at 7669) **builds** a fixed metric S — "the permutation
covariance of the relevant coefficient draws (replications plus the comparator row, listwise complete)"
(7674-7677) — and inverts it for W = d' S⁻¹ d in three families: the joint pre-trend test (mata
`wooldid_permqf` at 7749), the across-subgroup equality test (7839) and the adversarial line test
(`wooldid_permadv` at 7896). So an event-study-scoped covariance over the permutation draws demonstrably
exists and is consumed by other quadratic-form diagnostics **in the same output**, while the refusal says
randomization inference produces none.

**Why a user stumbles.** The refusal's stated reason is contradicted three lines of output later by a
permutation joint pre-trend test and a permutation adversarial line test, both built on exactly the
covariance the refusal says does not exist. A user comparing the two few-cluster remedies finds one delivers
Roth diagnostics and the other declares them impossible.

**Pathways.** (a) Feed the permutation-covariance metric to `wooldid_ptp` — the plumbing exists, but it costs
a design ruling on whether Roth's formulas are meaningful on a null-relabelling covariance (they assume a
sampling covariance), plus a new caveat paragraph, and it risks shipping a diagnostic whose calibration is
unmeasured. (b) Withdraw `pretrendpower()` under `bootstrap()` too — throws away a shipped feature for
symmetry alone. (c) Reword the refusal to say the permutation covariance is a null-relabelling metric that
Roth's construction is not defined on, rather than that no covariance exists — one string, no behaviour
change. **The package already has the vocabulary:** when the bootstrap event-study covariance cannot be
built, ado:9098/9360 caution that ptp is "skipped for the affected families rather than computed on an
analytic covariance under a bootstrap header", i.e. "the covariance is of the wrong kind" as distinct from
"no covariance exists".

**Defensible?** Yes on the refusal; no on its stated reason.

**Verification.** Confirmed. Corrections applied: `wooldid_permqf`/`wooldid_permadv` are Mata functions
(13212, 13381) called from `wooldid_permjnt`, with the calls spanning 7749-7924 rather than 7671-7800; and
the `impvalues` refusal at 9249 is parallel in reasoning but not the same string. See CS-19, CS-53, CS-92.

#### CS-96 — assignment-dependent cfx lattice types are blanked under `permute()` and bootstrapped without comment under `bootstrap()`

**X.** `wooldid_permcellsx` computes a `drift` flag for lattice types whose bins are functions of the
assignment — `baseTreated`, `baseTreatedPost`, `baseTreatedPre`, `baseTreatedRef`, `byTreated`,
`byTreatedByTime`, `byTtre`, `byTtreByTime` — and blanks their inference columns: the `drift` local at
ado:8039 (gated off when `cfxlattice()` is supplied), then `M[., 4..6] = J(rows(M), 3, .)` at 8050-8052 (plus
column 3 unless `permtstat`), skipping the `permdist_`/`permcmp_` posts too (8065-8067). The program header
comment at 8026-8031 states the rationale in the finding's own terms — "lattice_k is a different dose range
in every replication" — and sthlp:436 and 967 match.

**Y.** `wooldid_bootcellsx` has **no equivalent gate at all**. Its cfx loop (6345-6432) tests only
`colsof(e(nm)) < 6` and accumulator existence; there is no `lattype`/`cfxlat` parameter in its syntax line
(6353) and no lattice-type test anywhere in the program — and the call site at 438 passes neither, unlike the
permute call site at 236, so it could not gate even if it wanted to. Its only blanking is F067's bca→bc family
swap. sthlp:967 says the cfx matrices "follow bootdisplay() and boottstat exactly as the main tables do",
with no drift caveat, and no boot-side lattice caution exists.

**The premise holds.** `wooldid_bootrun`'s header (8136-8138) describes the rep loop as "nested fastnose
wooldid calls under reweighted [aw/pw]", so each replication re-enters the lattice construction at 1190-1235,
where `gquantiles`/`xtile` recomputes the bins on that replication's weights.

**Why a user stumbles.** The same `lattice(baseTreated 5) makecfxplots` request yields a plot with blank
bands and a caution under `permute()` and a plot with full bootstrap bands under `bootstrap()`, with nothing
telling the user the second one is ranking cells whose definitions moved between replications.

**Pathways.** (a) Port the drift list into `wooldid_bootcellsx` and blank those cells with the same caution —
a dozen lines, but it removes bands users currently get, and it needs a ruling on whether bootstrap
resampling really destroys the bin definition (the quantile is a statistic of the resample, so arguably the
bootstrap is estimating its variability correctly). (b) Drop the permute drift gate — cheap, but permute's
objection is stronger, since a reassignment changes which units are "treated" by construction rather than by
sampling. (c) Keep both and say why in the help's lattice-type table: mark each type as assignment-invariant
or not and state which engine acts on that — no code change; the divergence stays but stops being invisible.

**One accidental exception worth naming in any of the three.** `latticeignoreweights` drops the weights from
the lattice computation (1194-1196), so on a weight-carried bootstrap replication the bins are **fixed** at
their original values — a real, undocumented exception to the drift objection under bootstrap, with no
analogue under permute (a reassignment changes the `ttre`-conditioned `xtileif` sample regardless of
weights).

**Defensible?** No as it stands, though pathway (a) needs the ruling above.

**Verification.** Confirmed; the whole of `wooldid_bootcellsx` was read. Corrections applied: "resamples
clusters" is loose for the default engine — the standard and bayes modes carry the resample as weights rather
than physical rows except on the narrow re-id legs, and the bins still move because the lattice quantiles are
weighted; and the `latticeignoreweights` exception above. The verifier also noted the reinforcing asymmetry
on the sibling surface: permute blanks the `e(gteffects)` cells for **all** runs while bootstrap rewires them
unconditionally (6382-6404). See CS-13.

#### CS-97 — the continuous-treatment containment screen blanks ame inference analytically, and the rewiring fills it back in

*Filed high; verifier: medium, and `defensible` flips to true — the package's own caution recommends the
engines that refill, so this is a half-disclosed design rather than a silent bypass.*

**X.** With `contreat()` and a clustered VCE, the AMESE R1 containment predicate fires and every downstream
`ame` capture nulls its se/t/p/ci cells while keeping the point estimate (`amescrub` set at ado:1599-1610,
`amescrublo` defaulting to 2 at 1761, applied at 2558-2560, 2611-2613, 3510-3513, 3529-3532, 3556-3559,
3664-3667, 3708-3711; `e(V_es_ame*)` withheld at 4353). The harsher `amevoid` route additionally blanks the
estimate with a printed caution (1872-1877). The screen is skipped by design under `impvalues` (comment at
1594-1598) and under `fastnose`/`unconditionalse`.

**Y.** Both rewiring stages decide what to exclude by testing the **estimate** cell only —
`if `wm'[`rr', 1] >= . matrix `excl'[1, `kk'] = 1` at 8672/8681 (`wooldid_bootwire`) and 7442/7452
(`wooldid_permwire`) — and then write into every non-excluded row (8810-8846). Under the R1 route the estimate
is present, so the row is refilled: `bootstrap()` returns full se/t/p/ci for exactly the rows the analytic
engine declined to report. `permute()` refills the p-value and the shift interval only — the SE column is
left as the worker posted it, blank under the fast-mode default and analytic only under
`permtstat`/`permanalyticse`. The `amevoid` route is immune, because it blanks the estimate and so sets
`excl`.

**Why a user stumbles.** The user runs an analytic specification, is told the ame standard errors are being
withheld because the design cannot support them, adds `bootstrap()` — and gets a full ame inference column
back, with nothing at the point of refill saying which rows now carry resampled inference the analytic engine
declined to report.

**What is not true.** The screen is **not** silently bypassed. `wooldid_amescreen` prints an unconditional
caution at 10276 naming the contained cells and saying "Consider switching to use of bootstrap SEs, IM
p-values, or randomization inference if finer clustering levels are non-viable" — and under `bootstrap()` the
observed pass is the ordinary noisy worker run, so the user sees that caution in the same run. The package
therefore tells the user to reach for the engines that refill. The residual defect is disclosure quality, in
two places: that caution's text is stale when `bootstrap()`/`permute()` is already active ("consider
switching to" what you are already using), and nothing at the point of refill marks the affected rows.

**Pathways.** (a) Carry the containment flag into the rewiring's exclusion vector so blanked rows stay blank
under all three engines (post `e(wdid_amescrub)` and OR it into `excl`). Costs more than filed: it now
contradicts the package's own printed advice, and it withholds bootstrap ame numbers that may be sound, since
a cluster bootstrap does not rely on the within-cell variance the analytic screen protects against. (b) Rule
that a resampling engine legitimately supersedes the screen and say so where the screen fires — de facto the
package's stance already, so this reduces to making the 10276 string engine-aware. (c) Disclose at the point
of refill: count the refilled rows and print a caution naming them. Small, auditable, and the cleanest live
fix.

**Defensible?** Yes, on the verifier's reading — deliberate and half-disclosed, not a bypass.

**Verification.** Partly confirmed (code-verified; the box was at the cap). Corrections applied: the line
cites (8672/8681, 7442/7452, write at 8810-8846); permute refills only the p and the interval, not a full
inference column; `amevoid` is unaffected; and the 10276 caution refutes "no indication that a correctness
screen was bypassed" and "the only hint is a boottstat caution" (the 6433 caution is a different surface
entirely — gteffects/cfx cells with no analytic SE to studentize).

#### CS-98 — `contreatwithin` prints "wooldid will not produce these plots" and then produces them

**X.** The help presents `contreatwithin` as a narrow relaxation: "This option disables some of *safety*'s
features (the estimability-relaxation checks; the data-integrity screens, such as onset constancy within
units, still run unless *safety(off)* is given explicitly)" (sthlp:739).

**Y.** The code sets the working safety local to off wholesale — `local safety off` at ado:5157-5162, with
the comment "the above is required since we are intentionally making the base effects irrelevant" — and that
flipped value reaches the worker through `wooldid_validate`'s blanket `c_local` return list (5455). Only the
checks rewritten to gate on `safetyinput` (captured at 4989) survive. Two gates on the plain `safety` local
therefore change, and the second is the sharp one:

1. **The amevoid screen at 1760** (`if "`safety'" != "off"`), which carries the caution that every treated
   dose-slope term was dropped and "the within-cell average marginal effects (ame) that wooldid reports for
   continuous treatments are NOT IDENTIFIED in this design" (1872). Under `contreatwithin` the whole block is
   skipped, so `amevoid` stays at its pre-initialized 0 (1758-1759) and the ame results are reported as the
   structural zeros the help says only `safety(off)` users get — silently. `contreatwithin` is precisely the
   design in which within-cell dose variation *is* the identification, so it is the design where "every
   treated dose-slope term was dropped" most needs saying, and it is the one setting where wooldid will not
   say it.
2. **The non-`att` cfx plot gate at 2655** (`if "`makecfxplots'" != "" & ("`item'" == "att" | "`safety'" ==
   "off")`). This is worse than "no message connecting the two": the caution that announces the refusal sits
   at 4912-4920, **earlier in `wooldid_validate` than the flip at 5157**, so it still reads the user's raw
   safety value and still prints the refusal.

**The run, both ways.** On `test/data/dta_cont.dta`,
`wooldid y i t ttre, contreat(intensity) att_i contreatcontrols(basic) makecfxplots lattice(all 2)`:

- without `contreatwithin`: prints "As such, wooldid will not produce these plots unless you specify
  safety(off)."; `graph dir` afterwards is empty.
- with `contreatwithin`: prints the **same sentence**, and `graph dir` afterwards lists
  `cfxame_i_cfxpost  cfxatt_i_cfxpost`.

So wooldid tells the user in one breath that it will not draw these plots and draws them in the next.

**Why a user stumbles.** The user reads that `safety(off)` is what turns off the multicollinearity and
identification screens, never types it, and gets a run in which the exact identification failure the package
advertises it catches passes without a word — and a run whose printed message contradicts its own output.

**Pathways.** (a) Unify on `safetyinput`: move the remaining `safety`-gated checks onto it and give
`contreatwithin` a named, enumerated relaxation list. The audit is real work but smaller than it looks — the
20 plain-`safety` comparisons in the file span the shell and display subprograms where `safety` is the user's
raw value and the flip is irrelevant. (b) Keep the flip but announce it: print, under `contreatwithin`, the
list of screens it turned off — one display line built from a fixed list, no behaviour change. (c) Correct
the help, which currently implies only the estimability-relaxation checks are affected. **Independently of
all three:** gate the 4912-4920 caution on the post-rewrite value, or move the rewrite above it, so the
self-contradiction goes away.

**Defensible?** No.

**Verification.** Confirmed, with a live run
(`/Users/hegland/.claude/jobs/be9644d0/tmp/consist/skeptic-cw/t.log`) for the cfx-plot half; the amevoid half
is code-read only, since the fixture's ame is identified so the screen does not fire there. The verifier added
the propagation mechanism (the flip is in `wooldid_validate` and reaches the worker only via the `c_local`
loop at 5455) and softened the "about twenty gate sites" estimate in pathway (a). See CS-37, which verified
the *maxvar* half of the same flip and narrowed it, and CS-25 sense (7).

#### CS-99 — the `jointtests` screen is a hard error analytically and is waived without a word once an engine is added

Same phenomenon as **CS-41** (policy lens, partly confirmed) and **CS-15**; those citations apply. **Adds the
display half:** the whole screen is gated on `"`bootstrap'`permute'" == ""` (ado:1498), and **nothing is
printed** to say the screen was waived — the joint-test display block (9786-9801) prints only the p-values,
and no "exempt"/"screen" text appears anywhere in the display path. The waiver is documented only in the help
("That combination is exempt from the analytic joint-test screen described under jointtests - bootstrap
inference is what that screen recommends"). Replay gives no trace either. Pathway (c) is the cheap one: one
line under the joint-test block when an engine is on and the cluster nests `i` ("the analytic joint-test
screen does not apply here; these are bootstrap/permutation tests"), plus an `e()` flag for replay.
Defensible. *Not verified (cap)* for this framing.

#### CS-100 — the continuous-treatment twin of a stored result is named three different ways in one run

Same naming split as **CS-69** (stored-results lens); see there for the pathways and the F065 timing point.
**Adds the display contradiction, which is this lens's own contribution:** the line printed for
`e(joint_Choriz)` reads "p-value for joint test that all pre-treatment **ame** horizon estimates are 0"
(ado:9796) — the word `Choriz` never appears on screen. Meanwhile the across-subgroup equality tests post
`e(joint_ame_allsgssame)` via `local itemtext = subinstr("`item'","att","ame",.)` (4499-4501), the adversarial
tests post `e(adversarial_pval_ame)` (4458-4460), the event-study objects are `e(b_es_ame)`/`e(V_es_ame)`, and
the pairwise-contrast rows of `e(pwsgs)` are labelled `ame:` (4009). So a user who reads "ame" on screen and
writes `e(joint_ame_horiz)` gets nothing, while the sibling equality and adversarial tests in the same footer
do use the `ame` spelling. Cheapest mitigation short of renaming: print the `e()` name beside each footer
p-value. *Not verified (cap).*

#### CS-101 — the same joint test is an `e()` scalar overall and an `e()` macro per subgroup, while the adversarial p is a scalar in both

Same fact as **CS-60** (stored-results lens, confirmed there), whose correction matters for the pathway: the
display block at 9928-9942 dereferences the subgroup values as macros and must be edited under pathway (a),
while the `allsgssame` existence idioms at 6733/7862 need no change. **Adds the within-run juxtaposition:**
the adversarial p is posted three lines below in the same subgroup loop and is `ereturn scalar` at both
scopes (ado:4459, 4483), so inside one subgroup run `e(joint_horiz_s1)` must be read as a string and
`e(adversarial_pval_att_s1)` as a number, and the emptiness tests differ (`"`e(joint_horiz_s1)'" != ""` versus
`e(adversarial_pval_att_s1) < .`). Anyone looping over subgroups to build a results table writes one branch
and gets half the cells. *Not verified (cap)* for this framing.

#### CS-102 — `gteffects`, the cfx matrices and the cohort histogram get three different treatments under the same engine

**X.** Under `bootstrap()`, the `gteffects` cells and the `e(cfx*)` matrices are fully rewired to bootstrap
inference from per-cell accumulators (ado:6381-6401, 6404-6432), and the display note says so ("under
bootstrap() the se, p, lb and ub columns are the ... bootstrap family of each cell's own replications",
9990).

**Y.** Under `permute()` the same `gteffects` cells are blanked (7558-7572, notes at 9990 and 10064) while
the cfx matrices are permuted for invariant lattice types. Under `impvalues` or `unconditionalse` **neither**
surface changes at all — the display falls through to "standard errors are the stage-1 regression's
coefficient-level VCE" (9991) even though the main table is IM- or unconditional-based. And
`histogramcohorteffects` is untouched by every engine, with a caution under both ("histogram outputs are
drawn from the stage-1 ANALYTIC covariance", 274 and 475) — although the histogram margins calls all carry
`nose` (2168, 2172, 2198, 2205), so `e(histogramestimates)` has no inference columns for any covariance to be
wrong about; what actually differs is that its cohort points are the observed-run points.

**Why a user stumbles.** Four cell-level surfaces that all read "the same estimates, disaggregated" behave
four ways under one option, and only two of the four say which. A user who checked that `gteffects` was
bootstrapped assumes the histogram was too; a user under `impvalues` has no on-screen signal that the
disaggregated SEs are on a different footing from the aggregate ones.

**Pathways.** (a) Unify the **disclosure** rather than the behaviour: emit one "inference provenance" line
per auxiliary surface on every run (analytic / bootstrap-<family> / permutation / blank), driven from a single
table — modest display work, no numbers change, and it makes replay self-describing. (b) Unify the behaviour
where possible: bootstrap the histogram too (one more per-replication margins read) and blank it under
permute for the same key reason `gteffects` is blanked — real replication cost, and the histogram is a
descriptive plot that arguably needs no inference. (c) Keep as is and fix only the histogram caution's
wording, which names a covariance the histogram never consumed. Defensible. *Not verified (cap).* See CS-13,
CS-50.

#### CS-103 — `safety(off)` sometimes changes the answer and sometimes only silences the sentence about it

Same taxonomy as **CS-25** and **CS-40**, both verified; their corrections apply. **Adds the sharpest pairing
in one place:** with `subgroups()` and `esprelength() > 0`, `espostlength(0)` is silently raised to 1 and only
the caution is gated on safety (ado:4936-4941), so under `safety(off)` the estimand window is still widened
and only the notice disappears — while elsewhere `safety(off)` changes what is computed, as when a subgroup
variable that splits a unit proceeds under **share weighting**, so `att_i` and `att_itime` weight a
straddling unit by its weight share rather than equally (1118-1127) — a different estimand. The help's model
of `safety()` is "caution or halt" (sthlp entry), so a user reads `safety(off)` as "stop stopping me" and is
surprised either way: they set it to silence one message and unknowingly change the `att_i` estimand
elsewhere, or they set it expecting `espostlength(0)` to be honoured and it is not. Pathways as CS-40, plus:
give every safety-gated behaviour change its own explicit opt-in, as `espostlength(-999)` already is for this
one. *Not verified (cap)* for this framing. See CS-34, CS-51.

#### CS-104 — the cell collapse's degeneracy screen guards the binary subgroup reads and vanishes as soon as `contreat()` is added

**X.** Every binary subgroup margins call is followed by `wooldid_fcdeg`, which inspects the returned
variance for a not-estimable marker or a singular contrast direction and, under default safety, restores the
full data and re-issues the call. The five sites are ado:3094 (by-subgroup event-study horizons), 3303 and
3321 (the subgroup joint pre-trend tests) and 3364 and 3419 (the pairwise contrasts), each wrapped in
`forval fcpass = 1/2` with `restore, preserve`. The outcome is recorded in `e(wdid_fcguard)` and carried into
bootstrap replications.

**Y.** The nine continuous-treatment subgroup margins calls carry no screen at all: 3501, 3520, 3547, 3611,
3624, 3655, 3686, 3694 and 3732 all run `contreatelasticitytype(contreat) ... over(subgroups)` with no
`wooldid_fcdeg` after them. The help states the gap ("the continuous-treatment subgroup reads under contreat()
are not screened; see the design record") but the `fastcollapse` eligibility paragraph still lists
`subgroups()` as "eligible, but screened", which is true of the binary half only.

**Why a user stumbles.** A user who has read the long `fastcollapse` note comes away believing `subgroups()`
is protected, and the protection silently halves when they add a dose variable — in the branch where
cell-level singularity is most likely, since the ame reads depend on within-cell dose variation that the
collapse replaces with a cell mean.

**Pathways.** (a) Call `wooldid_fcdeg` after the nine dose reads too — a few matrix operations per call, so
the cost is small, but the F012 research record (2026-09-06) is a negative result: the ame-singular /
att-healthy case exists in 5 of 10 fixtures and is **inert**, no ame quantity moving beyond 7.5e-13 between
collapsed and full-data arms, with the caveats that the probe used two subgroups, no weights and balanced
data. So (a) would buy an unmeasured benefit and make some `contreat` subgroup runs slower for no
demonstrated gain. (b) Drop the binary screen as well — unacceptable; it was built on measured cases where
the collapse changed what was reported. (c) Say it in the eligibility paragraph rather than only in the
aside: change "subgroups() is eligible, but screened" to name the binary/continuous split, and record it in
`e(wdid_fcguard)`'s documentation. Defensible — and this is the standing choice the author already framed as
"document, or extend the screen for symmetry". *Not verified (cap)* here; the site list was independently
confirmed by the F012 record.

#### CS-105 — `impvalues` refuses eight partner features by error but not the one the help says gives incorrect results

**X.** Ibragimov-Muller inference is refused outright with `subgroups` (ado:4777-4780), `jointtests`
(4804-4806), `makecfxplots` (4792-4794), `adversarial` (4727-4730), `unconditionalse`/`customvce`
(4751-4754), no cluster or multiway cluster (4740-4748), `fastnose` (5372-5374), `pretrendpower` (9249) and
both replication engines (5733, 6180) — nine error sites in all.

**Y.** There is **no gate anywhere** on `impvalues` combined with `att_i`, `att_it` or `att_itime`, even
though the help states the combination is wrong: "wooldid does not adjust the weights used for computing
different estimands when IM-style inference is requested: in general, IM-style inference will yield incorrect
results for non-att estimands if clusters cross-cut cohorts or (worse) cohort-time period cells." A grep of
`impvalues` against `att_i`/`att_it`/`att_itime` finds only the shell pass-through at 4996. The run proceeds
silently and prints those estimands in the ordinary table.

**Why a user stumbles.** The refusal wall teaches that wooldid polices `impvalues` combinations hard, so
silence on `impvalues att_i` reads as endorsement. The condition that makes it wrong — clusters crossing
cohorts — is exactly what the code already computes elsewhere for the collapse key and the nesting checks, so
it looks checkable.

**Pathways.** (a) Add a nesting test (the cluster-nests-`i` assert already exists at 5904-5912 and
6238-6247) and refuse or caution `impvalues` with cell-weighted estimands when the cluster crosses `i` —
twenty lines of existing code, and it would break runs that are currently fine when the cluster does nest.
(b) Downgrade some of the nine refusals to cautions for symmetry — wrong direction; several guard
r(321)-class failures. (c) Put a caution where the estimand is printed, naming the crosscutting condition —
a single conditional caution, no behaviour change. *Not verified (cap).*

#### CS-106 — `subgroups()` screens are enforced over the treated sample analytically and over the whole design population under `permute()`

Same fact as **CS-52** (policy lens); see there for the pathways. **Adds the two permute-only cautions that
accompany the widening:** ado:5951 (share weighting under `safety(off)`) and 5954 (the shuffled assignments
also re-randomize which subgroups the treated clusters fall in, when `permstrata()` is absent), so the
engine both tightens the screens and adds disclosure the analytic path has no need for. The reason — any
cluster can be assigned treatment in a replication — is stated in every message, which makes this the
best-handled instance of a widened screen in the package; the stumble is expectational. Defensible. *Not
verified (cap).*

#### CS-107 — the two replication engines accept different partners for the same stated reason

**X.** `bootstrap()` refuses `customvce()` (ado:6192) on the ground that "the reported standard errors would
be bootstrap-based, contradicting the custom-VCE request", and refuses `boottstat` with `poisson` (6140)
because "the Poisson replications are estimated on sampling weights, so a replication's own standard error is
not the same functional of the resampled data as the observed standard error".

**Y.** `permute()` has no `customvce()` gate at all — `wooldid_permctl` (5722-6000) never mentions it and the
shell forwards it into every replication (129) — even though under `permute()` every reported p-value and
interval is a permutation quantity, so the custom VCE is overridden exactly as the bootstrap says it would
be. And `permtstat` is accepted with `poisson`: its only gates are `fastnose` (5799) and the cluster-at-`i`
caution (5902).

**Why a user stumbles.** The user reads the two "Not compatible with" lists side by side in the help and
cannot recover a rule: the bootstrap names `customvce()` and permute does not, though permute's shift
intervals replace the custom VCE just as thoroughly; the bootstrap forbids studentizing under poisson and
permute allows it. Switching engines to compare few-cluster remedies changes which of their other options
survive.

**Pathways.** (a) Unify on the bootstrap's list: refuse `customvce()` under `permute()` as well — one gate,
and it removes a combination that currently runs and is deliberately forwarded to replications. (b) Unify on
permute's: accept `customvce()` under `bootstrap()` with a caution that it governs only the hidden
`e(V_analytic)` — reverses a shipped refusal. (c) Give the reasons per engine: the poisson split is genuinely
principled (permute does not resample, so the same-functional argument does not apply) and deserves a
sentence saying so beside both option entries; the `customvce` split looks like an oversight and should at
minimum be stated. Defensible for the poisson half. *Not verified (cap).* See CS-48, CS-127.

#### CS-108 — multiplicity across the post-treatment horizons is adjusted under `permute()` and not under `bootstrap()`

**X.** `permute()` computes the max-|studentized| step-down of Westfall-Young / Romano-Wolf over the
post-treatment horizons for each family and subgroup, unconditionally whenever two or more post horizons
exist — "computed whenever at least two post horizons exist (jointtests is not required)" (ado:7783-7800) —
and posts `e(perm_stepdown)`, `e(perm_stepdown_i)`, `e(perm_Cstepdown)`, `e(perm_Cstepdown_i)` plus a footer
summary (10041-10044).

**Y.** `bootstrap()` has no step-down at all: `wooldid_bootjnt` (6503-6700) computes joint pre-trend tests,
equality tests and adversarial tests and nothing else, and no `e(boot_stepdown*)` exists. The bootstrap
replication draws are exactly the object the same construction needs, and the fixed-metric joint tests shipped
2026-09-06 already recentre and studentize them (6513-6567). The help's permute entry has a "Multiplicity
across the post-treatment horizons" paragraph with no bootstrap counterpart anywhere in the file.

**Why a user stumbles.** A referee asks for family-wise error control across event-study horizons. The user
finds it in the `permute()` help, switches to `bootstrap()` because their design has enough clusters, and the
column disappears — not because the bootstrap cannot support it but because it was built on one engine. The
two engines otherwise present near-identical joint-test menus, which makes the gap look like a mistake in the
output rather than a missing feature.

**Pathways.** (a) Build the same step-down from the bootstrap draws (recentred at the observed estimates,
studentized on the fixed metric `wooldid_bootjnt` already computes) — real work, and it needs a ruling on
whether a bootstrap step-down under resampling has the exactness the permutation version claims under the
sharp null. (b) Withdraw the permute step-down — loses a shipped feature. (c) One sentence in the
`bootstrap()` help entry stating that the horizon step-down is a permutation-only surface and why — it stops
the gap reading as an omission in a given run's output. Defensible. *Not verified (cap).* See CS-91.

#### CS-109 — `bootdisplay(bca)` gives the coefficient rows a BCa interval and the `gteffects`/cfx cells a bias-corrected one

**X.** The coefficient tables under `bootdisplay(bca)` carry the BCa family, with the acceleration from a
delete-one-cluster jackknife (`e(boot_ci_bca)`, `e(boot_accel)`, `e(boot_jk_completed)`); the footer explains
the jackknife (ado:10104).

**Y.** `wooldid_bootcellsx` rewrites `bdisp` to "bc" before it starts —
`if "`bdisp'" == "bca" { local bdisp bc; dis as text "Caution: the gteffects/cfx cells carry the
bias-corrected (bc) family under bootdisplay(bca)..." }` (6356-6359) — so within one run the aggregate rows
and the cells they aggregate carry different interval families, and `e(boot_gteffects_family)` records "bc"
while the main tables are BCa.

**Why a user stumbles.** The weighted-mean identity the `gteffects` help promises ("the weighted mean of
estimate under any of these weight columns reproduces the corresponding aggregate exactly") holds for the
point estimates but not for the intervals, and the reason is a family switch announced only in a caution that
scrolls past. A user reconstructing an aggregate CI from cells gets a different number.

**Pathways.** (a) Unify on BCa by jackknifing every cell — G extra fits per cell surface, which the F067
ruling already judged too expensive. (b) Unify on bc for everything under `bootdisplay(bca)` — defeats the
option. (c) Keep both but make it durable rather than transient: post the mixed state in a visible `e()`
macro and repeat it in the `gteffects` display note, not only as a one-shot caution — a few lines, and
replayed results then disclose it. Defensible (recorded as judgment call 2 in the ship record). *Not verified
(cap).* See CS-76, CS-85.

#### CS-110 — the binary treatment's transformation is a bare flag with one setting; the dose's is a four-valued option, with different labels and different collapse eligibility

**X.** For the binary treatment the only transformation is `semielasticity`, a bare flag resolving to `eydx`
(ado:5047-5052). It labels its output "Semi-elasticities" with no qualifier (2149) and always turns the cell
collapse off (5107).

**Y.** For the dose the same idea is `contreatelasticitytype(dydx|eydx|dyex|eyex)` (5046, domain check
5417-5420). Its labels carry the operator — "Semi-elasticities (eydx)", "Semi-elasticities (dyex)",
"Elasticities" (2150-2152) — and its collapse eligibility splits by operator: `eyex` and `eydx` turn the
collapse off (5075) while `dydx` and `dyex` stay eligible, `dyex` having left the ineligible list under
poisson on 2026-09-06 (5076-5078, F032). There is no `dyex` or `eyex` for the binary treatment at all.

**Why a user stumbles.** A user reporting `semielasticity` results who moves to a continuous treatment cannot
find the flag; one who learned `contreatelasticitytype(eydx)` cannot ask for the binary analogue by that
name. Reading a table, "Semi-elasticities" and "Semi-elasticities (eydx)" are the same quantity on the two
sides of the run, and the collapse silently applies to one and not the other.

**Pathways.** (a) Unify on the four-valued form: add `elasticitytype(dydx|eydx)` for the binary treatment
(only those two are meaningful for a 0/1 regressor) and keep `semielasticity` as a synonym for
`elasticitytype(eydx)` — modest parse work, and the label gains the "(eydx)" qualifier, changing printed
output. (b) Unify on the flag — worse; four operators do not reduce to flags. (c) Cross-reference in both
help entries and drop the "(eydx)" qualifier so the two sides print the same label — two help edits and one
string. Defensible. *Not verified (cap).*

#### CS-111 — a companion option supplied without its parent is an error in some families and an ignored caution in others

Same taxonomy as **CS-24** (options lens, partly confirmed) and **CS-78** (display lens, confirmed with a
live run); their corrections and evidence apply. **Adds the two clearest same-shape pairs:**
`poisexpresults` without `poisson` errors (ado:5433-5436, recorded at 5433 as author ruling SI-9) while
`esrelativeto()` without `esfixedbaseperiod` cautions (5354-5358), and both are inert in exactly the same
way; likewise `bootalpha()` with `bootstrap(standard)` errors (6092) while `cfxlattice()` without
`makecfxplots` cautions (5290-5293). A scripted sweep that tolerates cautions dies on one branch and not the
other. Pathways as CS-24; note that (b) would reverse the shipped SI-9, VAL-G19 and NFADV-07 rulings, all
made to close silent-accept paths. Defensible per site. *Not verified (cap)* for this framing.

#### CS-112 — the same `esfixedbaseperiod` prerequisite degrades one pre-trend diagnostic with a caution and refuses its sibling

**X.** `adversarial` under the pooled reference regime is not refused: "Caution: the adversarial pre-trend
p-value requires esfixedbaseperiod (the pooled reference has no single omitted period); the fitted line is
reported without a p-value" (ado:1055). The run continues and `e(adversarial_line_*)` is still posted
(4453-4461, whose comment reads "independent gates - the display line survives when inference is
unavailable").

**Y.** `pretrendpower()` with the same prerequisite unmet is a hard r(198): 9253-9255, "pretrendpower()
requires esfixedbaseperiod. Under the pooled reference the comparison period is a cohort-specific pool of
pre-periods, so the single omitted reference period that the Roth (2022) linear-trend normalization requires
is not defined." The reason given is the same one the adversarial caution gives.

**Why a user stumbles.** Both options are pre-trend diagnostics listed together in the Event Study section,
both need a single omitted reference period, and both are typically added in the same edit. Adding them
together to a pooled-reference run kills the run; adding only `adversarial` produces a partial answer.
Nothing in the help pairs the two prerequisites so the user can see they are the same condition.

**Pathways.** (a) Unify on the caution: let `pretrendpower()` print a caution and skip, as it already does
when the pre-period covariance block is unavailable (its own help says it is "skipped with a caution when the
pre-period variance-covariance block is unavailable") — a few lines, consistent with its other failure mode,
and it removes a fail-fast that currently protects the user from a meaningless request. (b) Unify on the
error: refuse `adversarial` under the pooled reference — loses the display-only fitted line the code
deliberately keeps. (c) Pair them in the help: one sentence in each entry saying the other option shares the
prerequisite and reacts differently — no code change. Defensible. *Not verified (cap).* See CS-38, which is
the same option's *other*, mute cancellation.

#### CS-113 — `saveplots()` writes a different set of files, drawn by a different command, depending on the inference engine

Same fact as **CS-86** (display lens); see there for the mechanics and pathways. **Adds where the user would
look:** the help admits the difference in the `gteffects` entry ("as twoway renderings of the matrices rather
than marginsplot output, so they can differ cosmetically from an analytic run's plots, and saveplots()
additionally writes one name_cfx.dta with every plotted point") and **not** in the `saveplots()` entry, which
is where a user assembling a paper would look. So an analytic and a bootstrap version of the same figure do
not sit side by side in a draft, and the file inventory changes with the engine, with no warning in the
option's own documentation. Cheapest fix: list per engine, under `saveplots()`, which files each writes and
note the renderer switch. Defensible. *Not verified (cap).*

### 2.7 Documentation

#### CS-114 — Example 13 still says bootstrap joint tests are absent without `boottstat`; the `bootstrap()` entry says they are the fixed-metric rank test

**X.** The Examples section, example 13, tells the user that `jointtests` under `bootstrap()` requires
`boottstat` and that without it the joint tests do not exist: sthlp:1296, "13) Cluster bootstrap standard
errors with 500 replications, reproducible, and joint pre-trend tests formed from the replications' own test
statistics ({cmd:boottstat}; without it the joint tests are absent under {cmd:bootstrap()})".

**Y.** The `bootstrap()` option entry, rewritten when F034 shipped (`bd8ee3c`, 2026-09-06), says the joint
pre-trend tests and the across-subgroup equality tests are bootstrap quantities in **both** forms, with
`e(boot_jt_family)` recording which (sthlp:485). The entry, not the example, matches the code: ado:6570 sets
`ereturn local boot_jt_family fixedmetric` on the non-`boottstat` path, 6639 sets `studentized`, and 10114
prints a footer describing the fixed-metric joint tests. The `boottstat` entry itself ends "Without
{cmd:boottstat} those tests are the fixed-metric rank tests described under {cmd:bootstrap()}", so example 13
is contradicted by two option entries, the stored-results section (sthlp:1061) and the code.

**Why a user stumbles.** The Examples section is where most users start. A user who wants bootstrap joint
tests reads example 13, concludes `boottstat` is mandatory, and pays the full studentized cost — physical
cluster re-identification plus per-replication variance estimation, both documented as materially slower —
for a test the default path now delivers. Worse, `boottstat` carries a documented degeneracy warning for
pairwise contrasts under clustering at or above `i` (CS-41), so the example steers users into a configuration
the `boottstat` entry itself warns against.

**Pathways.** (a) Make the option entry canonical and rewrite example 13's parenthetical to say `boottstat`
chooses the *studentized form* of the joint tests rather than enabling them — cheapest, no code change, but
it leaves example 13 without a reason to show `boottstat` at all, so the example needs a new motivation (e.g.
show it as the percentile-t route). (b) Make the example canonical by reverting F034 — costs the shipped
feature and contradicts the author's own ruling. (c) Delete the parenthetical and cross-reference
`bootstrap()`'s joint-test paragraph — least text, but it loses the genuinely useful signal that two different
joint tests exist.

**Defensible?** No — a stale doc claim.

**Verification.** Confirmed; both quotations verbatim, and `review_process/qr_ship_2026-09-06.md:28` records
`bd8ee3c` as shipping exactly this. Note the stumble is a needless cost, not a wrong number: a user following
example 13 gets correct output, just slowly.

#### CS-115 — `bootdisplay()` says every family is computed on every run, then says BCa is computed only when requested, and the code's own error text repeats the stale claim

**X.** The opening sentence of the `bootdisplay()` entry promises that all inference families are always
computed and stored, so the option only selects a display: sthlp:504, "Every family is computed and stored in
{it:e()} on every bootstrap run regardless; this option chooses only what is displayed".

**Y.** Four sentences later the same entry says "because of that jackknife cost the BCa family is computed
only when it is requested, unlike the other three." The Stored Results block agrees with the second version
("under {cmd:bootdisplay(bca)} also {it:e(boot_ci_bca)} ...", sthlp:1061), and the code gates the jackknife
on `bootdisplay(bca)`: ado:365 and 394 pass `jkacc(`=cond("`bootdisplay'" == "bca", "`wdbjk'", "")')`, 8573
builds the jackknife matrix under `if "`jkacc'" != ""`, 8726 sets `local wdbbca = ("`jkacc'" != "")`, and
every BCa `ereturn` (8865, 9009, 9136-9145, 7029) sits under `if `wdbbca''. So `e(boot_ci_bca)`,
`e(boot_p_bca)` and `e(boot_accel)` genuinely do not exist without `bootdisplay(bca)`.

**And the user-facing error text repeats the stale claim.** ado:6149: "The reported standard error is the
standard deviation of the bootstrap distribution under **all three**; every family is computed and stored in
e() whichever is displayed" — "all three" and "every family" against the four legal values accepted one line
above at 6148.

**Why a user stumbles.** A user who reads the first sentence writes a script that fits once with the default
display and then reads `e(boot_ci_bca)` for a BCa interval. The matrix is not there, with no error — the read
silently returns nothing. The error message they hit when they mistype `bootdisplay()` repeats the same wrong
promise, so the correction they are given is itself wrong.

**Pathways.** (a) Make the second sentence canonical: reword the opening to "the percentile, bias-corrected
and Wald families are computed and stored on every run; BCa is computed only when requested", and fix 6149 to
say "all four" and name the BCa exception — cheap, text-only, and it makes the `e()` coverage claim honest.
The sweep must cover **five** places, not two: sthlp:504, sthlp:1061 ("All three inference families are
stored on every bootstrap run", which itself undercounts once BCa is displayed), ado:6149, and the source
comments at 6123 and 8783. (b) Make the opening canonical by always running the jackknife — G extra fits on
every clustered bootstrap run, which the author explicitly declined. (c) Add an explicit "what is stored
when" table to the Stored Results bootstrap block — most complete, most text, and the option entry still
reads as self-contradictory in isolation.

**Defensible?** No, as text.

**Verification.** Confirmed; both halves verbatim and the code gates traced. Corrections applied: the stale
claim appears in five places, not the two the finding named, and sthlp:1061 — cited as agreeing with the
second version — itself undercounts. The verifier also spotted an independent stale remnant in the same
neighbourhood: the comment at 6126-6129 describes a `bootdisplay` value "sd" that is no longer legal at 6148.
See CS-11, CS-35, CS-76.

#### CS-116 — `gteffects` documents its interval columns as `lb_95ci`/`ub_95ci` while `level()` claims to govern them and the code names them `lb_<level>ci`

**X.** The `gteffects` option entry lists the returned table's columns with hardcoded 95% names: sthlp:957,
"Columns: ... {it:estimate}, {it:se}, {it:p}, {it:lb_95ci}, {it:ub_95ci}, {it:inwindow} ...".

**Y.** The `level()` entry, rewritten for F054/F044, says `level()` governs the `gteffects` cells and that the
interval columns are named for the level in force, with a documented decimal-point convention: sthlp:489,
"... the {it:gteffects} cells and the Roth slope interval ... The interval columns are labelled {it:lb_#ci} /
{it:ub_#ci} for the level in force (a non-integer level writes its decimal point as {it:p}:
{it:lb_92p5ci})". The mata builder writes exactly that: ado:10801-10802 sets `gtal` from the passed level and
`gtlt = subinstr(strofreal(100-200*gtal,"%9.0g"),".","p",.)`, and 10848 writes the colstripe as
`"lb_"+gtlt+"ci"` / `"ub_"+gtlt+"ci"`. The level is real, not incidental: 1939 passes `wdlvl`, set at 529 as
`c(level)` immediately after the `set level` block, so `level(90)` yields `lb_90ci` and `level(92.5)` yields
`lb_92p5ci`, exactly as the main tables already do.

**Why a user stumbles.** Scripts that post-process `e(gteffects)` by column name — `matrix G = e(gteffects)`
then `colnumb(G, "lb_95ci")` — break the moment the user adds `level(90)`, returning a missing column index
rather than an error. The `gteffects` entry is the natural place to look for that column layout, and it is
the one place that still asserts 95.

**Pathways.** (a) Make `level()` canonical: change the `gteffects` column list to `lb_#ci`/`ub_#ci` with a
pointer to `level()` — a one-line edit, no behaviour change, aligned with how the main results tables are
already described. (b) Make `gteffects` canonical by exempting the cells from `level()` — reverses part of
F054 and leaves the cells' intervals disagreeing with every other interval in the same output. (c) Keep
`lb_95ci` as the illustrative default but annotate it "(lb_#ci at other level() settings)" — smallest diff,
but users skimming a column list read the literal name.

**Defensible?** No.

**Verification.** Confirmed. sthlp:957 is the only hardcoded `lb_95ci` left anywhere in the help; every other
occurrence in the ado is either an internal plot-drawing alias (10556-10660, renamed back before saving) or a
`subinstr` that rewrites it to the level in force (6997, 7509, 7548, 7998, 8777, 8942). One clarification: the
`.025` fallback in `gtal = (args() >= 11 ? ... : .025)` is dead in practice — the sole caller at 1939 always
passes 11 arguments — so the columns read `lb_95ci` at the default only because 95 is the default level, not
via that branch. See CS-2, CS-62.

#### CS-117 — the Roth companion scalars are still documented as a 95% CI although `level()` is now passed into `wooldid_ptp`

**X.** Stored Results describes `e(roth_slope_ll_type)`/`e(roth_slope_ul_type)` as a 95% CI: sthlp:1103,
"Companion scalars: {it:e(roth_slope_{it:type})}, {it:e(roth_slope_se_{it:type})},
{it:e(roth_slope_ll/ul_{it:type})} (the GLS pre-trend slope and 95% CI)".

**Y.** The `level()` entry names "the Roth slope interval" among the things `level()` governs (sthlp:489),
and both `wooldid_ptp` call sites pass the user's level — ado:4595 (analytic leg) and 4596 (bootstrap leg,
`hbaronly`), each ending `level(`wdlvl')`, with the replay path at 418 passing
`level(`=cond(`level' == 0, c(level), `level')')`. The plumbing is real, not nominal: `wdlvl` is `c(level)`
after the F054 `set level` block (529), and inside `wooldid_ptp` the critical value is computed from it —
9294-9295, `if `ptpdfok' local ptpcrit = invttail(`dfr', (100 - `level')/200)` else
`invnormal(1 - (100 - `level')/200)` — and that `ptpcrit` is what produces the scalars ereturned at
9572-9573. So `level(90)` demonstrably moves the stored Roth slope interval while the help calls it a 95% CI.

**Why a user stumbles.** A user running `pretrendpower(80)` with `level(90)` reads a 90% slope interval out of
`e()` while the documentation tells them it is 95%. Roth slope intervals are exactly the kind of number that
gets pasted into a table footnote with its nominal coverage — a wrong coverage label is a reporting error,
not a cosmetic one.

**Pathways.** (a) Make `level()` canonical: change "95% CI" to "{cmd:level()} CI" — one word, zero risk. (b)
Make the Stored Results entry canonical by pinning the Roth slope to 95 regardless of `level()` — reverses
F044 and would make one interval in the output disagree with every other. (c) "level() CI (95% by default)" —
accurate and self-explaining, slightly longer.

**Defensible?** No.

**Verification.** Confirmed. A grep for "95%" over the whole help returns exactly two hits — this line and
sthlp:413, the permute shift-interval width **measurement**, which is not a claim about the option — so this
is the only surviving literal 95% attached to an interval as a statement of behaviour. The one residual risk,
a path reaching `wooldid_ptp` without `level()` and taking its `level(real 95)` syntax default (9279), is
unreachable: all three call sites pass it. That stray default is itself CS-62's third `level` convention.

#### CS-118 — `e(boot_stat)`'s documented value list omits `bootdisplay(bca)`

**X.** Stored Results enumerates `e(boot_stat)` over three display settings: sthlp:1061, "macros
{it:e(boot_method)}, {it:e(boot_stat)} ({it:z} or {it:t} under {cmd:bootdisplay(wald)}, {it:z_wald} or
{it:t_wald} under the percentile and bias-corrected displays, {it:t_pt} under {cmd:boottstat})".

**Y.** `bootdisplay()` is a four-valued option, and the code's rule is "wald or not":
ado:9156, `ereturn local boot_stat = cond(`wdbts', "t_pt", cond("`bdisp'" == "wald", "`wdblab'",
"`wdblab'_wald"))`, so `bca` falls into the `_wald` branch and posts `z_wald`/`t_wald` with nothing in the
help saying so. `bca` is reachable — 6148 accepts it and it is refused only with `boottstat` (6169). The same
omission runs through `e(boot_display)`, whose legal values are never enumerated at all.

**Why a user stumbles.** The stored-results text reads like a complete case list, so a user writing display
code that branches on `e(boot_stat)` has no entry for the BCa run and cannot tell whether its absence means
"behaves like percentile" or "not yet handled". It also plants the idea that BCa is not a rank-based display,
when it is. The omission reads as deliberate exclusion rather than oversight, because the same sentence-block
names `bca` two clauses earlier ("under bootdisplay(bca) also e(boot_ci_bca), e(boot_p_bca), e(boot_accel)").

**Pathways.** (a) Make the code canonical and rewrite the enumeration as "z/t under bootdisplay(wald);
z_wald/t_wald under the percentile, bias-corrected and BCa displays; t_pt under boottstat" — one clause, and
it also settles the reading of the `t_wald` label under BCa. Do the same pass on the two neighbouring stale
enumerations in the same paragraph: "All three inference families are stored on every bootstrap run" and "The
reported standard error is the bootstrap standard deviation under all three", both of which silently exclude
the BCa display. (b) Give BCa its own label and record the family in `e(boot_display)` instead — a code change
and a test re-pin for no user-visible gain. (c) Replace the inline enumerations for `e(boot_stat)` and
`e(boot_display)` with a two-column table of display value → statistic label — best for scanning, most
editing.

**Defensible?** No, though it is a contract hole rather than a wrong number.

**Verification.** Confirmed. Corrections applied: the two neighbouring stale enumerations in the same
paragraph, which pathway (a) should sweep together. See CS-63, CS-115.

#### CS-119 — `e(perm_pconv)` is documented with two values in Stored Results and five in the `permabsp` entry, and `e(perm_absp)` has no `permside`/`permrawp` twin

**X.** The `permabsp` option entry gives the full five-value vocabulary: sthlp:475, "The convention in force
is recorded in {it:e(perm_pconv)} ({it:equaltailed} or {it:absolute}; {it:raw}, {it:left} or {it:right} under
the two options below)."

**Y.** The Permutation-inference Stored Results block, in the same help file, enumerates only two: sthlp:1073,
"{it:e(perm_analyticse)}, {it:e(perm_absp)}; macros {it:e(perm_method)}, {it:e(perm_stat)} ({it:b} or
{it:t}), {it:e(perm_pconv)} ({it:equaltailed} or {it:absolute})". It also lists `e(perm_absp)` among the
scalars with **no gloss**, immediately beside a two-valued `e(perm_pconv)`, so the pair reads as a redundant
binary flag plus binary macro. And the same block still describes `e(perm_p)` as "the equal-tailed
permutation p, or the absolute-value p under {cmd:permabsp}" (sthlp:1065) — also unaware of `permside` and
`permrawp`. The code is five-valued (ado:7648) and the scalar is binary
(`ereturn scalar perm_absp = (`wdpabs' == 1)`, 7649), so `permrawp` (wdpabs 2) and `permside` (3/4) all leave
`e(perm_absp)` at 0; no `e(perm_side)` or `e(perm_rawp)` exists anywhere.

**Why a user stumbles.** A user auditing which p-convention produced an archived result reads the Stored
Results block, sees a binary macro plus a matching binary scalar, and concludes the run was equal-tailed when
it was `permside(right)`. The footer code branches on all five values (10009-10018), so the display is
complete and only the Stored Results block is stale.

**Pathways.** (a) Make the option entry canonical: copy the five-value list into Stored Results, extend the
`e(perm_p)` parenthetical the same way, and re-describe `e(perm_absp)` as "the legacy binary flag, 1 only
under permabsp; read e(perm_pconv) for the convention" — text-only, and it retires the trap. (b) Retire
`e(perm_absp)` and keep only `e(perm_pconv)` — cleanest surface, breaks any script already reading the
scalar. (c) Add `e(perm_side)` and `e(perm_rawp)` scalars so every convention has a flag — symmetric, but it
multiplies redundant state and needs cert pins.

**Defensible?** No — a plain documentation lag behind F064.

**Verification.** Confirmed. Corrections applied: "lists e(perm_absp) as if it were the convention flag" is
stronger than the text supports — the block gives no description at all, and the trap is the adjacency; and
the same block's `e(perm_p)` parenthetical is stale in the same way and should be fixed in the same pass. See
CS-61.

#### CS-120 — the cell collapse's four `e(wdid_fc*)` objects are split three ways: two documented, two undocumented, and a blanket sentence disowning all of them

**X.** The `fastcollapse` and subgroups-screen paragraphs document `e(wdid_fcol)` and `e(wdid_fcguard)` as
things the user can read: sthlp:281, "is recorded in the hidden macro {it:e(wdid_fcol)} (1 when it was
applied, 0 when it was not)"; sthlp:364, "The outcome is recorded in the hidden macro {it:e(wdid_fcguard)}:
empty when nothing was flagged, {it:fallback}, or {it:caution}" — a vocabulary that matches the code exactly
(ado:2100 empty init, 5708 `fallback`, 5712 `caution`, no fourth value); sthlp:318 cites `e(wdid_fcguard)` as
the user's audit handle for detecting that a `subgroups()` run switched to full data.

**Y.** The Replay paragraph in Working With Output says the whole namespace is off-limits: sthlp:1183,
"internally, the display state travels in hidden {it:e(wdid_*)} macros, which are implementation details and
not part of the documented {it:e()} surface". Meanwhile F061 (`de0d5bb`, 2026-09-06) added
`e(wdid_fcncell)` and `e(wdid_fcskip)`, which appear **nowhere** in the help — a grep of the four names over
`wooldid.sthlp` returns exactly three hits, all `fcol`/`fcguard`, while the ado posts all four at 4555-4561.

**Why a user stumbles.** The `fastcollapse` entry explicitly tells reproducibility-minded users to check
`e(wdid_fcguard)` before comparing archived output; the Replay paragraph tells them nothing named `e(wdid_*)`
is safe to depend on. A user cannot tell whether their audit script is sanctioned. And the two objects added
precisely to answer "how many cells, and why was it skipped" — the natural follow-ups to the two that are
documented — are invisible unless the user runs `ereturn list, all`.

**Pathways.** (a) Carve out a documented sub-family: add a short "Cell-collapse diagnostics" synopt to Stored
Results listing all four with their value vocabularies, and narrow the Replay sentence to "the display-state
{it:e(wdid_*)} macros" — makes the promise explicit, costs a paragraph, and commits the author to keeping the
four names stable. This is the better pathway for a further reason: `e(wdid_fcguard)` is read **internally**
to gate permute and bootstrap behaviour (ado:135/139, 305/310), so it is load-bearing across engines, not
merely a diagnostic. (b) Disown all four: strip the references from the `fastcollapse` entry and give users a
`verbose`-output route instead — smallest documented surface, but it removes the only machine-readable way to
check whether an archived run was collapsed, which the `fastcollapse` entry itself argues matters. (c)
Promote the four to non-hidden names (`e(fcollapse)`, `e(fcguard)`) and document them normally — clearest,
but it changes the `e()` surface and needs cert re-pins.

**Defensible?** No.

**Verification.** Confirmed; the grep counts and the `git log -S wdid_fcncell` attribution to `de0d5bb` both
check out. Two qualifications: the Replay sentence's relative clause is non-restrictive, so read strictly it
*could* be scoping the disowning to the display-state members — but the wildcard is written as the whole
family, so a user checking whether `e(wdid_fcguard)` is safe gets no signal that it is exempt; and `fcskip`
is not entirely invisible, since ado:2143 prints its reason on screen ("the default cell collapse does not
apply with `fcskip'"), so what is missing for that one is the machine-readable handle rather than any
disclosure. See CS-70.

#### CS-121 — three p-value conventions live in four help passages, each defended locally, and one split sits inside the bootstrap engine itself

The help half of **CS-1**, which the ado half has now partly outrun.

**X.** `permute()`'s default is the add-one equal-tailed rule, defended at length in the "Which convention,
and the floors" paragraph (sthlp:481) as the exact dual of the shift interval; and the joint, equality and
adversarial tests under **both** engines use the add-one share (sthlp:485 for the fixed-metric form).

**Y.** The bootstrap coefficient p-values use plain order-statistic counts, defended in the `bootdisplay()`
entry (sthlp:504) as what keeps each p dual to its own interval, "measured: 6 disagreements in 35,646
constructed cases, against 36 for the add-one alternative"; `boottstat`'s coefficient p is likewise
no-add-one (sthlp:499). So within a single bootstrap run the coefficient rows and the joint-test rows use
different conventions, and each entry justifies its own choice persuasively without naming the other. The F033
research record states the same thing from the code side and leaves it for the author.

**What changed on the ado side, and what did not.** `b6ef398` rewrote the footer to qualify the resolution
sentence — "The rank-based **coefficient** p-values have resolution 2/n ... The joint, equality and
adversarial p-values are add-one shares over their own listwise-complete replication counts" — so a user
reading the output is now told. The help was not touched: sthlp:504 still says "when it happens a line under
the table and the run's footer say that the 0 means p below 2/n" with no joint-test qualification, and no
help passage cross-references the other conventions.

**Why a user stumbles.** Two p-values printed in the same bootstrap output — an event-study horizon's and the
joint pre-trend test's — are computed by different rules, and a user comparing them at `reps(99)` can see the
joint test refuse to go below .01 while a coefficient p reads exactly 0.

**Pathways.** (a) Unify on add-one everywhere (the F033 record finds it conservative and never
anti-conservative) and rewrite the `bootdisplay()` duality defence — costs the measured duality (6 versus 36
disagreements) and re-pins every bootstrap p cert. (b) Unify on no-add-one for the joint tests — F033 measured
the no-add-one joint p exceeding alpha by 1/(n+1) where listwise n < R, trading a documented conservatism for
a documented anti-conservatism. (c) Carry the footer's new qualification into the help and add one
cross-referencing "which p-value where" block naming coefficient rows, joint/equality/adversarial rows and
each engine. Cheapest, loses nothing measured, and it is the half `b6ef398` did not do. *Not verified (cap)*
for this framing; CS-1 was verified.

#### CS-122 — the bootstrap prints and documents a "p = 0 means p < 2/n" line; `permrawp` can also print exactly 0 and neither the help nor the code offers the analogue

Same defect as **CS-8** (engines lens, partly confirmed) and **CS-80**; their corrections apply, notably that
sthlp:481 does state "0 under {cmd:permrawp}" among the floors, so the gap is on-screen rather than
documentary. **Adds the help-side asymmetry:** sthlp:504 documents both the under-table line and the footer
sentence for the bootstrap, while the `permrawp` entry (sthlp:479) says nothing about a printed 0 or any
resolution line, and a grep for the disclosure string finds only the bootstrap site (ado:9645, keyed on
`e(boot_completed)`). Because `permrawp` exists specifically to reproduce ritest/permute numbers for a
referee, this is the p-value most likely to be quoted verbatim. Pathways: (a) add the same under-table line
and footer sentence to the permutation display when any printed p is 0 — a small code addition and one new
caution to pin; (b) drop the bootstrap line and rely on the footer alone — loses a disclosure just shipped;
(c) documentation only, adding a sentence to the `permrawp` entry — free, but the disclosure stays absent
where it is needed, in the output. *Not verified (cap)* for this framing.

#### CS-123 — `reps()` imposes no minimum while `permreps()` has a hard floor of 20 and a resolution caution, and the sibling entries state opposite policies without cross-referencing

The help half of **CS-20** (confirmed there). **X:** sthlp:487, "No minimum is imposed, so small values can
be used to time a specification - but bootstrap inference needs many replications: hundreds for standard
errors, and a thousand or more before the percentile and bias-corrected intervals are worth reading." **Y:**
sthlp:441, "{cmd:permreps()} must be at least 20, and sets the resolution ... below 39 a 5% rejection is
unreachable under the default convention and a caution says so", backed by ado:5817. **Why it matters:** a
user who learned that wooldid protects them from an under-resolved randomization test assumes the same
protection on the bootstrap side, runs `reps(50)` to iterate, and reports a percentile p of 0.040 whose
resolution is 2/50 = .04 — the floor, not evidence. The bootstrap's resolution facts *are* stated, in the
`bootdisplay()` entry, but as a property rather than a guard, and nothing in `reps()` points at them.
Pathways: (a) add a reps-based caution and mirror the resolution ladder into `reps()` — reverses F030; (b)
drop the permute floor — reverses F036; (c) add a cross-reference sentence to each entry saying the other
engine's policy differs and why (a permutation p has a hard combinatorial floor, a bootstrap p only a
Monte-Carlo one) — cheapest, and the asymmetry is genuinely defensible, it just needs to be visible.
**Fix the 39/40 boundary either way:** the code's caution advises "at least 40" (`cond(absp,"20","40")`)
while its own trigger fires at n ≤ 38 and the help says "below 39", so the two documents disagree about
whether 39 is enough (CS-132). *Not verified (cap)* for this framing.

#### CS-124 — every permutation option carries a `perm*` prefix, but two of the bootstrap options are bare names, and `ties` never says it needs `bootstrap()`

**X.** The permutation surface is uniformly prefixed: `permute`, `permreps`, `permexact`, `permstrata`,
`permprobs`, `permtstat`, `permtci`, `permsave`, `permanalyticse`, `permabsp`, `permside`, `permrawp`.

**Y.** The bootstrap surface is `bootstrap`, `bootalpha`, `bootdf`, `boottstat`, `bootdisplay`, `bootsave` —
plus `reps()` and `ties`, which carry no prefix and sit in the same synopsis group as the `perm*` options
with nothing marking their engine (sthlp:72). `level()`, `dots()` and `rseed()` are genuinely shared by both
engines and are formatted identically to `reps()`/`ties` (CS-28). The `reps()` entry opens "required with
{cmd:bootstrap()}" (sthlp:487); the `ties` entry (sthlp:495) never states an engine requirement at all,
though ado:6033 refuses it.

**Why a user stumbles.** From the synopsis line a user cannot tell which of `reps`/`ties`/`level`/`dots`/
`rseed` belong to which engine. Typing `permute(shuffle) reps(500)` is the natural mistake — `reps` is the
Stata-wide name — and produces an error naming `bootstrap()`; typing `permute(shuffle) ties` produces an
error the `ties` entry gave no warning of.

**Pathways.** (a) Unify on prefixes: add `bootreps()` and `bootties` as the documented names with
`reps()`/`ties` kept as undocumented aliases — clean scheme, but it makes `reps()` (the Stata-wide
convention, which the help elsewhere leans on when comparing to the `bootstrap` prefix) a second-class name,
and every example and test needs updating. (b) Unify by dropping prefixes where unambiguous — worse; the
`permsave`/`bootsave` and `permtstat`/`boottstat` twins genuinely need them. (c) Split the synopsis "Standard
Errors" run into labelled sub-runs (analytic / cell collapse / permutation / bootstrap / shared by both
engines) and add "requires {cmd:bootstrap()}" to the `ties` entry as `reps()` already has — no code or name
change, resolves the actual confusion, cheapest of the three. Defensible. *Not verified (cap).* See CS-17,
CS-28.

#### CS-125 — engine-parallel `e()` objects diverge in name and in whether they are documented

**X.** The permutation Stored Results block documents the Monte-Carlo precision surface fully: sthlp:1066,
"{it:e(perm_p_se)} and {it:e(perm_p_ci)} (the Monte-Carlo standard error of the p-value's share and its
{cmd:level()} Clopper-Pearson interval)", plus `e(perm_n)` and `e(perm_level)`.

**Y.** The bootstrap block documents `e(boot_navail)` and `e(boot_level)` and never mentions
`e(boot_ptails)`, which F064 posts as a hidden matrix (ado:9133; grep count in the help: 0). And the level in
force is recorded three times under three names — `e(level)` (4554), `e(perm_level)` (7642),
`e(boot_level)` (9155) — documented separately at sthlp:489, 1061 and 1072 with no cross-reference and, per
CS-2, capable of disagreeing.

**Why a user stumbles.** A user who learned the permutation surface expects a bootstrap twin for each object
and goes looking for `e(boot_p_se)`/`e(boot_p_ci)`; the tail record that would let them compute it exists but
is undocumented, so they conclude the bootstrap has no Monte-Carlo precision surface. Separately, a script
reading `e(level)` on a bootstrap run gets a value and one reading `e(boot_level)` also gets a value, with
nothing saying whether they can ever differ — so which one is authoritative for a footnote is a guess.

**Pathways.** (a) Unify names: post `e(level)` alone and retire `e(boot_level)`/`e(perm_level)`, or make the
engine-specific ones documented aliases — cleanest, but it breaks existing scripts and the footer code at
10097-10103 reads `e(boot_level)` directly. (b) Unify coverage: document `e(boot_ptails)` in the bootstrap
block and add a sentence to each engine block naming its twin in the other ("e(boot_navail) is the bootstrap
twin of e(perm_n)") — no code change, and it turns the two blocks into a readable pair. (c) Add a short
"engine twins" table to Working With Output mapping each `perm_*` object to its `boot_*` counterpart and
marking the ones with no twin (step-down, `ties`, jackknife acceleration on one side; shift intervals, null
band, comparator on the other) — most useful, most to maintain as the engines evolve. Defensible in part.
*Not verified (cap).* See CS-33, CS-56, CS-71.

#### CS-126 — Working With Output declares `lincom` and `test` incompatible while `bootdf()` promises that a hand `test`/`lincom` reproduces the reported p-values

The help half of **CS-4** and **CS-77**, both confirmed. **X:** sthlp:1186, "{ul:Incompatible
Post-Estimation Commands:} {it:lincom}, {it:test}, and all other post-estimation commands that rely on
receiving a valid variance-covariance matrix ... these are for display purposes only." **Y:** sthlp:503 twice
promises the opposite for bootstrap runs, in bold, as the reason `e(df_r)` is aligned — "so a {cmd:test} or
{cmd:lincom} run by hand afterwards reproduces the reported p-values rather than diverging from them" and,
under the default percentile display, "a {cmd:test} or {cmd:lincom} on {it:e(V)} reproduces the Wald family
rather than the printed p". And `fastnose` adds a third position: "{it:e(V)} is not posted ... so {it:test},
{it:lincom}, {it:esttab} standard errors, and anything else needing a variance-covariance matrix will refuse"
(sthlp:267). **Why it matters:** three sections give three answers to "can I run `test` after wooldid?" —
never, yes-and-it-matches, and it-will-refuse — and the one that says "never" is the section titled for that
question, so it is where a user looks. The genuinely important nuance (single-coefficient `lincom` is fine on
the diagonal; cross-coefficient `test` is not, except on the reposted `e(b_es_*)`/`e(V_es_*)`) is stated in
the same paragraph but after the blanket refusal. Pathways: (a) make Working With Output canonical by
softening it into a three-case rule and noting the bootstrap addendum and the `fastnose` case — text only,
and it is the section users read; (b) delete the `bootdf()` promise — loses a real guarantee the F053/`bootdf`
work was built to provide; (c) cross-reference all three — least editing, but the section titled "Incompatible
Post-Estimation Commands" still opens with a claim the other two contradict. *Not verified (cap)* for this
framing.

#### CS-127 — the `permute()` and `bootstrap()` incompatibility lists are structured differently and name different options, so absences cannot be read

**X.** `permute()` ends with "Not compatible with impvalues, unconditionalse, contreatwithin,
contreatcontrols(untreatedZero) [with a reason], oldsyntax, suppresssgprimaryeffects [with a reason],
pretrendpower() [with a reason], or bootstrap()" (sthlp:441) — some entries carry parenthetical reasons, some
do not.

**Y.** `bootstrap()` ends with "Incompatible with permute(), impvalues, fastnose, unconditionalse,
customvce(), suppresssgprimaryeffects, oldsyntax, and multi-way cluster()" (sthlp:485) — no reasons at all,
naming three things permute's list does not (`fastnose`, `customvce()`, multi-way `cluster()`) while omitting
three that permute names (`contreatwithin`, `contreatcontrols(untreatedZero)`, `pretrendpower()`).
`pretrendpower()` is separately flagged in the bootstrap entry as now supported.

**The asymmetries have four different causes and the help distinguishes none of them.** `fastnose` is
permute's **default** (sthlp:267, 441) rather than a refusal; `customvce()` is **unreachable** under
`permute()` by construction, since permute requires `cluster()` and ado:4758 refuses `customvce()` with
`cluster()`; multi-way `cluster()` is likewise excluded by permute's single-cluster requirement; and
`contreatwithin` is **refused explicitly** at 5750.

**Why a user stumbles.** A user reading the two lists side by side sees `fastnose` refused by one engine and
unmentioned by the other, and cannot tell that under permute it is not merely allowed but the default. They
see `customvce()` refused by bootstrap and unmentioned by permute and reasonably infer permute supports it,
then hit the unrelated `customvce`/`cluster` error with no idea why. An absence in one list means four
different things — allowed, default, structurally unreachable, or simply unwritten.

**Pathways.** (a) Unify the format: give both engines the same list shape with a one-clause reason on every
entry, and mark structurally-unreachable combinations as such rather than omitting them — text only, and it
makes the two lists diffable, which is what a user actually does with them. (b) Unify membership by making
each engine name every option in the other's list with its own verdict — most complete, two long lists to keep
in sync. (c) Move both into one "Engine compatibility" table with a row per option and a column per engine
(allowed / default / refused / unreachable) — highest scanning value, and the natural home for the
fastnose-is-the-default fact; a new table to maintain. Defensible. *Not verified (cap).* See CS-48, CS-94,
CS-107.

#### CS-128 — `gteffects` is called "the one inference surface that survives fast mode" and then has its inference columns blanked under `permute()`, whose default *is* fast mode

**X.** sthlp:966: "Fully available under {it:fastnose} (the regression's own VCE is intact even though the
{it:margins} stage skips variances), making {it:gteffects} the one inference surface that survives fast
mode."

**Y.** sthlp:967, in the same entry: "Under {cmd:permute()} the {it:gteffects} cells are keyed by cohort
onset, which the replications reassign, so they have no permutation distribution and their inference columns
are blank (the points stay)" — and sthlp:441 says "Without {cmd:permtstat} the observed run and the
replications all run in fast (no-SE) mode". The blanking is recorded in `e(perm_gt_blank)` (sthlp:1072).

**Why a user stumbles.** A user reads the fast-mode sentence, adds `gteffects` to a `permute()` run expecting
to keep some inference, and gets an all-blank inference block with only a caution. The two sentences are true
of different fast paths — plain `fastnose` keeps the analytic VCE, permute's fast path additionally has no
permutation distribution for onset-keyed cells — but the entry states the first as an unconditional slogan.

**Pathways.** (a) Qualify the slogan: "...the one inference surface that survives plain {it:fastnose} (under
{cmd:permute()} the cells are blanked for a different reason - see below)" — one clause, keeps the useful
claim. (b) Drop the slogan and state the three cases plainly (fastnose: analytic cell SEs kept; bootstrap:
bootstrap cell SEs; permute: blank) — clearer, slightly flatter. (c) Harvest a permutation distribution for
the cells so the slogan becomes true — a real engineering change with no obvious meaning, since a replication
reassigns exactly the key the cells are indexed on. *Not verified (cap).* See CS-13, CS-50, CS-102.

#### CS-129 — three sibling "fewer than half of the replications" cautions exist in the code; only the `boottstat` one is promised in the help

**X, promised.** sthlp:501, the `boottstat` "Costs and limits" paragraph: "a caution names any coefficient
that used fewer than half of the completed replications", pointing at `e(boot_tstud)` row `n` (implemented at
ado:9169).

**Y, printed but unpromised.** ado:7661, "Caution: `wdplow' reported coefficient(s) were present in fewer
than half of the completed replications (a reassigned onset can make a horizon or a subgroup level
inestimable in a replication); their permutation p-values rest on the counts in e(perm_n)."; and ado:9231,
"Caution: the following estimates were available in fewer than half of the completed replications; their
bootstrap standard errors rest on the listed counts". Neither the `permute()` nor the `bootstrap()`/`reps()`
entries mention that a caution is printed — they mention only the counts (`e(perm_n)`, `e(boot_navail)`). The
`subgroups()`+`permute(shuffle)` caution at 5954 is a second such case: the help's Subgroups paragraph
(sthlp:441) describes the behaviour without saying a caution prints.

**The gap is one-directional.** Five promised cautions were spot-checked and all exist:
`cfxlattice`-without-`makecfxplots` (5293), `esrelativeto`-without-`esfixedbaseperiod` (5357),
fewer-than-6-treated-clusters (6279), `permtstat` at `cluster(i)` (5902), and `permexact` above 1,000
assignments / refused above 10,000 (10973-10980).

**Why a user stumbles.** A user who has read the `boottstat` entry expects wooldid to warn about thin
replication coverage and treats the absence of such a warning elsewhere as a clean bill of health. When the
permutation caution *does* fire on a staggered design they cannot find it documented anywhere, so they cannot
tell whether it is routine — it often is; a horizon simply drops out of some draws — or a sign of a broken
run.

**Pathways.** (a) Unify by documenting: add one sentence to `permute()` and to `bootstrap()`/`reps()`
promising the same caution the `boottstat` entry promises, naming the count object each points at — no code
change, and it makes the caution readable when it fires. (b) Unify by suppressing: keep only the `boottstat`
caution, since the other two are frequently benign — loses genuine disclosure on runs where coverage really
is thin. (c) Add a general "Cautions wooldid prints" section listing the whole family with one line each —
the most durable fix for a package with ~100 caution sites, and a substantial new section to write and keep
current. *Not verified (cap).*

#### CS-130 — the same option's argument is written differently in the syntax diagram and in its own entry, and the newest options lost their abbreviation markers

**X.** The "Options by Purpose" synopsis writes `contreatelasticitytype(dydx, eydx, dyex, or eyex)`,
`regtol(+int)`, `lattice(variable[+int] OR latticetype latticenum)`, `subgroups(variable[>=0 int])`
(sthlp:72, 82).

**Y.** The option entries write `contreatelasticitytype(dydx/eyex/eydx/dyex)` (sthlp:734) — different
separator and different order — `regtol` with no argument at all (sthlp:936), `lattice(latticevar OR
latticetype latticenum)` (sthlp:745), `subgroups()` (sthlp:694). Separately, `permside` and `permrawp` are
the only options in the standard-errors group written without a minimum-abbreviation colon in either place,
though ado:506 declares `PERMSide(string) PERMRAWp` — so `perms` and `permraw` work but are undocumented, and
`perms` resolves to `permside` while `permsa` is `permsave`, a near-collision the help gives the user no way
to see (CS-30). `bootdisplay(spec)`'s four legal values are never enumerated in the synopsis, unlike
`permside(left|right)` and `contreatelasticitytype` right beside it.

**Why a user stumbles.** The synopsis is the abbreviation reference: users learn from it that `permts` is
enough for `permtstat`, then guess wrongly at `permside` because it has no marker, or type `perms` meaning
`permsave` and silently get `permside`. The differing argument spellings make it unclear whether
`contreatelasticitytype` takes a comma-separated list (the synopsis reads like one) or a single value (the
entry's slashes).

**Pathways.** (a) Make the option entries canonical and regenerate the synopsis from them mechanically —
guarantees they cannot drift again, and is the standard Stata help discipline; costs one editing pass. (b)
Make the synopsis canonical and copy its spelling into each entry — same effect, but the bracket notation
reads badly as an entry heading. (c) Fix only the two that can mislead: add the abbreviation markers to
`permside`/`permrawp` and give `regtol` its argument. *Not verified (cap).*

#### CS-131 — the help names the same concepts several ways across sections

**X.** `bootstrap(mode)` is "mode" in the syntax diagram and entry heading, but the same paragraph calls it
"the Bayesian **flavor**" (sthlp:485), and sthlp:345 has "with every **flavour** of {it:contreatcontrols()}"
— a third sense and a spelling change. "Family" does triple duty: the display families
(percentile/bc/bca/wald), the joint-test construction in `e(boot_jt_family)` (studentized/fixedmetric), and
the estimand families in the step-down paragraph ("for each family and subgroup", sthlp:483).

**Y.** The event-study objects appear as `e(b_es_)` (fastnose entry, sthlp:267), `e(b_es_*)` (permute entries,
sthlp:441) and `e(b_es_{type})` (Stored Results, sthlp:1043). Fast mode appears as "Fast exploratory mode",
"fast mode", "fast (no-SE) mode", "this fast path", and — a different thing entirely — "G extra **fast fits**"
for the BCa jackknife (sthlp:504). The observed run is "the rep-0 comparator" and `e(perm_comparator_b)` on
the permutation side but "the observed estimates" on the bootstrap side (sthlp:506 vs 441), with no comparator
object documented there.

**Why a user stumbles.** Searching the help is how users navigate 1,375 lines of it. Searching "family"
returns three unrelated concepts; searching `e(b_es_att)` misses the two entries that most matter (`fastnose`
and `permute` both withhold it). A user grepping their own log for "fast mode" will not find the BCa
jackknife's "fast fits". And the comparator asymmetry hides a real fact recorded in the F018 research note —
the bootstrap stores no rep-0 cell comparator at all.

**Pathways.** (a) Fix the vocabulary: reserve "mode" for `bootstrap(standard|bayes)`, "family" for the display
families only (rename the `e(boot_jt_family)` prose to "construction" and the step-down "family" to "estimand
block"), settle on `e(b_es_{type})` everywhere with a note that `*` denotes the subgroup suffixes, and pick
one name for fast mode — text-only except the `e(boot_jt_family)` prose, and the highest yield per line
changed. (b) Add a short glossary near the top of Working With Output defining mode, family, fast mode,
comparator, shift interval, null band and cell collapse, and leave the prose alone — cheapest, but the
ambiguous words still appear unlabelled at their use sites. (c) Fix only the `e(b_es_)` spellings, which are
the ones that break a text search for a real object name — smallest defensible subset. Defensible in part.
*Not verified (cap).*

#### CS-132 — three numbers stated in the help disagree with the code that prints them or with the record they came from

**X.** The `permute()` entry sets the unreachable-5% threshold at "below 39" replications and gives the
resolution ladder as 50 → .039, 100 → .020, 199 → .010, 999 → .002 — all exactly 2/(n+1) (sthlp:441). The
help's collapse-under-`impvalues` note quotes a 1e-12 agreement (sthlp:340).

**Y.** The caution the code prints advises "at least 40" for a 5% test — one more than its own trigger
requires, since it fires on `wdp_floor > .05`, i.e. n ≤ 38, matching "below 39" — and rounds the ladder to
"roughly 100 for .02, 200 for .01 and 1000 for .002" (ado:5817). And the internal cert record for the same
`impvalues` measurement reports 7e-15 crossed and exact under nesting
(`review_process/qr_ship_2026-09-06.md`, F062 row), three orders tighter than the help's figure.

**Why a user stumbles.** A user who runs `permreps(39)` after reading the help expects no caution and gets
none — correct — but at `permreps(38)` the caution tells them to use 40, so the two documents disagree about
whether 39 is enough. The 1e-12 figure is not wrong (it is a bound, not the measurement) but it is three
orders looser than what was measured, which understates the case for the change the author just shipped and
makes the help's precision claims look inconsistent with its neighbours, which quote measured quanta rather
than bounds.

**Pathways.** (a) Make the help canonical: change the caution's advice to "at least 39" and its ladder to the
exact 2/(n+1) values, and tighten the `impvalues` figure to the measured 7e-15 (with "bounded by 1e-12" if
the bound is wanted) — text-only on both sides, small cert-text re-pin for the caution. (b) Make the code
canonical: change the help to "below 40" and "roughly 200/1000" — loses the exactness the ladder was written
to provide. (c) Add "roughly" to the help's ladder — hides the discrepancy rather than resolving it, and the
39/40 boundary is a yes/no question a user can hit. *Not verified (cap).* See CS-123.

## 3. Harmonization roadmap

Eight packages a maintainer could ship one at a time, each self-contained and each resolving a group of the
items above. Every one of CS-1 to CS-132 appears in exactly one package. Breakage counts are
`grep -rlF` file counts over `test/` at `2d59e70` and are stated per package rather than per item; they
include the build copy at `test/suite/out/build_new/wooldid.ado`, which is regenerated rather than edited.
Effort is S (an afternoon), M (a day or two with cert pins), L (a multi-commit project with a battery
re-bless).

### P1 — Resolve `level()` once, in the shell

**Resolves** CS-2, CS-36, CS-46, CS-62, CS-116, CS-117.

**User-visible change.** Under `set level 90` with no `level()` option, `bootstrap()` and `permute()`
intervals become 90% intervals instead of 95% intervals printed under 90-named columns; `e(level)`,
`e(boot_level)` and `e(perm_level)` agree on every run; a bad `level()` on a plain analytic call gets
wooldid's message instead of Stata's; and the `gteffects` and Roth-slope help text stops asserting 95.

**How.** One line in the shell before any engine call —
`local level = cond(`level' == 0, c(level), `level')` — which also removes the six subprogram defaults
(6369, 6454, 7417, 7681, 8033, 8639) and the stray `level(real 95)` in `wooldid_ptp` (9279). Then move the
`level()`/`dots()` range checks ahead of `bootctl`'s `exit 0`, and edit sthlp:957 and sthlp:1103.

**Breakage and migration.** Explicit `level()` behaves exactly as today, so nothing in the current battery
moves — no config runs at a non-95 session level, which is why this has never surfaced. Needs one new pin in
`test_level.do` (the L1-L9 file from a969126) covering `set level 90` plus each engine. Watch the worker's
`set level` / `WOOLDID_LEVEL0` restore pair at 68-70 and 525-528: resolving in the shell makes `level != 0`
always true, so the repmode branch must be checked.

**Effort: S.** Ship first — it is a help/code contradiction in a sentence shipped the same day (sthlp:489),
and the wrong number is currently mislabelled rather than merely undocumented.

### P2 — One p-value convention map, and a floor statement wherever a p can print 0

**Resolves** CS-1, CS-8, CS-9, CS-11, CS-20, CS-21, CS-44, CS-80, CS-81, CS-89, CS-121, CS-122, CS-123,
CS-132.

**User-visible change.** Every printed p carries a recoverable convention; a p of exactly 0 is explained
wherever it can occur, not only under the bootstrap; the low-replication caution names the floor the run
actually has; `reps()` and `permreps()` accept each other's spelling; and one help paragraph states the
floors, resolutions and requirement regimes of both engines side by side.

**How.** Text plus two small code changes. (i) Carry `b6ef398`'s new footer qualification into
sthlp:504 and add the "which p-value where" block (CS-121). (ii) Generalise `wooldid_pfloor` to key on "a
rank-based p is in force" rather than on `e(boot_completed)`, computing 2/completed for the bootstrap and
1/completed for `permrawp`, and finding the column as `p_perm` as well as `p`. (iii) Thread `rawp`/`side`
into the `cond()` at 5816 so the F036 caution stops telling `permrawp` and `permside` users a floor they do
not have. (iv) Accept `reps()` under `permute()` and `permreps()` under `bootstrap()` as synonyms. (v) Fix
the 39/40 boundary and the resolution ladder in the caution text.

**Breakage and migration.** No printed number moves unless the author also takes the add-one/no-add-one
decision the F033 record left open — which this package deliberately does **not** force. The caution text
change re-pins Q19; the new `pfloor` branch needs a permute pin. `p_perm` appears in 5 test files.

**Effort: S–M.** The one open question the author must answer is whether to keep three formulas and label
them (cheap) or unify (moves every pinned bootstrap p).

### P3 — One inference-disclosure contract: what `e(V)` is, what is absent, and what survives replay

**Resolves** CS-4, CS-5, CS-16, CS-42, CS-59, CS-65, CS-75, CS-77, CS-82, CS-84, CS-85, CS-90, CS-91, CS-93,
CS-99, CS-102, CS-126, CS-129.

**User-visible change.** Every run — analytic included — says in one line what `e(V)` holds and what
`test`/`lincom` will do with it. Every result the package declines to compute is an `e()` flag restated in
the footer, under both engines, so a replayed or archived run carries its own disclosure. No test line prints
a bare `.`; no footer sentence describes a feature the run did not use; and the joint/equality/adversarial
lines name the engine that produced them.

**How.** Adopt the rule the bootstrap already states in its own comment at 6573 — *cautions are for
conditions, flags are for absences* — and audit all three paths against it. Concretely: post `e(vsource)`
and `e(perm_v_posted)`/`e(perm_v_source)`; move `ereturn hidden scalar df_r_analytic` above the `bootdf()`
if/else (a one-line fix for a promise the help makes twice and keeps on one branch); post
`e(perm_jt_absent)`/`e(perm_sgjt_absent)`/`e(perm_adv_absent)` and choose one polarity against
`e(boot_adv_absent)`; add the `!= "."` guard at 9806/9948 (or, better, stop posting the analytic adversarial
p under an engine, after which the existing gate suffices); print `e(perm_failed)`; gate the bootstrap
contrasts sentence on `wordcount(e(wdid_sglist)) > 1`; make the adversarial wording branch on the engine
rather than on `e(perm_method)`; give the analytic path a short footer; move the bootstrap's standing
disclosures into `part(boot)`; and rewrite the Incompatible Post-Estimation paragraph per engine.

**Breakage and migration.** This is the package that moves the most output. Every analytic cert re-pins if
the analytic footer is added (that half can ship separately); the bootstrap and permute footer baselines move;
the two stray blank lines at shell 47/49 must go for the "byte-for-byte" replay claim to hold. `se_analytic`
appears in 9 test files and `perm_kernel` in 9.

**Effort: M**, or L if the analytic footer and the replay migration ship together. The `df_r_analytic`
one-liner and the adversarial dot guard are worth landing on their own first.

### P4 — Parallel stored results: one name, one shape, one visibility rule

**Resolves** CS-12, CS-17, CS-33, CS-54, CS-55, CS-56, CS-57, CS-58, CS-60, CS-61, CS-63, CS-64, CS-68,
CS-69, CS-70, CS-71, CS-72, CS-73, CS-83, CS-100, CS-101, CS-119, CS-125.

**User-visible change.** A script can move between the engines by substituting one token. The family token in
a hidden table name is derivable from a documented macro; `lb` means one thing; `e(permdist)` is either the
draw matrix or renamed; the joint-test p is a scalar at both scopes; `e(N_clust)` exists; and one "engine
twins" table maps every `perm_*` object to its `boot_*` counterpart, marking the ones with no twin.

**How, in shipping order.** (i) Free additions that break nothing: post `e(N_clust)` beside `e(Nclus)`,
`e(boot_family_tag)`, `e(boot_sd)`, `e(boot_p_wald)`, `e(bootdist_b)`, `e(perm_display)`, and
`e(boot_seed)`/`e(perm_seed)` unconditionally. (ii) Documentation-only: the engine-twins table; the five
values of `e(perm_pconv)` and the gloss on `e(perm_absp)`; the four `e(wdid_fc*)` names and the narrowed
Replay sentence; `e(boot_ptails)`. (iii) Renames with aliases: `e(perm_lb_<est>)` → `e(perm_bandlo_<est>)`
(the one that actually misleads); `e(boot_ci_normal)` → `e(boot_ci_wald)`; `e(permdist)` →
`e(permdist_est)`; `e(perm_stat)` → `e(perm_permuted)` with `e(perm_stat)` re-used for the column label.
(iv) Type and visibility: post the subgroup joint p as a scalar at all three sites **and** fix the display
block at 9928-9942 that dereferences it as a macro; adopt the hidden/visible rule and unhide
`e(boot_ptails)`, `e(perm_kernel)`, `e(boot_ase_obs)`. (v) The step-down family keys, which are cheap to move
**now** — 7 files and 2 help lines — and expensive later.

**Breakage and migration.** Sizing: `boot_pctl_mainresults` 6 files, `perm_kernel` 9, `boot_ci_normal` 5,
`e(joint_horiz_s` 3, `e(permdist)` 3, `perm_p_lo` 7, `boot_ptails` 2, `e(Nclus)` 1, `boot_navail` 1. Steps
(i) and (ii) break nothing; (iii) needs a one-release alias policy and a note in the help; (iv) moves any
baseline that snapshots `ereturn list`. Note that renaming `boot_pctl_*` to `boot_percentile_*` would tighten
the 32-character cfx budget from 22 characters to 16 (CS-54, CS-72) — prefer `e(boot_family_tag)`.

**Effort: L**, but it decomposes cleanly: (i) and (ii) alone are S and close about half the items.

### P5 — One column-stripe vocabulary across engines and surfaces

**Resolves** CS-6, CS-7, CS-13, CS-18, CS-26, CS-66, CS-67, CS-74, CS-79, CS-88.

**User-visible change.** Either every column that stopped being the analytic quantity says so — including
the bootstrap's `p` and `se`, and the `gteffects`/`e(cfx*)` cells under both engines — or none does and the
provenance lives in one documented place. Either way, one help table gives the stripe per engine, per display
family and per surface, and no column header is permanently blank.

**How.** The decision is (a) mark everything or (c) publish the map. If (a): `p` → `p_pctl`/`p_bc`/`p_bca`/
`p_wald`, `se` → `se_boot` on the bootstrap tables; extend the F053 relabelling to the cell matrices
(`t` → `t_wald`/`t_pt`/`t_perm`, `pvalue` → `p_perm`); add the level suffix to the cfx bounds; and drop or
rename permute's permanently blank contrast `t_perm`. If (c): leave the stripes and write the table, but
still fix the blank contrast column and add the level suffix, since those are defects under either rule.

**Breakage and migration.** The heaviest name surface in the package: `lb_95ci` appears in **47** files under
`test/`, `se_analytic` in 9, `p_perm` in 5. `wooldid_pfloor`'s `colnumb(pm,"p")` at 9641 falls back to
positional column 4, so it degrades rather than errors. The cfx colnames are pinned in the cfx cert and in
`saveplots` datasets users may already hold, which is the strongest argument for (c) over a full rename.

**Effort: M** for (c) plus the two defects; **L** for (a).

### P6 — One refuse / caution / skip policy, with `safety(off)` as the single documented downgrade

**Resolves** CS-24, CS-25, CS-34, CS-37, CS-38, CS-40, CS-41, CS-43, CS-47, CS-49, CS-50, CS-51, CS-78,
CS-98, CS-103, CS-105, CS-111, CS-112.

**User-visible change.** A published rule for what an inert or orphaned option does; no accepted option is a
silent no-op; `safety(off)` never silently changes a number or an estimand window; `contreatwithin` announces
what it turned off and stops contradicting itself about plots; and a diagnostic that cancels itself says so.

**How.** (i) Publish the tier rule and move the silent tier into the caution tier, so tier 3 disappears. (ii)
Publish the `safety()` table — screen, default behaviour, behaviour under `safety(off)`, and which macro the
site gates on — and print the disabled list at run time under `safety(off)`. (iii) The three targeted fixes
that need no policy decision: the one-token change at 1446 (`safety` → `safetyinput`) so the maxvar preflight
survives `contreatwithin`; gate the 4912-4920 cfx-plot caution on the post-rewrite value so the run stops
printing "will not produce these plots" and then producing them; and print the `espostlength` raise
regardless of `safety()`, since `safety(off)` should not hide a change to the estimand. (iv) Screen
`adversarial`'s prerequisites in `wooldid_validate` the way `wooldid_ptpcheck` does, or at minimum give
ado:1042 a caution. (v) One vocabulary — refused / reported without inference / fully supported — used
identically in messages, footers and help, which also fixes the self-contradictory 8320 string. (vi) One
message contract in the process checklist, applied opportunistically to the impvalues-era refusals.

**Breakage and migration.** (iii) and (iv) are the behaviour changes: re-arming the maxvar check turns some
`contreatwithin` runs into an early r(907)-preventing error, which is the point; screening `adversarial`
breaks scripts carrying a decorative `adversarial` on runs where it never did anything. Note the trap in
pathway (a) of CS-37: any `contreatwithin` exemption list **must** keep the multicollinearity relaxation,
whose own error text prescribes `contreatwithin` as the remedy. Battery re-bless if (iii) lands.

**Effort: M.** The three targeted fixes in (iii) are each one line and could ship immediately.

### P7 — Engine parity: gates, partners, defaults, and the surfaces only one engine has

**Resolves** CS-3, CS-10, CS-14, CS-15, CS-19, CS-22, CS-23, CS-27, CS-32, CS-35, CS-39, CS-45, CS-48,
CS-52, CS-53, CS-86, CS-87, CS-92, CS-94, CS-95, CS-96, CS-97, CS-104, CS-106, CS-107, CS-108, CS-109,
CS-110, CS-113, CS-127, CS-128.

**User-visible change.** The two engines' compatibility lists become diffable and say why; an option is not
refused by one engine and silently inert under the other; the same failure fraction draws the same reaction;
`cluster()` says which unit it is naming on each run; the cell and plot surfaces behave the same way or say
why they do not; and the features that exist on only one engine are labelled as such rather than looking like
omissions.

**How, cheapest first.** (i) Text: one "Engine compatibility" table with a row per option and a column per
engine (allowed / default / refused / unreachable), which subsumes both incompatibility lists; add
`bootstrap()` and `pretrendpower()` to the `fastnose` entry's list; reword the `permute()`+`pretrendpower()`
refusal to say the permutation covariance is the wrong *kind*, not absent — the package already has that
vocabulary at 9098/9360; note that the bootstrap's few-cluster referral leads to an engine with its own
exchangeability limit; say that the step-down and the `pretrendpower()` support are engine-specific and why.
(ii) One-line gates and cautions: move the `fastnose` refusal block after permute's own assignment (or exempt
permute runs), and close the `fastnose`+`permanalyticse` hole; widen the 5902 exchangeability caution to the
default statistic and extend the floor caution to `permexact`; refuse `customvce()` under permute's fast-mode
default; add a treated-cluster caution to `permctl`; harmonize the failure thresholds; pre-check
`saveplots()` names or wrap the graph saves in `cap`; move the worker's `espre + espl > 0` guard into both
shell `makeplots` blocks; fix the `bootdf()` error text. (iii) Behaviour, needing rulings: extend the
fixed-metric rank test to the adversarial statistic (noting the null changes); port the drift list into
`wooldid_bootcellsx`, or document the divergence and the accidental `latticeignoreweights` exception;
disclose at the point of ame refill; decide whether the histogram is bootstrapped or blanked; build the
bootstrap step-down; add `bootcluster()`/`permunit()`.

**Breakage and migration.** (i) breaks nothing. (ii) turns a handful of currently-silent combinations into
cautions or errors and needs new gate pins. (iii) is where the numbers move and where the battery re-blesses.
Verify CS-87 (`makeplots` with no horizons under an engine) **first** — it is the one item whose severity
depends on an unrun test, and if it errors it is a run-losing bug rather than a cosmetic divergence.

**Effort: L**, decomposing to S for (i), M for (ii).

### P8 — Documentation sync: the 2026-09-06 stragglers and the option surface

**Resolves** CS-28, CS-29, CS-30, CS-31, CS-76, CS-114, CS-115, CS-118, CS-120, CS-124, CS-130,
CS-131.

**User-visible change.** No help passage contradicts the code or another help passage; the synopsis is a
reliable abbreviation reference; and the shared inference companions are findable from the permutation entry.

**How.** (i) The four stale claims from the 2026-09-06 wave: Example 13's parenthetical (CS-114); the
"every family is computed regardless" sentence in **five** places — sthlp:504, sthlp:1061, the error text at
ado:6149 with its "all three", and the comments at 6123 and 8783 (CS-76, CS-115); the `e(boot_stat)` value
list and its two neighbouring enumerations (CS-118); the four `e(wdid_fc*)` names and the Replay sentence
(CS-120). (ii) The option surface: split the synopsis "Standard Errors" run into labelled sub-runs so
`level()`/`dots()`/`rseed()` are visibly shared and `reps()`/`ties` are visibly bootstrap-only; add "requires
`bootstrap()`" to the `ties` entry; add the abbreviation markers to `permside`/`permrawp`; give `regtol` its
argument; regenerate the synopsis from the entries so the two cannot drift again. (iii) Vocabulary: reserve
"mode" and "family" for one meaning each, settle on `e(b_es_{type})`, and pick one name for fast mode. (iv)
The `bootdf()` default message and the `saveplots()` timing note.

**Breakage and migration.** Text only, except the `bootdf()` message string and any test pinning it. The one
judgment call is whether to also fix the abbreviation collisions themselves (`fast` → `fastnose`,
`permt` → `permtci`), which is a breaking change and belongs in its own commit if taken — CS-30 pathway (a).

**Effort: S.**

## 4. Already consistent by design

These asymmetries came out of the sweep as findings and are, on the evidence, inherent to the methods rather
than to the implementation. They are listed so they are not re-raised, and so that any move toward
"making the engines look alike" does not trade them away. Each names the item it belongs to and the part of
that item that is **not** defensible, since in most cases the behaviour is sound and only the labelling or
the documentation is at issue.

1. **The bootstrap's observed run consumes an analytic standard error and permute's does not** (CS-3, CS-23,
   CS-94). `bootstrap()` refuses `fastnose` because the replications have nothing to attach inference to;
   `permute()` defaults to it because randomization inference does not consume the observed analytic SE. Two
   different states wearing one name. *Not defensible:* the ordering trap that makes typing the word break a
   legal run, and the silence about which `e()` matrices go away.

2. **`boottstat`'s poisson and bayes refusals have no `permtstat` analogue** (CS-22, CS-107). Bootstrap
   replications carry duplication or Dirichlet weights, so a replication's `se_r` is a different functional of
   the data than `se_obs`; permutation replications relabel onset without reweighting rows, so that failure
   mode cannot arise. Principled, not an oversight. *Not defensible:* that the two engines' incompatibility
   lists do not say so.

3. **`pretrendpower()` under `bootstrap()` and not under `permute()`** (CS-19, CS-92, CS-95). Roth (2022)
   consumes a **sampling** covariance. The bootstrap supplies one; the permutation covariance is a
   null-relabelling metric, which is why the joint and adversarial tests can use it and Roth's power
   calculation cannot without its own research. *Not defensible:* the refusal's stated reason, that
   randomization inference produces no covariance, which the same output contradicts three lines later.

4. **`gteffects` cells are blanked under `permute()`** (CS-13, CS-102, CS-128). The cells are keyed by cohort
   onset and the replications reassign exactly that key, so there is no permutation distribution to report.
   Blanking rather than filling from the analytic run is the right move. *Not defensible:* the "one inference
   surface that survives fast mode" slogan stated unconditionally in the same help entry.

5. **The engine exemption from the analytic joint-test screen** (CS-15, CS-41 site (i), CS-99). The screen
   refuses analytic clustered joint tests because the covariance is degenerate at cluster ≥ `i`, and its own
   text recommends the engines; the engines do not use that covariance for those tests. *Not defensible:* the
   silence when the waiver fires, and the caution-only treatment of the `boottstat` pairwise contrasts, which
   *do* use it (CS-41 site (ii)).

6. **`permute()`'s wider `subgroups()` screens** (CS-52, CS-106). Any cluster can be assigned treatment in a
   replication, so a subgroup variable that is only well-formed on today's treated units is genuinely
   ill-formed for the randomization. The messages state the reason. The stumble is expectational only.

7. **BCa computed only on request** (CS-11, CS-35, CS-76, CS-109, CS-115). The delete-one-cluster jackknife
   costs G extra fits, and the cell surface would need a jackknife per cell. Both are deliberate (F067).
   *Not defensible:* the four printed and written places that promise every family is computed regardless,
   and the un-flagged bc/bca disagreement between `e(boot_display)` and `e(boot_gteffects_family)`.

8. **The cell collapse's ineligibility refusals have no `safety(off)` escape** (CS-49, CS-104). Each entry is
   a *demonstrated wrong answer* rather than an accepted risk, which is a different category from every other
   screen in the package. *Not defensible:* that the same feature also skips silently by default, so one
   feature holds the package's strictest and its least communicative policy at once.

9. **The bootstrap reposts a full `e(V)` and aligns `e(df_r)`; the analytic path posts a diagonal**
   (CS-4, CS-59, CS-77, CS-126). The bootstrap has a genuine cross-estimand covariance and the analytic path
   does not; the `df_r` alignment is deliberate and measured. *Not defensible:* the single un-engine-aware
   help paragraph, the absent provenance marker, and `e(df_r_analytic)` going missing on the default clustered
   branch.

10. **`t` versus `z` is chosen by the estimator analytically and by the reference distribution under the
    bootstrap** (CS-18, CS-88). The two engines genuinely have different things to name — one an estimator,
    one a critical-value family. *Not defensible:* that nothing says which fact drives the letter.

11. **`cluster()`'s four meanings** (CS-27). Each is right for its feature, and the differing requirements
    follow from the mechanics — a variance can be clustered finer than `i`, a randomization unit cannot.
    *Not defensible:* the divergent error wording across the four sites, and that no run names the resampling
    or randomization unit it used.

12. **The `e(cfx*)` matrices keep margins' `r(table)` vocabulary** (CS-66, CS-79). `b se t pvalue ll ul` is a
    deliberate "this came from margins" signal, and renaming would break the `saveplots` dataset contract.
    *Not defensible:* that the cfx bounds alone carry no level marker while every sibling surface does.

13. **Footer blocks that exist under one engine only** (CS-91, CS-108). A permutation step-down has no
    bootstrap twin today and a percentile interval's asymmetry has no permutation twin; those are real
    differences between the methods. *Not defensible:* the bootstrap's contrasts sentence printing on runs
    with no `subgroups()`.

14. **`reps()` has no floor and `permreps()` does** (CS-20, CS-44, CS-89, CS-123). A permutation p has a hard
    combinatorial floor set by the design; a bootstrap p has only a Monte-Carlo one that shrinks with
    `reps()`. Both are standing author rulings (F030, F036) pointing in opposite directions for a reason.
    *Not defensible:* that neither entry mentions the other, and that the option names give no hint they are
    the same concept.

## 5. Refuted

**No finding was refuted outright.** All 132 items survived verification in substance; 21 of the 48 verified
items lost or narrowed a sub-claim — the 17 scored "partly", plus four that were confirmed with a
correction — and those are listed here so nothing is silently carried forward. Each line names the
claim as filed and why it does not stand; the corrections are already applied in section 2.

1. **CS-1** — "four different p-value rules". Three formulas, not four: the bootstrap's joint family borrows
   permute's sampled add-one rule rather than adding a fourth.
2. **CS-1** — "permute uses the add-one Monte-Carlo form". Not under `permexact`, which uses the bootstrap's
   own plain `2*min/n`; permute's exact floor comes from self-membership, not from add-one.
3. **CS-1** — "only the footer mentions add-one, and only for the joint family". Understated as filed
   (nothing printed named it under `boottstat`) and now overtaken: `b6ef398` added the qualification to the
   footer on 2026-09-06.
4. **CS-3** — "silently absent ... only discovered when a downstream command refuses". The permute header
   and footer both announce the fast-mode state; what is undisclosed is the lost `e()` surface.
5. **CS-5** — the citation ado:12218 as the fixed-metric adversarial substrate. That line is the
   `boottstat` studentized construction; the substrate is `wooldid_permadv` at 13377-13425.
6. **CS-5** — that extending the fixed-metric test to `adversarial` is purely additive. The permutation
   construction tests a different null (sharp no-effect-at-any-horizon), so it needs its own footer sentence.
7. **CS-8** — "no disclosure at all" for a `permrawp` zero. sthlp:481 states the floors including "0 under
   permrawp"; the gap is on-screen, not documentary.
8. **CS-26** — "`se` means three different quantities". Two under the label `se`; the third is
   separately marked `se_analytic`.
9. **CS-26** — "cannot detect percentile-rank inference". Not from the colname stripe, but
   `e(boot_display)` and `e(boot_stat)` are posted and documented.
10. **CS-22** — `permtci` implies `permtstat` "silently", and `bootdf()`/`ties` under `boottstat` "provably
    do nothing". The implication is documented at sthlp:457 and the run is labelled in header and footer;
    both options are *display*-inert, not inert — `bootdf()` still sets `e(df_r)`, `ties` still writes the
    stored BC tables.
11. **CS-24** — `emptycellsoverride` as an orphan modifier of `verbose`. It is neither a modifier of
    `verbose` nor inert: it suppresses `set emptycells drop` unconditionally, changing the design matrix.
12. **CS-25** — that `contreatwithin` lands the user in the collapse-fallback case as well as the
    withheld-plots case. Because of the `safety`/`safetyinput` split it reaches the plot case only.
13. **CS-37** — that `contreatwithin` disables the multicollinearity error against the help's rule. That
    relaxation is documented and intended, and the error's own text prescribes `contreatwithin` as the
    remedy; only the maxvar/granular-`i` preflight at 1446 falls on the help's "still runs" side.
14. **CS-39** — "permute() then says nothing at all". Two cautions exist (`permtstat` at `cluster() == i`,
    and the sampled-draw p-floor); both are gated too narrowly to reach the user who arrives from the
    bootstrap's referral.
15. **CS-40** — "sthlp 977-987 describes only the first sense". The suppress-the-check sense is documented in
    roughly nine further passages and in two of the cautions' own printed text; the defect is scatter plus the
    missing run-time list.
16. **CS-40** — that running the multicollinearity scan under `safety(off)` costs the interpreted loop. Since
    the 2026-09 rework the loop is short-circuited by `if strpos("`bvars'", "o.") != 0` when no coefficient
    carries an omission flag.
17. **CS-41** — that the `permtstat` caution at 5902 is the same claim as the `jointtests` refusal. It fires
    only at `cluster()` exactly equal to `i` and its stated harm is different; it is a related third
    treatment, not a third instance.
18. **CS-42** — "permute posts two blanking flags". Five, and two of them sit beside condition-type cautions,
    so the post-plus-footer pattern is not confined to display-blanking cases.
19. **CS-6** — that permute's blank contrast `t_perm` column is undocumented and uncautioned.
    sthlp:454 and sthlp:413 state the exception and ado:10062 half-signals it; the residual defect is the
    column label plus two generic statements that contradict it.
20. **CS-57** — that `e(perm_lb_att)` "prints beside the estimate" as a CI bound. The footer and help both
    label it as the null band; the risk is programmatic only, which is why the priority was softened.
21. **CS-59** — that the fast-mode `e()` absence is footer-prose only. sthlp:417 documents it; and
    `e(b_es_*)` is withheld too, which the finding did not claim.
22. **CS-60** — that the `allsgssame` existence idioms at 6733/7862 need rewriting under pathway (a). They
    reference a scalar and are unaffected; the site that does need editing is the display block at 9928-9942.
23. **CS-61 / CS-119** — "reading the help, they do not learn that `e(perm_pconv)` can say 'raw' at all".
    sthlp:475 lists all five values; the lag is in the stored-results block only.
24. **CS-97** — "no indication that a correctness screen was bypassed", and "the only hint is a `boottstat`
    caution". The containment screen prints its own caution recommending bootstrap/IM/RI in the same run, so
    this is a half-disclosed design rather than a silent bypass; the 6433 caution is a different surface. The
    item's `defensible` flag flips to true on this basis.
25. **CS-115** — that the stale "every family" claim sits in two places. Five, and sthlp:1061 — cited as the
    correct version — itself undercounts once BCa is displayed.

Several findings also gained sub-claims under verification, which are folded into section 2 rather than
listed here; the largest are the `e(df_r_analytic)` branch gap (CS-4), the `fastnose`+`permanalyticse` hole
(CS-23), the F036 caution's `absp`-only floor arithmetic (CS-8), `e(boot_adv_absent)` being posted and never
read (CS-75), the fourth stale "every family" site at ado:6149 (CS-76, CS-115), the third `level` default
inside `wooldid_ptp` (CS-62), the `latticeignoreweights` exception to bootstrap bin drift (CS-96), and the
display block at 9928-9942 that pathway (a) of CS-60 must touch.
