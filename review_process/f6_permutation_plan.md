# F6 implementation plan: permutation / randomization inference for wooldid
Drafted 2026-07-27 by the coordinator (no subagent delegation, per the author's brief).
Status: PLAN FOR AUTHOR REVIEW — no code has been written. Line anchors reference HEAD
a70990a and will drift; anchors are re-verified at implementation time.

## 0. What the feature is

A `permute(...)` inference mode: hold the estimation problem fixed, re-randomize which
clusters carry which treatment-onset dates, re-estimate the target estimand(s) under each
counterfactual assignment, and compare the observed estimate to the resulting null
distribution (Conley–Taber 2011; MacKinnon–Webb few-treated-clusters motivation). Three
randomization designs (author's brief): count-preserving permutation, within-strata
permutation, and per-cluster assignment probabilities. Outputs: permutation p-values,
percentile CIs (5th/95th of the permuted-effect distribution by default), the SD of the
permuted-effect distribution (labelled SD), an optional inverted-t CI, and always the
number of possible assignments plus the number of reps used.

## 1. User-facing surface (proposed; author to confirm names)

New options (all parsed in BOTH syntax lines + returned via the c_local blanket, per S5):

- `permute(string)` — master switch + method selector. Grammar:
  - `permute(shuffle)` — method 1 (count-preserving permutation across ttre levels).
  - `permute(shuffle, strata(varname))` — method 2.
  - `permute(probabilities, probstring(...))` — method 3. (The probability string is
    long; it gets its own option `permprobs(string asis)` rather than living inside
    permute() — parentheses inside nested options are fragile in Stata parsing. See §5.3.)
- `permreps(integer)` XOR `permexact` — reps with replacement vs full enumeration.
  Method 3 requires `permreps()` (author's brief).
- `permtstat` — randomize on t-statistics instead of raw coefficients (disables the
  fastnose fast path for the reps; author's brief).
- `permtci` — additionally report the inverted-t-stat-style CI (see §7.3).
- `permseed(integer)` (proposed addition) — convenience seed; if omitted we use the
  session RNG state as-is and record c(rngstate) start in a hidden e() macro so runs are
  reproducible after `set seed`. Cheap, avoids "why do my p-values move" support burden.
  DROP if the author prefers a leaner surface.

Estimand scope: permutation statistics are computed for EVERY requested estimand
(att/att_it/att_i/att_itime and their ame twins under contreat) — each rep's wooldid call
already produces all of them at once, so multi-estimand support is free. Display and e()
are per-estimand, mirroring how joint tests are reported.

## 2. Validation gates (wooldid_validate; all r(198) with informative messages)

1. `permute()` requires `cluster()` with a SINGLE variable (the randomization unit).
2. Cluster nesting check (author's brief): every cluster must sit inside exactly one
   ttre level — i.e., the effective onset (after anticipation) is constant within
   cluster, with never-treated (ttre==.) a level of its own. Implementation: gegen
   min/max of ttre by cluster over the estimation-eligible sample; refuse on any
   mismatch (message names an offending cluster). This is validate-time work on the raw
   variables, mirroring the F13 within-i-constancy screen.
3. Refuse `impvalues`, `unconditionalse`, `contreatwithin` (author's brief).
4. Refuse `permexact` + `permreps()` together; require exactly one; method 3 refuses
   `permexact`.
5. `strata()`: clusters must be strict subsets of strata (same gegen check, cluster
   within stratum); strata with zero assignment freedom (all clusters in one ttre
   level) are legal — they just contribute a factor of 1 to the count.
6. Sanity floor: at least 2 clusters and at least 2 distinct ttre levels among the
   clusters being shuffled; refuse otherwise ("nothing to permute").
7. Interactions I propose to ALLOW in v1 with documented semantics: subgroups (perm
   stats for sg-specific estimands are NOT computed in v1 — only the pooled estimands;
   the sg machinery still runs on the observed estimation), esf/espre/espost (window
   options apply to every rep), makeplots/saveplots (observed run only), gteffects
   (observed run only), timetrends, anticipation (the permuted variable is the RAW
   user ttre; each rep's wooldid call re-applies anticipation itself), poisson.
   Interactions to REFUSE in v1 pending demand: `oldsyntax` (cmdline rebuild
   complexity), `update`. contreat: ALLOWED — the dose variable stays with the unit
   and only onset timing is permuted; the plan flags for the author that under
   permute + contreat the null being tested is "timing is as-if random given dose",
   which the help file must state plainly.

## 3. Architecture: where the loop lives

**Orchestration in the SHELL (program wooldid), after the observed worker call.**
Rejected alternatives: (a) inside the worker — recursion through the shell from within
the worker collides with the WOOLDID_CLEANUP/WOOLDID_LASTRC globals mid-flight and would
put a 3,000-line loop inside the size-critical program; (b) a bespoke lean re-estimation
loop that bypasses wooldid — maximum speed but duplicates the weight/margins pipeline,
the exact thing this codebase's history says not to do.

Flow in the shell when `permute()` is on:

1. Run the observed estimation exactly as today (full inference unless the user asked
   for fastnose themselves). Its e(), display, plots, gteffects all land normally.
2. Store the observed results: `tempname obsres; _estimates hold`? NO — `est sto
   `obsres'` on the posted e() (eclass store of wooldid's own post works; F11 relies on
   e() being complete). Also capture the observed statistic vector (per estimand:
   coefficient, and its se when permtstat).
3. Build the cluster→onset map ONCE (one row per cluster: cluster id, observed ttre
   level, stratum, probabilities row for method 3) in a tempname frame.
4. Rep loop (r = 1..R): write the rep's assignment into a tempvar onset column in the
   MAIN data (joinby/merge from the map frame — one `qui merge m:1` or a frame-linked
   recast; implementation detail benchmarked at build time), then
   `qui wooldid y i t <ttre_rep> [weights], <rebuilt options> fastnose` (or full when
   `permtstat`), read the per-estimand stats from e(), append to a results frame.
   Failures (r(459)-class: a rep assignment with no control group, collinear cells,
   etc.): `cap`, count as failed, continue; failed reps are EXCLUDED from the null
   distribution and reported (see §7.5).
5. Restore: `qui est restore `obsres'' (observed e() back), clear WOOLDID_LASTRC to the
   observed run's success state, then append permutation results via `ereturn` (pure
   post-restore adds, F2/F8/F3 precedent) and print the permutation display block.
6. The rebuilt rep command line: the shell owns `0' (the original invocation). Rebuild
   = replace the ttre token (4th word of the varlist... NOT robust under `if/in/weights`
   syntax — instead re-parse `0' with the same syntax line the shell already shares
   with validate is NOT available in the shell (it is deliberately thin). Robust plan:
   the WORKER, on a permute run, c_locals back the fully-normalized component pieces
   (varlist tokens, weight expression, option string minus permute*/plots/gteffects/
   summarystats/verbose) for the shell to reuse. This adds one c_local block at the end
   of the worker, built from locals it already holds. Rep calls also strip:
   summarystats, jointtests, adversarial, pretrendpower, makeplots/makecfxplots/
   saveplots, gteffects, histogramcohorteffects, verbose (each is pure waste per rep).
7. Nested-shell globals: rep calls run AFTER the observed run is complete and stored,
   so WOOLDID_CLEANUP churn is sequential and self-contained per rep; the loop clears
   WOOLDID_LASTRC after the last rep (observed run succeeded, so the end state is
   empty). The obsres store must NOT collide with the reserved names (tempname).
8. New subprograms (file-level, keeps the worker and shell small):
   - `wooldid_permctl` — validate-time deep checks + probability-string parsing
     (called from wooldid_validate; returns parsed structures via c_local/matrices).
   - `wooldid_permrun` — the shell-side engine: builds the map, enumerates/samples
     assignments, runs the rep loop, computes the summary statistics, ereturns, prints.
     (Shell stays thin: one gated call.)
   Both new programs join the rename recipe in agent_brief.md.

## 4. The three randomization engines

Common representation: the map frame has one row per cluster c with observed level
L(c) ∈ {l1..lK} (K includes the never-treated level when present among clusters).

### 4.1 Method 1 — count-preserving shuffle
The assignment vector is a permutation of the observed multiset {L(c)}. Draws: Fisher–
Yates shuffle of the level labels across cluster rows (gtools/`sort` on runiform is
O(C log C); C = #clusters is small by assumption of this method's use case).
Exact: the set of DISTINCT assignments is the multiset permutations, count
M = C! / (n1! n2! ... nK!), computed in logs via lnfactorial() (exact in double for the
counts involved; presented via %21.0g and a string macro when huge). Enumeration when
`permexact`: standard multiset-permutation recursion in Mata over the cluster rows
(deterministic lexicographic order); refuse-with-caution path per the brief: if M >
1000, print the caution and — decision needed — either proceed (could be slow) or
refuse. PROPOSAL: proceed when M <= 10,000 with the >1000 caution; refuse above 10,000
with a message directing to permreps() (author to confirm both thresholds).
Observed-assignment handling: the observed assignment is included in the enumeration
(standard practice; makes exact p-values conservative-valid: p >= 1/M).

### 4.2 Method 2 — within-strata shuffle
Identical machinery per stratum: shuffle level labels among the clusters of each
stratum independently; M = product over strata of the per-stratum multiset counts
(sum of lnfactorial terms). Exact enumeration: nested lexicographic product across
strata, same 1000-caution/10,000-refusal proposal. This forces each ttre level's
composition with respect to strata to be preserved (the author's stated intent).

### 4.3 Method 3 — per-cluster probabilities
Parse `permprobs()` per the author's grammar: first paren = ttre varname + level order
(levels validated against the data's actual levels, "." accepted); each later paren =
(cluster_value p1..pK) with sum-to-1 (tolerance 1e-6) and every cluster of the data
present exactly once. Storage: a C x K probability matrix. Draws: per cluster,
categorical draw via runiform() against the row's cumulative probabilities.
"Possible assignments" for the caution/report: the count of support combinations
Π_c |{k : p_ck > 0}| in logs. `permreps()` mandatory. NOTE for the author: under
method 3 the assignment distribution is generally NOT exchangeable with the observed
assignment, so the p-value is a Monte-Carlo tail probability under the DESIGN
distribution, not a classical randomization test unless the observed assignment is in
the support; we will caution when the observed assignment has zero design probability.

## 5. Statistics, p-values, and outputs

### 5.1 Rep statistic
Default: the rep's point estimate b_r per estimand (fastnose rep runs). With
`permtstat`: t_r = b_r / se_r from full rep runs (margins SEs), skipping fastnose.

### 5.2 P-values
Two-sided, centered at the null of no effect, standard permutation form:
p = ( 1 + #{r : |s_r| >= |s_obs|} ) / (R + 1) for sampled draws (add-one both sides,
includes the observed as one draw — standard Monte-Carlo correction), and
p = #{a : |s_a| >= |s_obs|} / M for exact enumeration (observed assignment is a member).
(Author call: also report one-sided p's? PROPOSAL: no; keep the surface small.)

### 5.3 CI and SD (author's brief)
- `perm CI` = 5th and 95th percentiles of {s_r} (the permuted-effect distribution),
  labelled exactly as what it is: the middle 90% of the null distribution of the
  estimate. Displayed with a one-line note that this is a null-band, not a CI for the
  parameter (the author's "of sorts" acknowledged in-text).
- `perm SD` = sd of {s_r}, labelled SD.
- `permtci`: inverted-t-style CI: b_obs ± q*(se-analog), where the se-analog is the
  perm SD and q is the (2.5,97.5) percentile of the STUDENTIZED permutation
  distribution when permtstat is on, else the normal 1.96 — DESIGN DECISION for the
  author: the cleanest "inverted t-stat CI" is b_obs ± t*_(95th pctile of |t_r|) * se_obs,
  which requires se_obs (available: the observed run has full inference) and the t_r
  distribution (requires permtstat-style rep SEs). PROPOSAL: permtci implies permtstat.
### 5.4 e() surface (pure-add, post-restore)
Per estimand X in the run's estimand list (att, att_it, ..., ame variants):
`e(perm_p_X)`, `e(perm_lb_X)`, `e(perm_ub_X)`, `e(perm_sd_X)` (+ `e(perm_tlb_X)`/
`e(perm_tub_X)` under permtci). Run-level: `e(perm_method)`, `e(perm_reps)` (requested),
`e(perm_completed)`, `e(perm_failed)`, `e(perm_npossible)` (double; may be +inf-ish) and
`e(perm_npossible_str)` (exact string from the log-count), `e(perm_seedstate)` (hidden,
c(rngstate) at loop start), `e(permdist_X)` matrix of the draws (CAPPED: returned only
when R <= 10,000; author to confirm cap) so users can build their own summaries.
Display: a compact block after the joint tests in part(2), existence-gated on
e(perm_method) (replay reproduces it for free).

### 5.5 Determinism and RNG
The loop consumes the session RNG stream (documented; reproducible via set seed or
permseed()). Battery/goldens: any committed battery config for permute must set the
seed inside the config (the battery harness already sets a global seed at the top, and
config ORDER matters — same lesson as the F3 RNG-ordering rule: the config will set
its own seed immediately before the call).

## 6. Performance envelope (to verify empirically at build time)

Per-rep cost = one `qui wooldid ... fastnose` on the user's data (~1-2s on battery-size
panels; F9 measured 2-4x speedups and far larger on poisson ES). 1000 reps ≈ 20-40 min
on battery-size data — acceptable for an explicitly opt-in inference mode, but the plan
includes: a progress dot every 10 reps (suppressible?), a note in the help about
runtime scaling, and rep-loop stripping of every non-essential feature (§3.6). If build
benchmarks show validate re-parse overhead is material (>20% of rep time), option B is
a `permfast` internal flag that lets the rep call skip validate via a pre-normalized
cmdline — deferred unless measured to matter.

## 7. Failure & edge handling

- Rep failures: cap'd, counted, excluded; if failed/completed > 10% print a loud
  Caution (assignment design frequently produces inestimable configurations); if
  completed < max(20, R/2), refuse to report p-values (absent e(perm_*), Caution).
- Degenerate null distribution (all s_r identical, e.g. 2 clusters): p reported but
  with a Caution; SD = 0 handled (no divide-by-zero: percentile CI still defined).
- User Break mid-loop: the loop runs inside the shell AFTER the observed run; a Break
  kills the rep call (its own shell epilogue cleans up) and the outer cap sees rc 1;
  the orchestrator then restores the observed e() and EXITS with rc 1 (no partial perm
  results posted). The observed run's results thus survive a Break during reps.
- The permuted tempvar and the map/results frames are tempnames (auto-cleaned; the
  s6-verified epilogue semantics cover the rep calls themselves).

## 8. Do-no-harm posture

- permute() OFF: zero behavioral change anywhere — every addition is gated on the
  option; battery 28/28 and goldens must be bit-identical (by construction; verified).
- permute() ON: the observed run's e()/display are IDENTICAL to the same call without
  permute() except for the added perm block/e(perm_*) (additivity test mirrors F3's:
  full-surface dump diff, only perm additions + cmdline differ).
- Size: shell gains ~10 lines; worker gains the c_local echo block (~10 lines, must
  respect the ~6-7 kB compiled slack — measured at build); the engine lives in new
  programs with their own budgets.
- S5: two syntax lines + c_local list extended for every new option; checker must pass.

## 9. Acceptance plan (suites to be written with the feature)

test_permute.do:
- P1 method-1 sanity on a constructed panel with a KNOWN strong effect: p small; with
  a placebo (no effect), p uniform-ish (single-run sanity: p > 0.05 in a seeded draw).
- P2 exact enumeration on a tiny panel (C=6, two levels, M=20): enumeration count
  matches the closed form; p matches a hand-computed value bit-for-bit; observed
  assignment appears exactly once.
- P3 exact-vs-reps consistency: R=2000 sampled p within Monte-Carlo tolerance of the
  exact p on the same panel.
- P4 strata: within-strata counts preserved in every draw (assert per-draw); M equals
  the product formula; a stratum with no freedom contributes factor 1.
- P5 probability grammar: round-trip parse of the author's example string; sum-to-1
  violation refused; unknown cluster refused; zero-probability observed assignment
  cautioned.
- P6 gates: each refusal (impvalues/unconditionalse/contreatwithin/cluster-not-nested/
  missing cluster()/permexact+permreps/method3+permexact) fires r(198).
- P7 additivity: same call ± permute() — observed e() identical except perm additions.
- P8 cautions: >1000-assignments caution; reps>10%-of-M caution; failed-reps caution.
- P9 replay: bare wooldid reprints the perm block; F11 parity from the banner.
- P10 determinism: same seed → identical e(perm_*); different seed → different draws.
- Battery: one pinned permute config (small R, seeded) added post-feature (C1 pattern).

## 10. Risks and kill criteria

- R1 (correctness of reuse): rep calls exercise the full wooldid pipeline; any
  path-dependent state leak between reps (globals, RNG, frames) corrupts the null
  distribution. Mitigation: P10 determinism test + a back-to-back-identity probe (run
  the same assignment twice mid-loop, must match bitwise). KILL if rep runs are not
  bit-reproducible under fixed assignment + seed.
- R2 (runtime): if the per-rep cost on realistic panels makes R=1000 impractical
  (>2h), escalate to the author before building the lean path (§6 option B).
- R3 (contreat semantics): if the author judges the "timing as-if random given dose"
  null misleading, gate permute+contreat behind an explicit acknowledgment option or
  refuse in v1.
- R4 (scope creep): sg-specific perm stats, one-sided p's, per-horizon perm bands are
  all OUT of v1; listed for future graduation.

## 11. Build sequence (each step gated as usual)

1. Commit A: validate gates + option plumbing (S5 both lines, c_local, permctl parser
   with its own unit probes) — inert without permute(); battery bit-identical.
2. Commit B: the engine (wooldid_permrun; methods 1+2 with reps; exact enumeration;
   p/CI/SD; e() + display + replay) + test_permute.do P1-P4, P6-P10.
3. Commit C: method 3 (probability grammar + draws) + P5.
4. Commit D: permtstat + permtci.
5. Battery config + sthlp section + solved-issues/plan-doc migration.
6. Verifier wave (3: correctness-of-inference lens incl. hand-computed exact case;
   do-no-harm sweep; edge/failure-handling lens).
