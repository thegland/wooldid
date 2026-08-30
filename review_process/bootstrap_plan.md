# Bootstrap SEs for wooldid — implementation plan (2026-08-29)

Author chunk: traditional cluster bootstrap + Bayesian (Dirichlet-reweighting) bootstrap,
borrowing reps()/ties/level()/dots()/rseed() from Stata's `bootstrap`/`bayesboot`;
non-cluster bootstrap if easy; EVERYTHING bootstrapped incl. joint-test p-values; built-in
compatibility ideal, DIY acceptable if provably matching. Plan authored by the coordinator
(F6 precedent); grounded in a 5-agent empirical research wave (artifacts in the session
scratchpad) plus the coordinator's own probes. Line numbers = HEAD (032d05a).

## 1. Research findings the plan stands on (all empirically verified)

**bayesboot reverse-engineered bitwise (r1).** bayesboot = `rwgen bayes` (generates reps()
Dirichlet weight VARIABLES) + `bootstrap, iweights(w1...wR)`. Weight algorithm, reproduced
max|diff|==0: rows grouped into atoms by id() (ascending sorted order); per rep, per atom
draw e_k ~ Gamma(n_k + priorpowers + 1, 1) from the seeded mt64 stream (default
priorpowers(-1) => alpha = n_k); atom weight p_k = e_k/sum; EACH ROW gets p_k/n_k.
rseed() == `set seed`. Cluster-level BB is reachable only via id(clustvar); its DEFAULT is
alpha_c = n_c (cluster-size Dirichlet), NOT the equal-alpha Rubin BB; equal-alpha needs
priorpowers(pvar = -n_c). **`ties` is NOT a resampling mode**: it only adjusts the BC/BCa
median-bias z0 for ties between replicates and the observed value (a no-op for continuous
stats); implementing it as a weight-mode switch would contradict the built-in.
Replicate files store FLOAT unless saving(,double) (~1e-8 noise in e(V)).

**bootstrap prefix on wooldid: blocker root-caused and patched in scratch (r2, r3).**
The prefix's rep loop runs under `set coeftabresults off`, which silently removes
margins r(table); wooldid's 16 r(table) harvests then die r(503) "conformability error"
every rep. Two one-line patches make the OFFICIAL prefix work end to end (proven:
45-50/50 reps, estat bootstrap all, Wald joint tests on bootstrap V): (P1) `set
coeftabresults on` at worker entry (program-scoped, auto-restored); (P2) declare the
shell `eclass properties(...)` so bootstrap's post-rep Display routes to _prefix_display
instead of tripping the replay gate. Limits of the prefix route: **rejects ALL user
weights** (r101 even with force), bayesboot: wooldid dead (injects [iw=]; wooldid takes
only aw/pw; cannot compose with existing weights for ANY command), r322 abort if any
exp_list stat is missing on the full sample, and nonconformable (shrunken-horizon) reps
are silently DISCARDED wholesale.

**e(V)/CI conventions to match (r1, r2).** e(V) = covariance of replicates, R-1 divisor,
deviations from the REPLICATE MEAN (mse option = deviations from observed, divisor R);
e(b_bs) = replicate means; e(bias); z-based table (e(df_r) missing); e(ci_normal),
e(ci_percentile), e(ci_bc) (+bca only if requested); e(z0) = BC median bias (ties adjusts
here); N_misreps counted, "complete replications only" policy. `bstat` reconstructs the
whole estat/test suite from a replicate matrix — our DIY matched it exactly.

**Rep capture (r5).** A nested fastnose rep posts the COMPLETE named e(b) vector — pooled,
pretreatment, ES horizons, AND per-SG blocks, no missing entries (probe: 18/18 named
cols) — with e(V) absent. Column names (ES_att_SG2_pre1 etc.) are stable and relyr-based;
width is DATA-DEPENDENT (a resample can extinguish a horizon/SG/cohort) => reps MUST be
stacked BY COLUMN NAME with an explicit missing-cell policy. The official prefix's
wholesale discard of shrunken reps (measured 5/20 extra discards) is a bias we can avoid.

**Every reported joint test is a Wald functional of e(b) elements (r5)** plus a covariance
wooldid never posts: joint_horiz(_i) = {ES_*_pre# jointly 0} (the over() set includes the
structurally-omitted _t==0 reference); joint_Choriz twins; per-sg variants (ereturn
LOCALS, not scalars — posting must mirror the split); allsgssame = equality chain across
SG blocks (df = #sg-1, rank-guard precedent at 2830-2849); pwsgs rows = differences of SG
elements. So a bootstrap covariance of name-aligned rep e(b) vectors reproduces the
entire joint-test surface. NOT functions of e(b) (outside "everything"): impvalues IM
inference (different estimand under impvalues — refuse), permute, pretrendpower, gteffects
cell SEs; adversarial line p is recomputable from ES columns + boot V (defer, note).
Posted e(V) is today DELIBERATELY diagonal (3653) — posting a full bootstrap V is a
surface change to make consciously (it newly enables cross-estimand test/lincom: a
feature, but a decision).

**The estimand fork (r3; the one real statistical decision).** For `att` (and `ame`:
obs-/w_att-weighted), count-weight injection == physical duplication == old-id pooling ==
official prefix, all EXACT (reldifs <= 9e-16 incl. ppmlhdfe pw; coordinator's independent
probe agrees). For the CELL-WEIGHTED estimands (att_it, att_i, att_itime) the weight path
and the re-ided path are DIFFERENT STATISTICS (measured mreldif ~.028 on a 4-cluster
probe): w_att_it renormalizes within cell, so a within-cluster-constant multiplier cancels
out of the aggregation — the weight path holds the cell set fixed and misses the
cell-sampling variance component (risk: plausibly-wrong UNDERSTATED SEs when effects vary
by cluster). The textbook cluster bootstrap of those estimands duplicates cells (re-id);
official bootstrap does the weight path unless the user adds group() (both official
variants exist). The Bayesian bootstrap CANNOT duplicate cells by construction — for
cell-weighted estimands BB reweighting reaches only the counterfactual fit. => Section 7.

**Known adjacent defect (filed, not this feature):** the F6 rep line re-embeds the
POST-validate esprelength, which validate increments AGAIN per rep under
esfixedbaseperiod when |esrelativeto| <= esprelength (the mini-assessment radar's "esf
rep esprelength+1 parity latency"). The bootstrap rep line must emit the ORIGINAL user
value (reverse the mutation, esposttrip-style); the permute twin stays a separate ticket.

## 2. Architecture

F6-sibling, DIY engine (the built-ins cannot do weighted runs, Bayesian composition,
name-aligned stacking, or SG/joint coverage — but they remain the MATCHING ORACLES, and
two tiny hardening patches make the naked prefix work for unweighted users as a
compatibility bonus):

- Options on BOTH byte-identical syntax lines (174/4023) + the 4688 c_local blanket +
  a stray-companion gate (4311-4323 pattern); check_validate_sync.py enforces the first two.
- `wooldid_bootctl` (gates; called from validate beside the permctl call at 4324-4327).
- `wooldid_bootrun` (engine; free-standing program after wooldid_permrun's end at 5207):
  cluster map frame with frozen __wdbrow ids (NFPERM-02 frlink re-sort trap), per-rep
  draw refresh via frget, ONE materialized product weight variable (reghdfe REJECTS
  weight expressions — measured r111), nested `cap qui wooldid <repcall> [wt=prodvar]`
  fastnose calls, name-aligned accumulation (mata: rep vectors keyed by colname into a
  master matrix), failure/partial policy, Break (rc 1) PROPAGATES out of the loop
  (fixing the permrun precedent's eat-Break-per-rep wart for this engine).
- Mata helper(s) beside wooldid_permx (with mata drop guards): name-aligned stacking,
  V/CI/BC construction, Wald tests.
- Shell orchestration block parallel to F6's (86-166): re-run validate to harvest
  normalized locals; build rep option string (keep everything shaping point estimates;
  KEEP subgroups — unlike permute — so SG estimands and SG joint tests get bootstrapped;
  strip display/plots/adversarial/pretrendpower/jointtests/gteffects/permute/bootstrap
  options; fastnose always); est sto/restore observed results; append e(boot_*);
  wooldid_display part(boot); replay reprint (existence-gated => no replayv bump).
- Worker touch: widen the wdp_* echo gate (217-225) to `permute!="" | bootstrap-on` +
  one wdb_on flag (worker is at its size ceiling — nothing else goes in it). Plus the
  P1 coeftabresults defense line at worker entry (independent hardening, own commit).
- DEFERRED DISPLAY (replace-mode mechanism, settled 2026-08-30): the worker calls
  wooldid_display part(1)/(roth)/(2) itself at estimation end — under bootstrap those
  calls would print the ANALYTIC tables before the reps ever run. Fix: gate the
  worker's display calls on bootstrap-off (3 one-line conditionals, ~60B in the
  worker); after the rep loop the shell rewires e() (wdidmainresults se/t/p/ci
  columns, e(V), joint-test scalars/macros, header SE-text macro) and THEN calls
  part(1)/(roth)/(2)/(boot). Because F11 display is e()-driven only, the rewired
  tables print bootstrap inference with identical formatting, and replay works
  unchanged for free. SI-A1 falls out naturally: with bootstrap on, the worker runs
  its analytic joint tests suppressed (never trips the screen); the shell posts the
  bootstrap Wald versions under the same e() names.

## 3. The two modes

- **bayes**: per-cluster draws g_c ~ Gamma(alpha_c, 1), normalized p = g/sum; row draw
  multiplies the user weight: prodvar = user_w * draw_c (draw joined cluster->rows via
  the map frame). Alpha convention: Section 11 Q3 (recommended: equal-alpha
  Dirichlet(1,...,1), the author's literal spec; bitwise-provable against
  `rwgen bayes id(clus) priorpowers(pvar=-n_c)`). All clusters present every rep =>
  estimability-stable; failed reps ~impossible on sane data.
- **standard**: per-cluster multinomial counts (how often each cluster is drawn),
  generated by CALLING `bsample, cluster() idcluster()` on the map frame (or count
  aggregation of its output) so same-seed draw sequences can be proven against the
  official prefix; prodvar = user_w * count_c; zero-count clusters drop from the rep
  sample (resampling semantics). Under Section 7 option A, cell-weighted estimands
  instead run the re-id engine leg (frame copy + bsample + i rebuilt as group(newclus,i)).
- **non-cluster**: same engines with cluster := observation (bayes: row-level Dirichlet =
  rwgen bayes row atoms; standard: row-level multinomial). Cheap because the engine is
  cluster-var-agnostic. Cluster default: cluster() variable when present (first var if
  multiway cluster() is REFUSED v1 - Section 6 governs); no cluster() => observation-level.
- Weight-type mapping for the injected [wt=prodvar]: keep the user's declared type;
  unweighted runs use aw (pw for poisson runs — poisson+aw is refused by wooldid).
  Point estimates are aw/pw-invariant (measured).

## 4. Inference construction ("everything gets bootstrapped")

From the name-aligned rep matrix (complete reps + per-column availability counts):
- se_boot per column = sd of available reps (R-1 divisor, deviations from rep mean —
  Stata convention); e(V_boot) = full covariance (nonsingular only if reps >> dims;
  posted regardless, tests guard rank).
- Table columns se/t(z)/p/ci recomputed: z-based (df_r missing convention) at level().
  CI flavors: normal-based in the table; percentile and BC stored in e() matrices
  (estat-bootstrap trio); `ties` = BC z0 tie adjustment exactly as the built-in.
- Joint tests recomputed as chi2 Wald with V_boot on the same element sets, SAME e()
  names (scalar-vs-macro split mirrored), z/chi2 metric documented; rank policy: if the
  relevant V_boot subblock is rank-deficient or a column has fewer than 2 completed reps, post
  missing + one caution (never a silent 0.0000). The SI-A1 jointtests screen (whose own
  error text recommends bootstrap) gets a carve-out for bootstrap-based joint tests.
- Per-column missing policy: a column absent in a rep (extinct horizon/SG) contributes
  nothing to that column's moments; per-column n_reps recorded; columns with
  n_reps < 2 get missing inference (arithmetic, not a gate - Q4) + a named caution listing them. Whole-rep
  failures (r459 no-controls draws etc.) counted like the official 'x' reps
  (e(boot_failed)), never aborting the run. NO completion floor (Q4 ruling; deliberate contrast with permreps>=20
  precedent).
- e() additions: boot_method, boot_reps/completed/failed, boot_level, boot_stat=z,
  b_boot (rep means), bias_boot, V_boot, ci_normal/percentile/bc matrices, per-column
  rep counts, bootdist matrix (reps x columns, cap like permdist 10k), rngstate/seed
  record. Observed point estimates NEVER change.

## 5. Option surface (borrowed names, as directed)

Master: `BOOTstrap(standard|bayes)`. Companions (gated on it via the stray-companion
screen): `reps(#)` (mandatory explicit, NO floor — Q4 ruling), `ties`, `level(cilevel)` (10-99.99 incl. decimals,
default 95; scoped to bootstrap outputs), `dots(#)` (progress every # reps; 0=silent;
default dots(1) with a Stata-style header), `rseed(#)` (= set seed; session RNG advances,
documented like the built-in). NOT colliding with perm* names. `update` unaffected.

## 6. Gates (bootctl) and interactions

Refuse: bootstrap+permute (two rep engines); bootstrap+impvalues (different estimand;
permctl precedent); bootstrap+fastnose (user-level — nothing left to bootstrap onto...
actually reps use fastnose internally; user fastnose means no SEs are wanted: refuse);
bootstrap+update; oldsyntax; multiway cluster() v1 (single numeric cluster var, permctl
precedent); reps() must be present (no floor); >=2 clusters hard error + few-cluster caution (<10);
suppresssgprimaryeffects+bootstrap (SG tests' inputs absent from rep e(b) — refuse or
auto-plain; recommend refuse v1). unconditionalse/customvce: the OBSERVED run keeps its
VCE for the analytic table; bootstrap replaces reported inference — v1 recommendation:
refuse (contradictory SE concepts, matches permctl), revisit later. contreatwithin:
allow (weights flow through; no onset machinery involved) unless review finds a snag.
Containment/amescreen: the observed run's screens/cautions still fire and STAY; bootstrap
computes rep dispersion regardless (that is R2-adjacent dispersion — do NOT silently
present it as a validated cure: keep the containment caution alongside the boot SEs).
esprelength re-normalization: rep line emits the ORIGINAL user value (Section 1 last
finding). Nested-call hygiene: WOOLDID_LASTRC/CLEANUP end-state restored (F6 precedent).

## 7. THE DECISION: cell-weighted estimands (att_it, att_i, att_itime, ame_it)

Options (author picks; Section 1 fork):
- **A (recommended): split engines by estimand class.** att/ame legs: weight path (exact,
  provable against the built-ins bitwise). Cell-weighted legs under standard mode:
  physical re-id resampling (frame copy + bsample cluster() idcluster() + i rebuilt as
  group(newclus, i)) — the textbook bootstrap, matches official-prefix-with-group();
  under bayes mode: no cell-resampling exists — post the weight-path numbers WITH a
  documented caution that BB reweighting does not propagate into fixed equal-cell
  weights (or refuse bayes for cell-weighted estimands: sub-choice).
  Cost: a second engine leg (~frame copy per rep for those estimands only).
- **B: v1 = att/ame only**; cell-weighted estimands + bootstrap refused with an
  informative message; graduate later. Leanest, nothing wrong ships.
- **C: weight path for everything + caution.** One engine; cell-weighted SEs understate
  cell-sampling variance (measured divergence ~3% on a tiny probe; worse when effects
  vary across clusters). Risks the exact "plausible wrong numbers" class.

## 8. Prefix compatibility (bonus deliverables, independent commits)

(P1) `set coeftabresults on` defense at worker entry — makes wooldid survive EVERY
Stata prefix's replication environment (bootstrap/jackknife/simulate/permute all set it
off). (P2) shell `properties()` so bootstrap's Display works (exact property set to be
minimized at implementation; svy compat NOT claimed without testing). + an sthlp note
with the correct naked-prefix incantations (unweighted only; idcluster/group() forms).

## 9. Validation ("prove we match")

- standard/att: same-seed rep-by-rep equality vs official `bootstrap _b, cluster()
  idcluster(): wooldid-wrapper` (patched build), and count-path == duplication identity
  re-proven on wooldid output (coordinator probe already: <=9e-16).
- standard/cell-weighted (if A): vs official prefix WITH group(i).
- bayes: our draws == `rwgen bayes` bitwise under the same seed (r1 recipe replicated
  the generator exactly, incl. priorpowers mapping); end-to-end vs bayesboot on an
  unweighted twin using generate()+saving(,double).
- V/CI construction: `bstat` on our replicate matrix must reproduce our posted
  se/ci_normal/percentile/bc exactly (r2 proved the bstat route reconstructs the
  official suite).
- Joint tests: our Wald chi2 == `test` after the official prefix on the same replicates.
- Size sanity MC (feature-scoped, small): boot SE vs true sampling SD + CI coverage on
  a known DGP, both modes, G in {10, 30}; no full-suite run (author directive).
- Targeted cert suite test_bootstrap.do (gates, both modes, weighted+unweighted,
  poisson pw, subgroups, esf window, failed-rep counting, replay, determinism, Break).
- Suite/battery: new-only configs rc_old=198 pins + one battery config blessed on the
  feature build; data_classes rows for any new fixtures; patho-mode config per the
  new doctrine (bootstrap on dta_path_onepercell: observable behavior).

## 10. Commit sequence

A: hardening P1+P2 (independent, battery-neutral, own tests) ->
B: option surface + c_local + bootctl gates (loud "engine not yet enabled" tails; F6-A
   precedent) ->
C: bayes engine + name-aligned collection + e(boot_*) + table rewiring + display/replay ->
D: standard engine (count path; + re-id leg if Section 7 = A) ->
E: joint tests + CI trio (level/ties/BC) + dots polish ->
F: sthlp + configs + matching-proof suite + MC.
Assessor/verifier waves at C, D, E (high-risk feature: R1 review -> plan update -> blind
R2 per process); the author's standing Stata cap (8) and edit-serialization rules apply.

## 11. Author decisions — RULED 2026-08-30 (AskUserQuestion)

Q1 presentation: **REPLACE** reported inference vce(bootstrap)-style; point estimates
    unchanged; analytic V preserved as hidden e().
Q2 cell-weighted estimands: **Split, bayes cautioned — with the author's safety
    modification**: standard mode uses physical re-id resampling whenever a
    cell-weighted estimand is requested (weight fast-path when only att/ame);
    bootstrap(bayes) + cell-weighted estimands under DEFAULT safety **refuses** with a
    message saying it will run under safety(off); under safety(off) it **runs with
    cautions, not refusals** (gate keys on `safetyinput`, the F13/S5 precedent, so
    contreatwithin's internal safety rewrite cannot void it).
Q3 bayes alpha: **default = bayesboot's convention (alpha_c = n_c, cluster-size
    Dirichlet)**, plus an option to switch to equal-alpha Dirichlet(1,...,1).
Q4 reps: **mandatory explicit reps(), NO floor** — low-rep runs are legitimate (speed
    testing etc.); the many-reps reminder lives in the documentation. Where a statistic
    is arithmetically undefined at tiny rep counts (SE needs >=2 completed reps) it is
    posted missing, which is arithmetic, not a gate. (Deliberate contrast with
    permreps' floor-20.)

## 12. R1 adversarial-review amendments (2026-08-30; R1 verdict: amend-then-implement)

R1 confirmed the core machinery empirically (weight-injection == physical resampling
BITWISE incl. screens and horizon extinction — the zero-weight markout at L275 removes
zero-count clusters from the SCREENED sample so min_t/controlcheck/refperiod/horizon
grids recompute exactly as under physical removal; joint_horiz == Wald-F from
e(b_es)/e(V_es) to 15-16 digits; ES names relyr-anchored => name-meaning stable under
resamples; rngstate untouched by nested wooldid calls; admin-scale timing ~0.5-0.8 s/rep,
frame copy ~1%). Amendments adopted:

- **A1 (F1 BLOCKER — nesting gate).** The re-id leg is defined ONLY when the resampling
  cluster nests i (i constant within cluster — gegen minmax shape, permctl 4912-4920
  precedent). bootctl gate: bootstrap(standard) + any cell-weighted estimand requires
  cluster-nests-i; otherwise (incl. observation-level bootstrap) REFUSE v1 with an
  informative message (suggest att/ame, or cluster at/above i). bayes + cell-weighted
  keeps the Q2 safety-keyed refusal/caution regardless of nesting.
- **A2 (F3).** Boot V columns with <2 completed reps or extinct horizons: post via the
  -1e100 missing-diagonal sentinel machinery (3667-3699 precedent) — a raw repost of a
  missing-diagonal V is r(504) (measured). Sentinel columns are excluded from every
  Wald subblock.
- **A3 (F4).** BC z0 and `ties` compare replicates against a REP-0 FASTNOSE run of the
  observed assignment (NFPERM-10 precedent), transported as tempname scalars — never
  against full-run values (fastnose differs in final digits; ties would silently no-op).
- **A4 (F5).** The SI-A1 jointtests screen (worker ~1126-1165) gets a bootstrap carve-out
  keyed on the bootstrap flag: with bootstrap on, the worker's analytic joint tests are
  suppressed internally (never computed, never screened) while the jointtests REQUEST is
  recorded for the shell's Wald stage; the 1165 safety(off) caution twin is suppressed
  under bootstrap (its oversized-p warning describes the analytic tests being replaced).
- **A5 (F6).** v1 plots policy = document-and-caution: makeplots/makecfxplots/saveplots
  under bootstrap print one caution that plot CIs and saved plot data are analytic while
  the tables carry bootstrap inference; the esplots re-render-from-rewired-e() idea is
  FILED as an enhancement (feasible-looking, untested; cfx plots can never re-render
  from e()).
- **A6 (F7).** Rewired surfaces list is EXPLICIT: e(wdidmainresults) + e(wdidsgresults_*)
  se/t(z)/p/ci columns, e(V) (full boot covariance; analytic kept as hidden
  e(V_analytic)), e(pwsgs)/e(prepwsgs) (diffs of e(b) elements under boot V), all
  e(joint_*) names, e(V_es_*) REPLACED by name-aligned boot subblocks (analytic kept
  hidden) so honestdid/pretrendpower consume bootstrap Sigma_pre automatically.
  adversarial + bootstrap: re-run wooldid_advstat with the boot V_es post-rewiring
  (bounded effort; fallback = documented analytic provenance). STRETCH item, commit E.
- **A7 (F8).** n_c (bayes alpha) and the standard-mode draw universe are defined on the
  permrun smp concept (user-if sample + nonmissing cluster; byte-parallel bootctl/
  bootrun rule, NFPERM-03 pattern); the rwgen-bitwise proof pins it.
- **A8 (F9 — measured 3->4->5).** Fix form: validate c_locals the RAW user esprelength in
  a NEW blanket-listed local (esprelengthraw; safetyinput precedent) and the rep line
  emits that. Cert test must use DEEP pre-history (min_t > user espre + 1) — the defect
  is latent on short-pre panels.
- **A9 (F10).** The safetyinput capture (4332) moves ABOVE the bootctl call site, or
  bootctl receives `safety' (still raw at 4325). One-line trap on the Q2 ruling itself.
- **A10 (F11/F12/F14).** bootrun reads rep e(b) only when rep rc==0 (staleness measured);
  t->z column rename for non-poisson, lb/ub column names follow level(), header
  e(wdid_vcespec) says bootstrap; re-id leg does NOT re-id user covariates keyed to unit
  identity (official idcluster shares the limitation) — sthlp note + caution when fe()/
  interactivecontrols name the i variable under the re-id leg.

## 13. R2 blind-review amendments (2026-08-30; R2 verdict: amend findings 1-6, then
## proceed; commit B passes as a faithful A1/A8/A9 implementation, no re-review)

All R2 findings adopted (coordinator agrees with every fix as recommended; no
disagreement round needed):

- **A11 (R2-1 BLOCKER — covariance policy).** MEASURED: a pairwise-complete covariance
  over availability-varying columns can be full-rank INDEFINITE (min eig -0.042);
  invsym() then silently zeroes a negative pivot and returns a confident wrong Wald p;
  chi2tail on a negative statistic returns 1 silently. POLICY: every inverted quadratic
  form (joint-test subblocks, posted e(V_es_*) subblocks) uses COMPLETE-CASE (listwise
  within the subblock) covariance - PSD by construction - plus an explicit
  PSD/eigenvalue guard (never rank alone) with sentinel+caution on failure. The full
  posted e(V) is complete-case listwise across the observed column set with
  e(boot_N_V) recording the listwise rep count (bayes mode: == all completed reps by
  construction). DOCUMENTED consequence: a listwise subblock diagonal need not equal
  the table SE^2 (per-column availability vs listwise) - stated in sthlp, not left for
  a verifier to find.
- **A12 (R2-2 — structural-zero contamination).** (a) Every column whose OBSERVED
  reported b is missing in e(wdidmainresults) OR e(wdidsgresults_*) (amevoid and
  friends; SG blocks void independently) is EXCLUDED from bootstrap inference
  (sentinel + caution). (b) Rep-level structural zeros: in ame-family columns only
  (the mechanism = resamples losing dose variation; fastnose skips the amescreen so
  reps return exact 0.0 silently - MEASURED), a rep contributing an exact 0.0 to a
  column whose observed b is nonzero is blanked as unavailable for that column and
  counted; one caution names affected columns.
- **A13 (R2-3 — pretrendpower).** bootstrap + pretrendpower REFUSED v1 (wooldid_ptp
  computes/displays/posts inside the worker from analytic V_es before any rewiring;
  "automatic bootstrap Sigma_pre" was unfulfillable in this ordering). Message points
  users at feeding the rewired e(b_es_*)/e(V_es_*) to Roth's pretrends package
  directly. A6's wording amended accordingly.
- **A14 (R2-4 — suppression channel).** The worker's analytic joint tests are
  suppressed under bootstrap by gating the COMPUTATION SITES on the bootstrap flag -
  NEVER by setting the suppressjointtests local (that would post
  e(wdid_suppressjt)=suppressjointtests and the display would silently never print the
  boot joint tests, live or on replay).
- **A15 (R2-5 — rep line).** contreatwithin JOINS the rep keep-list (the only
  point-estimate-shaping option the F6 template omits that bootstrap allows; omission
  is silent - validate's internal rewrite emits safety(off), so reps would run the
  non-within model with no error). Cert case added.
- **A16 (R2-6 — cleanmatrices).** bootstrap + cleanmatrices REFUSED v1: the worker-end
  clear mata under cleanmatrices destroys the ado's compiled mata functions, killing
  any post-worker engine r(3499). MEASURED LIVE on permute at HEAD - filed as its own
  pre-existing defect in wooldid_standing_issues.md (predates this feature).
- **A17 (R2-7/9/10 + notes).** e(df_r): keep + document (user test/lincom after
  bootstrap are t/F-based on the rewired V and will not match the reported z/chi2 p's
  - sthlp note). Plan text corrected to the update-WINS convention (worker executes
  update first; a refusal is unreachable). Sentinel columns get zeroed off-diagonals
  in every posted V; the master column set is the OBSERVED e(b) name set (rep-only
  columns ignored). Low-availability columns (completed in a small share of reps) get
  a named disclosure caution (not a floor - Q4-consistent). Determinism recipe pinned:
  draws via bsample, cluster() idcluster() ON THE MAP FRAME (plain bsample diverges;
  map-frame draws == full-data draws same-seed, MEASURED); matching proofs compare
  fastnose-to-fastnose at project tolerance (the count-vs-duplication identity is
  ~1e-15-exact, not bitwise-0). Re-id reps that trip the SI-8 maxvar screen (r907)
  count as failed reps (sthlp sentence). A10's i-naming caution scans
  controls/coarsecohortcontrols/ccc_absorb too. A5's plots caution also names
  gteffects/histogramcohorteffects SEs as analytic. Commit-B token nits (reps(0)
  indistinguishable from unset in the stray screen; level(0)/dots(-1) sentinel
  overlap; rseed integer-only vs set seed state strings) documented in bootctl
  comments rather than churned.
