# Issue-and-Suggestion Files: Rejected Tickets

Responses to the external review tickets in `wooldid_external/issues_and_solutions/` that were
**rejected** (the issue is not real, the behavior is deliberate adjudicated design, the ticket was
already fixed, or the proposed cure is worse than the disease). Accepted tickets were implemented
in the `bugfixes` commit chain (see `git log` from b42b71f onward) and are not listed here.

Process: every rejection below was triaged against HEAD, checked against the project's
adjudication records (`wooldid_solved_issues.md`, `review_process/f6_permutation_plan.md`, the
help file), and — for medium/high-risk tickets — probed by one or two independent adversarial
reviewers per the process brief (high-risk rejections got a second, blind round). Where a
reviewer corrected my rationale, the corrected rationale is what appears here. Several
rejections adopt a small severable improvement ("adopted residue") without accepting the
ticket's proposal; those are noted.

A note on a systematic pattern: many tickets wrap a small real observation in heavy machinery —
versioned manifests, audit matrices, transactional snapshot/rollback publication, hash
registries, status enums (`NFPERM-12`, `NFRPL-01/03`, `NFES-02/03`, `NFPTP-05/07`, `CORE-10`
parts 2-3, and others reference each other's unbuilt contracts as dependencies). This apparatus
was rejected wholesale as contrary to the project's lean-fix standard; each underlying
observation was still judged on its own merits, and several were accepted in lean form.

---

## CORE-06 — IM degeneracy guard rescaling (medium)

**Ticket:** the Ibragimov–Müller degeneracy guard (`se < 1e-10*|b|` at the six ttest sites) is
wrong: it misses the exact all-zero case, can null legitimately precise results, and scales with
cluster count. Proposes a spread-based classifier with a `16*c(epsdouble)` tolerance.

**Rejected: the guard is measured, adjudicated design; the replacement is worse in both
directions.** The exact predicate was set by a three-assessor wave (commit 16263da) from
measurement: roundoff dispersion floor ≈1.6e-15 relative (worst ≈1e-12 under loosened `regtol`),
smallest realistic dispersion ≈1e-7, threshold placed between. The "missed all-zero case" is the
deliberately protected class — structural zeros (est=0, se=0) are endemic to
continuous+impvalues pre/AME horizons on ordinary data, and the no-absolute-floor design exists
precisely so they are never flagged (in-code comment at the guard). The ticket's spread==0 rule
would fire the caution on clean realistic runs and flip full b/V posts onto the 47defa0
sentinel path (moving verified baselines); its 3.55e-15 tolerance sits ~4.5 orders below the
adjudicated threshold, so the *verified* realistic-degenerate band (se/|b| ≈ 1.4e-15 at the
continuous sites; roundoff band up to ~1e-12) would silently pass and report t≈1e14, p≈0 —
reinstating the exact absurd-confidence artifact the author asked the guard to kill. The
ticket's own verification table concedes the point: it classifies a 1.3e-13 relative spread
(inside the measured roundoff band) as non-degenerate. Residual all-zero exposure is loud-ish
(ttest returns t=./p=.) and was accepted at adjudication. Reviewer: clean AGREE.

## VAL-G06 — weight-sign screen domain/bypass (medium)

**Ticket:** the weight-sign screen runs on the raw `if()` domain before markout (a negative
weight on a never-estimable row errors falsely) and is bypassable via `safety(off)` (negative
rows silently dropped). Proposes deleting the validator screen and hard-checking in the worker
after a split markout.

**Rejected on the small-benefit leg.** The cured defect class is a loud, narrow, `if()`-escapable
false positive: negative weights can never enter estimation on any path (default safety errors
r(411) at parse; under `safety(off)` the worker markout removes w<=0 rows before any
onset/reference/support construction, and the drops are counted in the displayed N line and
`e(N_dropped)`). No silent wrong result exists. The ticket's specific fix also carries verified
regressions: it moves the screen after `cleanmatrices` wipes user matrices/mata (a
pre-registered property of the validator design), changes multi-error first-message precedence,
adds bytes to the size-capped worker, and converts documented opt-out runs into hard errors with
no escape. Precedent note (reviewer-corrected): B23 (8d7bd01) rejected the `ifvar==1` domain
under the *existing* markout order (vacuity trap); the ticket's split-markout variant is a new
mechanism B23 never adjudicated — the rejection rests on the benefit test plus this screen's own
explicit bypass adjudication (f187d7d, in-code comment), not on precedent coverage or on any
claim that all screens share one gating convention (some input screens are unconditional).

## VAL-G08B — commensurability check non-bypassable duplicate (low)

**Ticket:** gating the t-ttre grid check on `safetyinput` is a defect because `safety(off)` can
let a common fractional offset run with the onset period mislabeled; wants an unconditional
duplicate assert in the worker.

**Rejected: informed-waiver design, twice adjudicated.** The screen was added and then
deliberately re-gated on `safetyinput` (f187d7d: internal relaxations must not bypass, explicit
user opt-out does), refined again with pinned tests (3151f6e), and the help documents the bypass
and its consequence in the same sentence. Deliberately-pathological input under an explicit
documented waiver is out of scope by project standards; the duplicate would also add worker
bytes and contradict the documented contract.

## VAL-G09B — onset-constancy hard invariant (medium)

**Ticket:** the F13 onset-constancy check should be unconditional and should validate the final
retained sample rather than the candidate domain.

**Rejected.** (1) The candidate domain is the *correct* surface: the row-by-row `_t`/`_i`
construction and the trim logic itself read candidate-row `ttre`, and the tranche-12
mint-staleness bug proved pre-trim violations drive downstream state — by final-retained-sample
time the construction has already consumed the data. (2) The bypass is triply deliberate (the
error text instructs `safety(off)`; the help documents it; tranche-12 invested a unit-level-
correct pooled branch plus `test_mint.do` pins to make the sentinel-coded `safety(off)` workflow
produce provably correct results). (3) The hard-invariant half would break the documented
retained-flagged-group escape (mixed missing/nonmissing recorded-from-onset data, two units
sharing an `i`, float-ID pooling — all currently run to completion under the documented waiver),
and the deferred-check half would let currently-erroring default-safety sentinel runs newly
succeed silently on the esf branch. Reviewer-corrected rationale recorded here.

## VAL-G13B — cfxplottypes canonicalization (low)

**Ticket:** numerically equivalent horizon spellings (02, 2.0, 2e0, -0) are not canonicalized;
proposes display-format canonicalization before dedup.

**Rejected: explicit recent author adjudication.** The same proposal was scoped by the author in
commit c48bdaa ("lightweight per the author's scoping"): the leading `+` is stripped; "the
02/2.0 spellings were deliberately left alone per the note." The current failure mode for 2.0 is
a loud crash, not silent wrongness. (The sibling ticket VAL-G13A — garbage keyword tokens and
fractional entries silently dropped — was real, unadjudicated, and was accepted and implemented
in 270854d with the `+` strip retained and the 2.0 status quo pinned in tests.)

## VAL-G15B — subgroup mixed defined/missing cells (high)

**Ticket:** the within-cell subgroup purity check cannot detect mixed defined/missing membership
({1,1,.} passes because gegen min/max ignore missings); proposes an unconditional structural
error on any active mixed cell.

**Rejected (two blind review rounds, both agree).** Missing-subgroup rows belonging to no
subgroup is the documented mechanism of the rule, not an oversight — the in-code comment at the
check states the exclusion, the F13 adjudication recorded "Missing sg stays allowed (help-file
design)", and the help documents it twice. {1,1,.} passing is the direct logical consequence of
that adjudicated domain decision, and the resulting estimand (missing rows ignored at margins)
is exactly what the help describes — no silent wrong result. Converting a documented allowance
into an unconditional error 459 is an estimand-policy change out of triage scope, and the
proposal is unsound on its own record: its appended adversarial review found a fatal
activation-predicate error (positive weight on a missing-membership row falsely activates cells
contributing to no returned subgroup), and it depends on the G15A union machinery whose domain
blocker its own addendum concedes. If the author ever wants signal here, the lean shape would be
a Caution counting mixed cells — filed as an author question, not adopted.

## VAL-G15C — subgroup purity check safety bypass (high)

**Ticket:** `safety(off)` bypassing the within-cell purity error is a high-severity defect; the
check should be unconditional.

**Rejected (two blind review rounds).** The bypass is author-adjudicated design advertised in
the error text itself ("To ignore this error, specify safety(off).") and documented in the help;
the gate is on the pre-overwrite `safetyinput` copy, so only an explicit user opt-out bypasses.
The bypass is never numerically silent: the unconditional across-time tier shares the cell
tier's exact row domain and missing-fold, so every within-cell defined-value mixing instance
also trips the mover caution under `safety(off)`. Unconditionalizing would strand a documented,
deliberately-escapable input class (repeated cross-sections carrying row-level subgroups under a
coarser `i` — adjudicated escapable, solved_issues F13 record). The ticket's own addendum calls
its proposal "a breaking estimand-policy change presented as a validation correction" and finds
two fatal design blockers. **Accepted-risk residual (recorded for the author):** under
`safety(off)`, the mover caution's "by-sg att/att_it are unaffected" phrasing is imprecise for
the within-cell mixing class; a one-line caution-text improvement is a below-the-line candidate
(the current text is an author-directive artifact and was left untouched).

## NFPERM-01 — failed-assignment accounting contract (high)

**Ticket:** excluding failed replications and dividing by the successful count silently converts
the advertised p-value into a conditional-on-success tail share; proposes a two-state result
contract, redraw-until-B-successes sampling with per-draw RNG isolation, and manifest/transaction
plumbing.

**Rejected; two lean residues adopted.** The failure accounting for sampled mode (excluded,
counted, displayed on every run, cautioned above 10%, floored at max(20, half)) is the
documented, planned design. The ticket's redraw scheme still samples the conditional-on-success
law — its own adversarial section concedes this — so it replaces documented semantics without
fixing the conditioning it complains about, and restructuring generation away from the single
upfront mata call would break the bitwise RNG-scheme pin (`ri_crosscheck`). Precision the review
process added (and which this report states rather than overclaims): the F6 plan's exact-p
formula is nex/M over the full enumeration with no exact-mode floor exemption; the shipped
nex/nuse-with-exemption behavior entered at implementation with documentation-grade author
linkage only. Accordingly: (a) adopted — exact-mode failures now always draw the loud caution
(any failure makes an enumeration incomplete) and the help states that the plain-share /
observed-is-a-member / 2-over-M-floor reading no longer applies when exact-mode failures occur;
(b) escalated to the author — whether exact enumerations with failed assignments should publish
the completed-subset share (status quo + caution) or refuse loudly (the plan-of-record reading).

## NFPERM-07 — zero-probability observed assignment (low)

**Ticket:** when the observed assignment has zero probability under the `permprobs()` design,
wooldid cautions and continues; demands a fatal refusal.

**Rejected: verbatim-planned, documented, test-pinned design.** The F6 plan states the
caution-and-continue decision and its rationale word-for-word (the p-value becomes a Monte-Carlo
tail probability under the design distribution); the help documents it; committed test P5e pins
rc==0 plus the caution text. The residual corner the ticket identifies (observed level absent
from the parsed level set entirely) was reachable only through the NFPERM-05 parser gap, which
was accepted and fixed (df70e81) — closing that route.

## NFPERM-08 — permtci is not test inversion (high)

**Ticket:** the inverted-t CI (b_obs ± q95(|t_r|)·se_obs) is not Fisher sharp-null test
inversion; demands deprecating it and building a sharp-null inversion engine with set-valued
confidence regions.

**Rejected; honesty hedges adopted.** The formula, name, and e() surface are the F6 plan's
explicit resolution of an author design decision (plan §5.3), documented exactly as computed.
The ticket's engine is a large new feature (root-finding over candidate nulls, breaking e()
renames that kill committed test P13, new options through both byte-synced syntax lines) with
hard dependencies on unbuilt contracts and on the anticipation removal — a feature proposal
masquerading as a bug fix. Adopted instead: the help and the runtime display line now state the
interval is a fixed-zero-null pivotality approximation relying on the permuted t-ratio
distribution approximating the t-ratio's distribution under shifted nulls — not a sharp-null
randomization-test inversion.

## NFPERM-13 — normalize tagged-missing onsets (medium)

**Ticket:** ordinary wooldid treats extended missings (.a–.z) in the onset variable as
never-treated, but permutation mode refuses them; proposes internal canonicalization via a keyed
level-map design.

**Rejected: deliberate, suite-pinned refusal with a one-line user cure.** The in-code comment
records the reason (level accounting treats all missing codes as one never-treated group, but
enumeration/probability matching distinguishes extended missings — naive acceptance silently
corrupts the design); the F6 landing record lists the refusal among pinned correctness
subtleties. A loud r(459) with "recode .a to ." beats silent design corruption. The proposal's
machinery is unanchorable at HEAD. Adopted residue: the plain-missing convention is now stated
in the `permute()` help paragraph itself (df70e81), so shuffle-mode users learn it before the
runtime error.

## NFPERM-14A — raw-mode unconditionalse relaxation (medium)

**Ticket:** the blanket `permute()` × `unconditionalse` refusal should be narrowed to
studentized modes, since raw permutation consumes only point estimates.

**Rejected: author-briefed, verified, documented decision.** The F6 plan records the refusal in
the author's brief verbatim; the permctl header comment repeats it; verification certified
brief compliance; the help documents the incompatibility. The refusal is loud with a trivial
user cure (drop `unconditionalse` for the permutation run — analytic SEs are not the object of
randomization inference). The ticket's guard condition also reads options that do not exist
(`permci()`). Relaxing an author-briefed refusal is an author question, not triage; it is listed
among the escalations in the final report.

## NFPTP-02 — rank-deficient PTP hard precondition (high)

**Ticket:** on rank-deficient pre-period covariance blocks the generalized inverse silently
selects an order-dependent coefficient subset; demands full rank as a hard precondition (error
459) with no `invsym()` on the PTP path.

**Rejected: explicitly adjudicated publish-with-full-disclosure design.** The F7 wave chose
disclosure over refusal (642084e added "depends on which columns are retained"; cce3236 added
"one representative of an underdetermined family"), because rank deficiency is common on
realistic shapes: rank(V) ≤ G under cluster-robust VCV, so `esprelength` ≥ #clusters yields
exact deficiency on ordinary few-cluster designs — precisely where these diagnostics matter. A
hard 459 would strip diagnostics there, and the ticket's atomicity rule would let one deficient
type kill all four types' diagnostics. The result is loud (caution + `e(roth_keff_*)` for
machine detection), not silent. Adopted residue (folded into the accepted NFPTP-07 fix): replay
now prints the rank-deficiency caution too (previously the disclosure existed only live).

## NFPTP-05 — PSD/eigenvalue validation of PTP inputs (high)

**Ticket:** PTP validates only missingness/diagonals before inverting; an injected indefinite or
asymmetric pre block yields plausible wrong numbers; proposes a versioned audit contract with an
eigen-validator and `cholinv()`.

**Rejected: the threat model is hand-injected e() matrices** — deliberately-pathological input by
project standards. Organic `V_es` pre blocks are PSD by construction, and the shipped guard set
(missing screens, non-positive diagonals, symmetrization, g-inverse with reduced df and caution)
was itself the adjudicated selection from the F7 pathology matrix. `cholinv` would change
numerics only below the 1e-10 indifference threshold where it agrees, while hard-failing exactly
where the project decided to disclose.

## NFES-03 — event-study geometry schema (high)

**Ticket:** the event-study exports lack a numeric, versioned, row-aligned geometry contract;
proposes capture maps, geometry matrices, reference/scale enums, a keyed manifest with a
validator wired into replay, and migration of internal consumers off label parsing.

**Rejected: the premise is materially false and the remedy is an enterprise-API wishlist.**
The numeric horizon map already ships: `e(escomp)`/`e(escomp_s#)` carry a documented numeric
`relyr` column row-aligned with the ES blocks (F8-adjudicated) — even on no-VCE/fastnose runs —
`e(estimate_type)` publicly distinguishes pooled from fixed reference, the fixed reference value
and the timetrends two-period normalization are recoverable from documented `e(cmdline)` and
visible as numeric gaps in the relyr column, and the timetrends ES trio is typically withheld by
the non-identified-VCV posting guard anyway. No demonstrated realistic failure exists; the proposal
self-declares a fatal integration blocker against the anticipation feature and an enum conflict
with the composition contract, and it would migrate pretrendpower off verified, 10-decimal-
tested label parsing for zero behavioral gain. The genuinely useful sliver — documenting which
geometries the honestdid handoff supports — is delivered by the accepted NFES-01 doc fix.

## NFCOMP-03 — sample units versus identified cells (high)

**Ticket:** `escomp`'s nunits counts sample units rather than identified coefficients, so under
partial collinearity omission the table silently presents a partially-identified horizon as
fully supported.

**Rejected: the central claim is refuted at default settings.** Partial omission of a treated
cohort-period dummy dies loudly with error 322 at the dropped-dummy screen before `escomp`
exists; the interacted-term case prints a loud Caution. The residual exposure requires
`safety(off)` — the documented opt-out — or `contreatwithin` (which internally relaxes that
screen); in continuous/interacted designs the screen cautions at most. The remedy is also
self-disqualifying: a stripe-parsing cell mapper and 15-column schema built on rejected column
renames (which would invalidate the additively-blessed escomp baselines), with its own
adversarial review conceding "a fatal schema gap remains" for fully-omitted horizons. Adopted
residue: a help hedge stating the nunits-equals-coefficients identity is enforced for plain
binary treated-cell coefficients under default safety and can fail under `safety(off)` or
`contreatwithin`. Routed to the author: a lean caution-demotion option (running the existing
plain-cell scan as a Caution under user-typed `safety(off)`, ~5 lines, zero numeric movement).

## NFRPL-01 — permutation failure atomic result (high)

**Ticket:** a permutation-engine failure after successful estimation clears WOOLDID_LASTRC,
restores the observed-only result, and leaves it replayable; proposes an atomic e-class
publication transaction (snapshot/rollback shell, new eclass helper, completion marker,
version-3 manifest).

**Rejected: the attacked behavior is on record as deliberate, verifier-audited design.** The F6
record states verbatim that an engine-stage failure clears the marker because "the observed
results ARE valid and restored, so replaying them is not misleading" — the user received a loud
nonzero rc for the permutation request, and result-selection semantics match the Stata norm.
Every failure path was failure-injection verified at landing. The proposal contradicts the
adjudicated thin-shell design, the 10k permdist cap, and the worker-size rationale for
shell-side orchestration. Adopted residue: one help clarification that a post-estimation
permutation-phase failure leaves the valid observed estimation replayable (the existing "parse
stage or later" sentence slightly overstated refusal).

## NFRPL-02 — remove the WOOLDID_LASTRC session global (high)

**Ticket:** the session global can collide with caller state, describes the session timeline
rather than the active result, and should be replaced by a result-bound completion marker plus
the NFRPL-01 transaction.

**Rejected: the marker's design is author-directed** (d349bb1, "author-directed replay UX
hardening" — the primary basis; the reserved WOOLDID_* namespace convention shared with
WOOLDID_CLEANUP additionally answers the collision mode, which requires the caller to write into
wooldid's reserved namespace — deliberately-pathological input). Decisively, the ticket itself
concedes standalone removal is unsafe and requires the entire rejected NFRPL-01 transaction plus
a replay-contract bump that would invalidate every stored .ster and the byte-parity replay
corpus — there is no lean version of this ticket. Adopted residue: the help's replay paragraph
now states the refusal is session-scoped (it also blocks replay of models explicitly selected
with `estimates restore` until any wooldid call completes successfully).
**Knowingly-retained limitation (adjudicated by this review, recorded so future reviews do not
re-litigate it as unexamined):** after a failed wooldid call, restoring an older complete
wooldid estimate still leaves bare-`wooldid` replay refused until some call completes; the
`macro drop WOOLDID_LASTRC` manual escape and a message-text extension are author options.

## PKG-03 — bundle README (low)

**Ticket:** the install bundle lacks a package-level overview; proposes authoring
`wooldid_README.md`, installing it via the .pkg inventory, and verifying it with a PowerShell
checker.

**Rejected: the premise is an artifact of the colleagues' limited view.** The repo root carries
a maintained README.md — the package-level overview the ticket asks for — where GitHub visitors
and archive readers actually land; installed users have the help file. Bundling a second README
adds a nonstandard installed file and a new version/date sync surface whose main purpose is
managing the problem it introduces (and the checker is wrong-platform with nonexistent paths).
Legitimate residue folded into the release-time PKG commit: refresh the root README's feature
list (permute/fastnose/replay/pretrendpower/gteffects, anticipation per the removal ruling) and
reconcile the Wooldridge citation year with the .pkg's. The ticket's one sound sub-point — if a
README were ever bundled, namespace it `wooldid_README.md` — is noted for the record.

## NFANT-00 — remove the anticipation() feature (high; author-ruled)

**Ticket:** complete removal of `anticipation()` (F4) — parsers, onset shift, diagnostics,
helpers, validator terms, permutation carry, `e(anticipation)`, help, tests — framed as
executing an already-made product decision that the feature is not worth its burden.

**Rejected — author ruling (2026-08-19): keep the feature.** The ticket itself concedes the
feature is numerically correct (verified bit-identical to hand-editing the onset variable on
~100 configurations at rtol 0). The burden case is dominated by the colleagues' own proposed
manifest/replay machinery (NFRPL-03/NFES-03/NFPERM-12), which this project rejected wholesale;
what remains is a small, verified feature matching `did_imputation`'s `shift()` convention that
provides drop-attribution diagnostics the documented manual workaround lacks. Removal would be
a breaking API change with a large coordination bill (roughly ten tickets' proposals assume the
feature exists or doesn't) one month after a six-agent verification wave landed it. The keep
ruling also un-blocks the tickets whose fixes touched anticipation-adjacent code (VAL-G10B's
reorder amended the anticesf caution wording; the band bounds keep their `+anticipation` terms).

## NFGTE-00 — remove the gteffects feature (high; author-ruled)

**Ticket:** complete removal of `gteffects` and `e(gteffects)`, replaced only by a private
cell mapper for NFCOMP-03.

**Rejected — author ruling (2026-08-19): keep the feature; soften the help.** The table is
verified exactly self-aggregating (~1e-15) in its documented domain, is the only inference
surface that survives `fastnose` (a deliberate synergy), and provides csdid-comparable
group-time output; the limitations the ticket cites are disclosed (contreat refusal, Poisson
xb-scale note, no cross-cell VCV by design) or were factually contradicted by the verification
record (the time-varying interactive-controls identity was measured HOLDING at 1.8e-15). The
two honest removal arguments were addressed without deleting the feature: the help's
"building blocks that every aggregate is built from" headline was softened to name the exact
domain, and a partial-omission caveat was added (a dropped coefficient has no row; under
default safety the dropped-coefficient screen errors first). The ~1.5kB of worker headroom the
removal would reclaim stays reclaimable later by subprogramming the extraction block. The
NFCOMP-03 mapper is moot (that ticket was rejected on a refuted premise).

---

*Process note: the NFPERM-02 ticket, initially rejected in triage, was RECLASSIFIED to accepted
after adversarial review refuted the rejection — its frlink-reorder mechanism was empirically
confirmed (stratified exact designs returned p=0 with the observed assignment absent from their
own enumeration) and the keyed-row fix was implemented (5848422/a785318). It is recorded here
for transparency about the process, not as a rejection. Author rulings 2026-08-19 also settled:
the adversarial pre-trend statistics are being FIXED (both tests, new subprogram, esf-scoped) —
so NFADV-01/02/05 are accepted, not rejected; and exact-mode permutation enumerations with
failed assignments keep the publish-with-caution behavior (the plan-of-record refusal reading
was declined; the always-on incompleteness caution and help sentence landed in a785318).*

*Escalations recorded for the author (no action taken): NFPERM-03's markout-domain question —
should deterministically marked-out rows also leave the randomization design (medium lockstep
change with the landed nesting check); NFCOMP-03's caution-demotion option — run the plain-cell
collinearity scan as a Caution under user-typed safety(off) (~5 lines); VAL-G15C's mover-caution
text is imprecise for the within-cell mixing class under safety(off) (one-line text improvement,
author-directive provenance).*
