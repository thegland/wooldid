# Adversarial referee report — attack on the AME-SE synthesis draft (2026-08-24)

Referee: adversarial-referee agent. Target: `~/.claude/jobs/640b9294/tmp/amese/synthesis_draft.md`.
Read: prereg, mc_results, theoryA/theoryB (locked), blind_theorist, code_audit, INV-1, plus
repo cross-checks (`wooldid_standing_issues.md` SI-1..SI-15 + family-A remedy menu,
`wooldid_old_bug_report.md` §5.2-5.4 incl. the NF-1/NF-3 measured rejection tables,
feasibility-verdicts memory E1/F1/F2/F7, `wooldid.ado` L905-1107 lattice/poisson/regress
injection sites). One probe run (refwork/refprobe1.do/.log; prediction written in-file before
running): Part A = mixed-cluster false-negative demo; Part B = byTime equality on a second
panel shape. Second allowed probe held in reserve, not needed.

---

## A. G1 bucket-(iii) call + "bucket (ii) is EMPTY" — CONFIRMED, with a required lens edit

Argued both ways, as instructed:

*For (ii):* the user's own bucket-(ii) parenthetical — "degenerate/one value" — describes G1's
surface signature exactly (one identical SE on all 10 quantities). And SE_cl(phihat) IS a
correct, consistent SE — of phihat. If "correct SE" means "the requested method's correctly
computed output" (Stata ships this CRVE), G1 is (ii): correct-per-method, totally undesirable.

*For (iii):* both locked theories define "correct SE" = consistent for the true sampling
variance of the POSTED estimator. SE_cl(phihat) is a correct SE of a coefficient the output
never reports; attributed to ame_h it understates variance 4-25x. Mis-attribution is an
incorrectness of the posted (number, label) pair, not an undesirable property of a correct
number. And the decisive separator between (ii) and (iii) is (iii)'s own clause: a correct SE
of the labeled estimand IS estimable (cross-cohort dispersion / cluster bootstrap — blind §2,
theory A §1, theory B §1 independently). Both halves of (iii) are satisfied; (ii)'s first half
("correct SEs") is not.

**Ruling: (iii) stands; (ii)-empty stands — but the draft currently earns it inconsistently.**
It uses the estimand lens to empty (ii) while using the within-method lens to float "(iv)
within the requested method" for G1. One lens must govern. Required edit: (1) state the
estimand lens ONCE as the operative definition (it is both theories' locked definition);
(2) keep within-method observations as parenthetical annotations symmetric across (ii) AND
(iv) ("per-method the number is the method's exact, degenerate answer — (ii)-flavored;
per-method no cohort-CRVE ame SE exists — (iv)-flavored; per the estimand lens both collapse
to (iii)"); (3) name the two nearest misses to (ii) explicitly: d4/d5_within_ci main/pre ame
(collapsed value ≈ correct for those two aggregates, ratios 1.00-1.08 — compositional
coincidence, wrong for the other 8 quantities; blind's committed nearest-miss) and d3's
missing SE read against the constant-0 statistic actually computed (theory B's "(ii)-flavored"
note). No group moves to (iv): G3's ame is (iv) only within the CRVE class and in-principle
estimable (audit's own "(iii)-capable in the wider sense"); (iv) is reserved for d3, where the
design — not the VCE — kills every estimator class (confirmed: the L918/L1047 lattice
construction only REFINES cells, so no contreatcontrols/lattice option coarsens the treated
slope; and the injection is shared by the reghdfe, ppmlhdfe, and regress paths).

## B. G5 split verdict — WEAKENED: the split is right, the "(i)" label is a fudge

The REMEDY split is well-supported and I uphold it (blind §4 rules the same split; a t-interval
degrades gracefully, an 8-dim whitened form does not). The BUCKET label is not:

1. The preregistered band does NOT support "(i)": [0.85,1.15] is the committed line for
   "approximately correct"; every G5 case sits mostly in the FLAG band (theory adjudicates)
   and one (d2_bycoh_cu, 0.68) below it. The draft may not cite the band as support — calling
   0.68-0.91 "(i)" is an adjudication OVERRIDE of the prereg's own middle band and must be
   labeled as such if kept.
2. Single-standard problem: the draft empties (ii) via "consistent for the truth" but admits
   G5 into (i) via "textbook behavior of the conventional estimator." Under the first standard
   G5 is (iii): the CR1/HC1 leverage bias is a FIXED proportional bias in the G-grows,
   cell-size-fixed asymptotics (blind: "does NOT vanish as G grows"), and better estimators
   exist (CR2/HC3/jackknife/bootstrap) — (iii)'s definition verbatim.
3. Unflagged 2-vs-1 disagreement: the blind theorist classified ALL these cases (iii)-B —
   his bucket (i) contains only the four d4/d5 cu/rb cases — while theories A/B said
   (i)-with-caveat. The draft silently adopts A/B. Protocol violation by omission (see I).
4. Numbers: the group's true span is 0.68-0.97 (d2_bytime_rb reaches 0.97), not "0.68-0.91";
   under-coverage at the 0.68 edge is ~13pt (coverage ~82%), not "8-12pt" — 8-12pt corresponds
   to the central 0.74-0.80 horizon ratios.

**Required change:** relabel G5 as **(iii)** (leverage-biased plug-in; correct SE possible)
with the remedy split UNCHANGED (displayed SEs: caveat + runtime small-cell note + opt-in
CR2/HC3, no silent default change; joint tests: mandatory fix per F3) — or, if the author's
(i) is deliberately "standard-practice-correct," write that definition down and apply it
grid-wide (which then re-opens (ii) for G1). Add the shape qualifier: the 0.68-0.91 magnitudes
are 6-obs-cell magnitudes; bias shrinks with cell size (see G4 below). Fix the two numbers.

## C. byTime spoofing / "rank guards cannot work" — CONFIRMED (probe-corroborated); stands

The claim never rested on one probe alone: the FWL identity gammahat_gh = s_gh − phihat_t is
algebraic; both locked theories predicted full-rank horizon-varying V from the cross-cluster
locus of phihat_t BEFORE the MC; and 512 MC reps (d1/d2 bytime_ci, ratios 0.38-0.69, WRONG)
independently confirm a full-rank V empty of within-cell noise. My probe Part B re-tested the
exact equality on a DIFFERENT shape and seed (20 cohorts x 9 units x 12 periods, seed
99110001): all 8 reported horizon SEs equal the visible-refit SE_cl(phi_t) to every displayed
digit (.1462556/.0936050/.1108707/.1519299 | .1703972/.1071221/.1586515/.0927140). Not a
draw artifact. The design conclusion (structural predicate, never a rank/eigenvalue guard)
follows from the MC alone and is now doubly probe-anchored. Edit (minor): cite the full
corroboration chain in the G2 row rather than "one blind probe," so the support is transparent.

## D. F1 predicate — WEAKENED on two counts: a demonstrated FALSE NEGATIVE, and the att side

**False negative (probe Part A, prediction pre-registered in the do-file).** d1 data, mixed
cluster variable cmix = cohort for i<=20, unit for i>20 (80 clusters): 10 of 15 treated
cohorts have every cell contained in one cluster, 5 do not. F1-as-drafted ("EVERY treated
(i,relyr) cell inside a single cluster") does NOT fire. Result: full-rank V_es_ame
(eigenvalues .0005-.0245), horizon-varying SEs, ZERO cautions — and the SEs are 0.37-0.77
(mean ~0.50) of the cluster(unitid) SEs on the SAME data; the same-data cluster(i) anchor
reproduces the rank-1 collapse (.0397219 x8, eig {.0126, ~1e-16}). The true-null run prints
p = 1.2e-08 on one horizon. So partial containment silently reproduces most of the
annihilation while defeating the predicate. Realistic analogues: cluster variables that
cross-cut cohorts (region/industry), multi-way clustering with one containing dimension.
**Required change:** predicate fires on ANY treated cell contained in a single cluster
(report the contained fraction; full containment → refuse, partial → strong caution at
minimum). Note this phrasing is inherited from §5.2 remedy (b) — the standing-issues menu
needs the same one-word fix.

**False positives: hunted, none found.** Every contreat run saturates the treated cells by
construction (the L1099 injection is shared by the reghdfe, ppmlhdfe, and regress paths, and
weighted estimation preserves the annihilation — the weighted normal equations zero the same
weighted within-cell moments the weighted cluster scores are built from), so whenever the
containment condition holds, the score-annihilation is guaranteed and the plug-in CRVE ame
SEs are never fine. The only near-miss class is d3-type data: the predicate fires CORRECTLY
(refusal is still right) but F1's finer-clustering advice would be wrong there — which is the
ordering fix above, a message-accuracy issue, not a false positive.

Also under-specified, must be stated in F1: (a) scope across estimation paths — the same
injection feeds reghdfe, ppmlhdfe (poisson; ML within-cell moment conditions annihilate the
same way), and the regress/`unconditionalse`/`customvce` path (INV-1 probe5: rank-1 under
plain regress too); (b) ordering — the d3/G7 unidentified-refusal must run BEFORE F1, else
F1's "use finer clustering" message gives wrong advice on d3-type data; (c) blast radius —
any battery/oldnew config with contreat + containing cluster changes displayed output;
F1 currently carries no re-bless line.

**Att side — the draft over-claims "VALID att variances" and under-scopes by silence.**
The only evidence for the G3/G9 att blocks is that plain regress posts positive,
normal-looking SEs (.152-.214). But the repo's own record contradicts "valid": §5.2 measured
att-side diagonal variance understatement of 1.5-2.3x on multi-unit designs under cohort
clustering (SE ratios ~0.66-0.82 — the same FLAG band) with contrast directions 7-20x off,
and INV-1 found V_es_att rank 5-of-8 under cluster(i) — "not a separate sound path." This
campaign's MC recorded NO att quantities, so att correctness is assumed, not measured — and
an att-side MC (`mc_amese_att.do`, "so the att-side SE correctness under each VCE is measured
rather than assumed") is RUNNING right now. **Required change:** delete "valid" (say
"positive, per-method plug-in att variances, themselves family-A-affected per §5.2/INV-1");
state F1's att scope explicitly (at minimum: the F1 caution also names V_es_att joint
surfaces); hold the final att scoping open until mc_amese_att lands and cite it.

## E. G7 (iv) for d3 — CONFIRMED; refusal is right, not over-paternalistic; two edits

"Genuinely unidentified" is strictly true for the labeled estimand within wooldid's design
space: dose constant within cohort ⇒ constant within every cell ⇒ no within-cell contrast
exists under ANY contreatcontrols, VCE, or lattice (lattice only refines cells — verified at
L918/L1047; the drop is X-structural and VCE-invariant, audit F). A pooled cross-cohort
dose-response IS conceivable from d3-type data but only outside this machinery under strictly
stronger assumptions — that makes it steer material, not a bucket change. The named refusal is
the opposite of paternalism given the current output: hard 0s sit in e(b) indistinguishable
from estimates, behind a generic caution — a live hazard for pipelines. Edits: (1) refusal
scope = the ame block only (att side stays; no abort); (2) the message should pair refusal
with the misspecification steer ("dose has no within-cell variation; ame not identified;
cohort-level dose designs need a different estimator class / stronger assumptions"), matching
the author's say-not-estimable + steer response types.

## F. d1c "never steer to robust" — WEAKENED: make the steer conditional, and flag the dissent

Steelman upheld in part. The blind theorist ruled the OPPOSITE of the draft ("the correct VCE
advice for d1c is robust (or unit clustering) — or a design-based cross-cohort estimator"),
and the absorption of cell-level shocks by the saturated cell dummies is STRUCTURAL, not a
lucky draw: any shock constant within (g,t) cells is conditioned out of every treated-cell
slope by construction (blind's table proof: est-SD(d1c)/est-SD(d1) ≈ 0.60-0.66 = sigma_idio
across all 10 quantities). Robust's 0.75-1.05 is therefore not an accident — for CELL-level
dependence. What robust cannot cover: unit-serial or other unabsorbed dependence, and the
control-side phihat channel (shocks in control cells are NOT absorbed; diluted here by ~1080
control obs — a shape/sigma_nu-fragile grace, see G4). **Ruling:** the PACKAGE default remedy
stands (bootstrap/dispersion — the only estimator valid under the user's own stated dependence
model), but "NEVER a steer to robust" must become a CONDITIONAL steer: the caution/docs should
say "if the dependence you fear is cell-level, the design absorbs it for the ame and
robust/unit clustering is approximately valid (small-cell caveat applies); if unit-serial or
broader, only the cluster-level replacement estimator is valid." This also reconciles the
draft with theory A's locked remedy for G1 ("steer to robust/unit-level cluster when iid
plausible") and with the empirical record it must not bury (robust 0.75-1.05 vs ci 0.18-0.46).
Required: flag the blind-theorist dissent and this resolution in the draft (see I).

## G. Methodology threats — one attribution error to fix; otherwise survivable with qualifiers

1. **Harness nondeterminism / "distributionally harmless": CONFIRMED with a stated premise.**
   The exchangeability argument is airtight IF the perturbation is a within-cell permutation
   of iid unit-level components — which is what the evidence shows (att rows bitwise invariant
   in the one dissected divergence, p_d1c.log; fresh-session builds bitwise reproducible,
   probe_cef). d1c is the EASY case, not the hard one: its shock is cell-constant, hence
   permutation-invariant by construction; all unit-level components (dose draws, noise) are
   iid within cells in every DGP. Edit: state the premise + evidence explicitly; add the cheap
   guard (assert att-row invariance in future harness runs) and the latent-7 `cap confirm
   matrix` fix to the suite-side ticket list.
2. **MC noise: no group verdict flips**, but two labels need honesty notes: d2_bycoh_cu's
   "WRONG" rests on a single 0.68 against the prereg's own classify-on-the-pattern rule
   (9 of 10 quantities are FLAG-band; should read FLAG/band-edge — the draft already treats it
   so, and should say mc_results' label deviates); d2_bytime_ci's top ratio 0.69 abuts the
   band. **Attribution error (required fix):** the G5 row's "joint tests over-reject 20-35%"
   sits in the "MC evidence" column, but THIS campaign measured no rejection rates — the
   numbers are NF-1 (M=1000, cluster-saturated designs, 20.5-48.4%) and NF-3 (M=504, A4 =
   the same 30x6 shape, robust: adversarial-ame 31.5%, joint_Choriz 22.0%). Cite §5.2/§5.4;
   note cluster(unitid) joint rates are inferred, not measured.
3. **Band circularity: minor and disclosed-by-design, but one use is illegitimate.** Bands
   were committed post-structure/pre-MC; the FLAG band explicitly delegates to theory — a
   deferral device, fine. Whether band edges were placed with knowledge of the theorists'
   predicted 0.65-0.85 range is unverifiable from the record; add one disclosure line. What
   the draft may NOT do is cite the FLAG band as supporting "(i)" for G5 (it says the
   opposite — see B).
4. **Shape fragility — add qualifiers to:** G5/G6 magnitudes and the under-coverage figures
   (6-obs cells; leverage bias ~2p/n_cell shrinks toward 1 with more units per cell — with
   realistic cell sizes the "material caveat" largely dissolves, which is exactly why F5's
   runtime note must key on measured cell size, as drafted); the joint-test rates (dimension
   x cell size); G4/G9 blackout FREQUENCIES (numerical knife-edge rates are shape- and
   noise-specific); d1c's "robust nearly right" (sigma_nu- and control-share-dependent via the
   unabsorbed control-side channel). Shape-ROBUST (algebraic, no qualifier needed): G1/G2/G3
   annihilation signatures, within≡basic, d3 unidentification — and my Part B replicates the
   G2 signature on a second shape.

## H. Fix list — order largely right; five required repairs

1. **ID collision (mechanical but mandatory):** the draft's fix IDs F1/F2/F7 collide with the
   repo's established feature IDs (feasibility campaign: F1 = boottest, F2 = the SHIPPED
   e(b_es)/e(V_es) capture, F7 = pretrendpower) — and the draft itself uses "F2" in BOTH
   senses in one document (fix list vs auditor-findings section). Rename fixes (e.g.
   AME-1..AME-7 or SI-4a..) before this goes anywhere near the author.
2. **F6 resurrection risk: REAL — gate it.** The regress fallback on G3/G9/G4 data restores
   the regime-A rank-1 ame SEs that the blackout was accidentally suppressing (audit regimes
   B/C/E: the regress discriminator posts .0324-type collapsed ame SEs next to the positive
   att SEs). F6 must be explicitly subordinate to F1: restore the att block only; whatever V
   the fallback produces, the ame side stays under F1's suppression. Also: "repair/split the
   V" is impossible post-zeroing (the information is destroyed — POSTEXP) — the wording must
   say re-estimate, not repair. And e(V)==0 detection alone is insufficient: margins carries
   its own copy of the screen and blanks whole se rows on healthy e(V) (audit latent 5, d4
   regress path) — detection must also check r(table)/r(V) after the margins calls.
3. **F3 scope:** yes, extend to pretrendpower under the predicate — roth_prep/nissize are
   already listed (SI-3), and the mdt columns, though their FORMULA validated (SI-3 "need no
   change"), consume Sigma_pre: garbage-in under the predicate. §5.2 remedy (b) already names
   pretrendpower + the honestdid handoff in the caution list; F3/F1 should match it. No
   contradiction with SI-3 once phrased as input-validity, not formula-validity.
4. **Missing items:** un-swallow the reghdfe warnings around L1099 (audit latent 2 — capture
   + re-display; one line, high diagnostic value); doc-grade notes for latent 8 (e(V) =
   diag(se^2) with fabricated zero off-diagonals — cross-estimate lincom/test on e() silently
   wrong in ALL regimes; fold into F7 with "F2-matrices are the sanctioned joint surface")
   and latent 9 (the -1e100 sentinel reaching raw consumers); suite-side tickets for latents
   6/7. F6 needs a priority tag (its detection+message half is P2-level user harm: G3/G9/G4
   users currently lose ALL inference with an unexplained caution; the att-restoration half
   waits on mc_amese_att per D).
5. **Prereg output-contract gap:** the prereg promises fixes "mapped to the author's three
   response types ... each with its re-bless blast radius." The draft maps neither
   systematically (blast radius given only for F2/F3; response-type labels absent). Add both
   columns. Consistency with established rulings otherwise good: F1's predicate = §5.2
   remedy (b) (with my "any" fix), F2 = remedy (c), F5/F7 = remedies (a)/(b) — but F2 must
   cite the boottest feasibility verdict (wild bootstrap via boottest is DEAD on the stored
   model — multi-FE, continuous excluded; so: pairs cluster re-estimation bootstrap or the
   analytic IM/dispersion estimator as primary, with a cost note).

## I. Internal inconsistencies — the structural one is the missing theory scorecard

1. **No adjudication of the locked theories anywhere in the draft.** The prereg locked A/B
   precisely so the synthesis could score them; the draft never does — which is exactly how
   the two material disagreements survived unflagged: (a) blind (iii)-B vs A/B (i)-with-caveat
   on G5 (draft silently picks A/B); (b) blind "steer to robust" vs B "never robust" on d1c
   (draft silently picks B). Required: a short scorecard (A/B's discriminating predictions P1-
   P6/rows: overwhelmingly confirmed; A's "d1c SE ≈ unchanged" vs B's "strictly larger" both
   effectively unresolved by the ratio-only grid; A's "byCohort reproduces every rep" upheld
   at n=256 where the audit's knife-edge framing predicted near-determinism) + explicit
   resolution notes for (a) and (b).
2. Number/label fixes: G5 span 0.68-0.97 not 0.68-0.91; G9 blackout rate "~13% of MC reps
   (33/256); 1/16 in the audit's fresh-session sweep" — do not present "1-in-16" as THE rate
   when an n=256 estimate exists; G4's "95/256" (37%) vs "3 of 4 audited" (75%) — compatible
   at n=4 but flag the small-n tension; G1's ratio line (~0.2/0.4/0.5) describes d1/d2/d1c
   only — d4-posting-reps and d5_within_ci in the same row have the compositional profile
   (main/pre ≈ 1.0, horizons 0.6-0.75) and must be footnoted or moved; d1c_basic_ci appears
   in both G1's case list and G9 (say the double-listing is deliberate); "Old wooldid shares
   every pathology" over-claims — verified for the rank-1 V + the g-inverse p (INV-1 probe8,
   SI-4); the blackout/leverage regimes are same-construction inferences, label them so.
3. Under-coverage "8-12pt" vs the band edge (~13pt at 0.68) — see B.
4. mc_results' d2_bycoh_cu "WRONG" vs the prereg's pattern rule — see G.2; the draft should
   note the deviation rather than inherit the label silently.

---

## Probe documentation (refwork/)

- refprobe1.do/.log — Part A: d1 seed 99000001, cmix = i(<=20)/100+unitid(>20), wooldid
  basic x {cluster(cmix), cluster(unitid), cluster(i)} on identical data (build-once,
  save/reload). Results: cmix full-rank (eig .0005-.0245), NO caution, SE ratios vs cu
  0.37-0.77, spurious p=1.2e-08 on a true null; ci anchor rank-1 (.0397219 x8). Part B:
  20x9x12 d1-style panel, seed 99110001, byTime+cluster(i); all 8 horizon SEs = visible-refit
  SE_cl(phi_t) to every displayed digit. refA.dta/refB.dta frozen.

## The three most consequential required changes

1. **Att-side re-scope (D):** delete "VALID att variances"; the att block under cluster(i) is
   family-A-affected per the repo's own §5.2 (1.5-2.3x) and INV-1 (rank 5/8); gate F6 to
   att-only restoration under the family-A caution, subordinate to F1's ame suppression
   (else it resurrects the rank-1 SEs the blackout was suppressing); hold final att scoping
   for the in-flight mc_amese_att results.
2. **F1 predicate "every" → "any" + scope/ordering (D):** probe-demonstrated silent false
   negative at partial containment (full-rank V, 0.37-0.77 ratios, no caution); state
   estimation-path scope (reghdfe/ppmlhdfe/regress), d3-check-before-F1 ordering, margins-
   level r(table) detection, and F1's blast radius; propagate the one-word fix to §5.2
   remedy (b).
3. **Single correctness standard + theory scorecard (A/B/I):** fix the estimand lens once,
   relabel G5 as (iii) with the remedy split unchanged (or write down the override), and add
   the scorecard that surfaces the two silently-resolved disagreements (G5 bucket; d1c robust
   steer — resolve the latter as a CONDITIONAL steer per F). Plus mechanical: rename the
   fix IDs (F1/F2/F7 collide with shipped feature IDs, "F2" used in two senses in one
   document) and re-attribute the 20-35% joint-test rates to NF-1/NF-3 measurements.
