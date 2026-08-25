# Code-path audit — AME standard-error regimes A-F (2026-08-24)

Auditor: code-path auditor agent (multi-agent AME-SE diagnosis). Line cites are into
/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado (bugfixes branch,
v1.4 2026-08-19) unless prefixed reghdfe5.mata (= ~/Library/Application Support/Stata/ado/
plus/r/reghdfe5.mata, reghdfe 6.12.3 08aug2023; Stata 19.0 MP). Probes + frozen data live in
this directory (probe_bd/probe_cef/probe_fix/probe_fix2/probe_fix3/probe_fix4 .do/.log;
d1_r1.dta, d1c_black.dta, d1c_clean.dta, d3_r1.dta, d4_r1.dta, d5_r1.dta). Repo untouched.

## Headline findings

1. **wooldid never wipes or fabricates an SE on these paths.** Every displayed/posted SE is
   a verbatim relay of margins r(table); every caution is triggered by missingness alone.
2. **The B/C/E "total 20-of-20 blackouts" are a three-layer event, and the decisive layer
   is NOT wooldid and NOT the CRVE math:** reghdfe computes a raw clustered V whose
   degenerate directions carry roundoff; **Stata's built-in `ereturn post` screen then
   declares it "variance matrix is nonsymmetric or highly singular", stores e(V)
   identically ZERO (bitwise, all 334x334 entries), and displays every SE as missing.**
   A controlled experiment (POSTEXP, 8 cases) brackets the screen's trigger class:
   asymmetry beyond a tolerance in (1e-12, 1e-6] or indefiniteness carried by
   off-diagonals ("correlations" > 1) — while PSD rank-deficiency, exact-zero blocks,
   and pure-diagonal negativity all pass.
   Plain `regress, vce(cluster i)` + the identical margins calls on the SAME data post
   positive, normal-looking ATT SEs (and the usual rank-1 ame SEs) — so the raw plug-in
   CRVE for the ATT rows exists and is positive; it is destroyed in transit.
3. wooldid runs reghdfe under `qui` (L1099), so the only mechanistic warnings a user could
   see ("variance matrix is nonsymmetric or highly singular", "missing F statistic;
   dropped variables due to collinearity or too few clusters") are swallowed; the user
   gets only the downstream "standard errors were not identified for 20 of 20" caution.
4. Structure-pass fact C is mis-stated: the d1c blackout is NOT a deterministic property
   of the cohort-shock error process. It is a numerical knife-edge hit by 1 of 16
   realizations (the seed-60000001 fresh-session realization among them). 15/16 d1c
   realizations behave exactly like regime A (identical positive ame SEs); d1 never
   blacked out in 16 seeds. Regime E is also a knife-edge, at high frequency (3 of 4
   realizations trip). Regime B (byCohort+ci) tripped in 4/4 observed realizations —
   effectively deterministic (its raw slope block is an exact 165-dim null space, so
   the relative asymmetry there is always O(1)).
5. Regime D's posted e(V_es_ame) is a bit-exact faithful J·e(V)·J′ of the model CRVE
   (max abs deviation 0.0): the code computes what the method says; any defect there is
   the method's, not the code's.

## 0. The missing-SE machinery map (task 1)

Pipeline from model to symptom, with the exact decision predicates:

1. **Model.** L1099 `reghdfe y i._i#i._t#i.targetflag##c.dose ..., absorb(i.i i.t
   <contreat controls>) <vce> tol(1e-9)` run **quietly**. contreat controls are engineered
   into absorb(): basic = `i.null##c.(ringer)` (L1017-1023), byCohort = `i.i#c.(ringer)`
   (L1030), byTime = `i.t#c.(ringer)` (L1034); ringer = double copy of dose (L906),
   invisible to margins. df_r captured L1128; threaded to margins as `df()` via
   `marginsdof` (L1609).
2. **Inside reghdfe (cluster case).** reghdfe5.mata `reghdfe_cluster` builds
   V = D·M·D·q, D = invsym of the partialled X'X, M = Σ_g s_g s_g', q =
   (N−1)/(N−nested−K)·G/(G−1) (L2896-2911; q = 1.249 in regime B — healthy). Its
   missing-V branch (L2720-2724: `if (missing(V)) V = J(.,.,0)`, message only at
   verbose>0) did NOT fire in our traces. reghdfe_fix_psd runs only for multi-way
   clustering (L2909-2912) — not in play.
3. **Stata's `ereturn post` screen.** reghdfe's ado posts b/V via plain `ereturn post`
   (reghdfe5.ado L457, RegressOLS). In regimes B/E and the d1c-blackout realization,
   Stata's built-in emits "warning: variance matrix is nonsymmetric or highly singular."
   and the STORED e(V) comes back identically zero; the estimation table shows all SEs
   missing. (Same warning semantics documented in Stata whatsnew13/16/18 entries:
   warning + "missing values for all the standard errors".) A controlled mini-experiment
   (probe_fix4 POSTEXP, eclass `ereturn post` of hand-built 3x3 Vs) pins the trigger:
   nonsymmetry beyond a tolerance in (1e-12, 1e-6] — see the POSTEXP section.
4. **Margins calls.** Binary att: L1657 (main), L1709 (pre), L1771 (ES). Continuous ame:
   L1915, L1962, L2186. All `[iw=w_att]`, w_att=1 here (L1209). With the zeroed e(V),
   margins computes r(V) = J·0·J′ = 0 exactly (verified: max|r(V) − J e(V) J′| = 0) and
   reports se = "." (missing) for the affected quantities in r(table). Margins also runs
   its own copy of the degeneracy screen: on a HEALTHY posted e(V) whose delta-method
   r(V) is near-degenerate it prints the same "nonsymmetric or highly singular" warning
   and blanks the ENTIRE r(table) se row while r(V) keeps its (partly positive) values
   (observed on the d4 regress path, regime E). So "e(V)=0" — or a knife-edge r(V) —
   surfaces as "all SEs missing", never as "SEs 0".
5. **Relay.** r(table) rows 1-6 go verbatim into the `___C*` matrices (e.g. L1916-1921,
   L2187-2192) -> `wdidpooled` (L2820-2898) -> display (wooldid_display L5690) and
   e(wdidmainresults) (L3523). No wooldid-side transformation of the se column. (The only
   code that actively sets se to missing is the impvalues IM-degeneracy guard L1640-1645/
   L1840-1858 — not on this path.)
6. **e(b)/e(V) posting + caution predicate.** L3279-3280 builds e(V) = diag(se²) from the
   displayed-SE column (off-diagonal covariances are fabricated zeros by design,
   L3266-3268). L3291 `ereturn post` rejects any missing diagonal; the rescue logic
   L3293-3327 counts missing diagonals (`evpostmiss`, L3307-3309):
   - some good & some missing → missing cells replaced by sentinel −1e100, repost, caution
     "...not identified for X of Y ... e(V) contains placeholder entries" (L3311-3319).
     [Regime F's "10 of 20".]
   - ALL missing (`evpostgood`==0 → L3311 gate fails) → b-only post (L3322) + caution
     "...not identified for X of Y ..., so no e(V) was posted" (L3324). [B/C/E "20 of 20".]
   **The predicate is exactly `missing(diag)` — nothing tests zero, identical, tiny, or
   rank-deficient variances.**
7. **F2 capture & suppression (the e(V_es_ame) path).** Capture at L2194-2207 (ame) and
   L1780-1795 (att; grabs the second half of r(b) via `matrixstartgrab` because
   dydx(targetflag) returns base+level columns) is gated ONLY on `"marginsdof" != "nose"`
   — raw r(b)/r(V) captured even when all-zero. Posting gate L3352-3368: `f2ok` = at least
   one diagonal of the captured r(V) is >0 and <. . **e(b_es_ame), e(V_es_ame), and
   e(numpre_ame) are all inside `if f2ok` (L3359-3363), so blackout regimes return
   NEITHER b nor V** (empirically confirmed: B run has b_es_ame/V_es_ame/b_es_att/V_es_att
   all ABSENT while e(b) exists 1x20). Partially-degenerate families post with the L3365
   caution + e(es_ndegen_*).

## A. basic + cluster(i) (d1/d2/d5-within): silent rank-1 collapse

- Reconfirmed on the frozen file d1_r1.dta (probe_fix PART 4): model e(V) has 331
  positive diagonals; the 165 live treated-cell slope diagonals are ALL equal to
  .000961874838 (= Var_cl(phihat)) to 11 digits; margins ame r(V) diag = the same
  constant x8 (identical SEs .0310141071); att r(V) diag positive and horizon-varying
  (.0224-.0442). e(V_es_ame) posts silently: 8x8, eigenvalues {.0076949987, ~8.7e-17...}
  — exact rank 1. Mechanism as INV-1: each treated (cohort x rel-period) cell is one
  calendar period inside one cluster and carries its own visible dummy+slope; the CRVE
  score of every cell slope is annihilated (only the pooled absorbed-ringer moment
  survives), so the slope block of e(V) = Var_cl(phihat)·11′.
- Faithful-relay verdict: the delta-method chain is exact; the defect is the reported
  covariance object (CRVE cell-saturation), reported with **no caution anywhere** because
  every diagnostic tests only missingness (see latent issue 1).
- Bucket input: the ingredients for a defensible SE are all in hand at L2186: per-cell
  slope estimates (in e(b)), horizon weights w_gh (rows of r(Jacobian)), within-cell
  residual variances (from the stored model ___wooldidfullmodel, L1117), Var_cl(phihat)
  (the constant itself). Var(b_h) = Σ_g w_gh² Var(s_gh) + Var_cl(phihat) is assemblable;
  a cohort-cluster bootstrap of the L1099+L2186 functional is the design-based
  alternative. So: correct-SE-possible; code posts a self-refuting V silently.
- Fix sites: detection after L1915/L1962/L2186 (rank/dispersion test on r(V));
  F2 gate L3352-3368 (extend beyond missingness); display caution near L3318;
  documentation steer at L1030-1034 / help file (cluster at a level finer than the
  contreat cells, e.g. unitid — struct pass shows cu/rb configs are full-rank).

## B. byCohort + cluster(i) (d1/d2): total 20-of-20 blackout — DETERMINISTIC, and the
## ATT half of the blackout is a computation-layer artifact

- wooldid surface (probe_bd, frozen d1_r1): rc 0; table shows all 20 point estimates with
  se "."; L3324 caution; e(b) 1x20 posted; e(V) ABSENT; e(b_es_ame)/e(V_es_ame)/
  e(b_es_att)/e(V_es_att) ALL ABSENT (f2ok=0).
- Model level (probe_bd/probe_fix, manual `reghdfe y i._i#i._t#i.targetflag##c.dose,
  absorb(i.i i.t i.i#c.(ringer)) cluster(i)`): **e(V) is identically zero — all 334
  diagonals and all 111,556 entries bitwise 0, no missing entries.** df_r=29, df_m=330,
  df_a=41, only 2 omitted coefficients (o.dose and the base-cell slope). reghdfe's own
  display shows every SE missing and prints (swallowed by wooldid's qui):
  "note: dose is probably collinear with the fixed effects", "warning: missing F
  statistic; dropped variables due to collinearity or too few clusters", "warning:
  variance matrix is nonsymmetric or highly singular." Its DoF table marks i (30/30
  redundant — nested in cluster) and i#c.ringer (30 coefs, '?').
- verbose(1) trace (probe_fix2): reghdfe's internal q = (2160−1)/(2160−372)·30/29 =
  1.2491 — healthy; the reghdfe-side missing-V zeroing branch (reghdfe5.mata L2720-2724)
  did NOT print. ⇒ the zero is created at **Stata's `ereturn post` screen** (layer 3 in
  the machinery map).
- Raw-CRVE discriminator (probe_fix PART 2, plain `regress y i.i i.t i.i#c.ringer
  i._i#..##c.dose, vce(cluster i)` on the same file): live slope-block diagonals
  1.7e-28..9.9e-22 (the structural zero of the byCohort design, realized as roundoff —
  with the cohort-specific baseline slope absorbed, every within-cohort calendar period
  is a saturated (1,dose) piece, so every cluster score vanishes in the slope
  directions); **live dummy-block diagonals positive-real up to .0456; margins att SEs
  .152-.214 post normally; ame SEs ~6.4e-12.** So the plug-in CRVE for the ATT rows is
  positive and well-defined on this data; only the ame rows are structurally degenerate.
- Determinism: 4/4 observed realizations (d1_r1, seed 60000002, struct-pass d1 and d2)
  show the identical 334/334-zero + full blackout. Effectively structural: the raw slope
  block is an exact 165-dimensional null space (diag ~1e-22 under regress), so the
  relative asymmetry of reghdfe's D·M·D there is always O(1) and the POSTEXP nonsymmetry
  screen trips every time.
- **Verdict:** two defects stacked. (i) AME rows: genuine CRVE annihilation — worse than
  regime A (the pooled-ringer channel is also nested in the cluster, so even the rank-1
  constant dies; the raw ame variance is exactly 0). Given cohort-level plug-in CRVE, no
  ame SE exists — the missing SEs are, accidentally, the honest answer. (ii) ATT rows
  (and ES_att): the raw CRVE is positive; the blackout of these 10 estimates is
  collateral damage from Stata's post screen zeroing the whole matrix. wooldid faithfully
  relays the wreckage and never learns why.
- Task-2 hypothesis adjudication: (a) "reghdfe itself returns missing/zero V for
  everything" — TRUE in effect (e(V) ≡ 0 as stored), but NOT via cluster-dof exhaustion
  (df_r=29, q=1.249, both healthy); the zero is stamped at Stata's ereturn-post screen,
  not computed. (b) "margins propagates a singular V into all-missing" — TRUE
  mechanically: J·0·J′ = 0 exactly and margins renders zero variance as se=".".
  (c) "a wooldid-side reaction over-broadly wipes SEs" — FALSE: no wooldid code wipes
  SEs on this path; wooldid's output is a faithful relay of (a)+(b).
- Latent observation (normalization dependence): the regress and reghdfe parameterizations
  drop different collinear columns, and the resulting AME POINT estimates differ (.554 vs
  .723 at the first horizon on identical data) while ATT points are identical — under
  byCohort+contreat the displayed ame is not invariant to the collinearity-drop choice.
- Bucket input: ATT rows — bucket (iii) shape: a correct(-per-method) SE exists (regress
  shows it) and nothing is posted. AME rows — no cohort-clustered plug-in SE exists
  (bucket (iv) flavor for this VCE), but a correct SE is computable from other
  information the code holds (cross-cohort dispersion of the per-cohort cell slopes
  gamma_gh, which are in e(b), around the horizon mean; or cluster bootstrap), so
  (iii)-capable in the wider sense.
- Fix sites: post-L1099 detection (one line: e(V) all-zero / mreldif(e(V),0)==0 ⇒
  mechanistic error message naming the byCohort+cluster interaction); un-swallow the
  reghdfe warnings around L1098-1099 (capture + re-display); wooldid's own reg-based
  path (`unconditionalse`/`customvce`, L1102-1108) already builds the same model through
  plain regress and is the natural internal fallback for the ATT rows; margins sites
  L1657/1709/1771/1915/1962/2186 unchanged (they are faithful); documentation for
  contreatcontrols(byCohort) + cluster(cohort) as a refused combination.

## C. d1c basic + cluster(i): blackout is a REALIZATION-SPECIFIC numerical knife-edge,
## not an error-process property

- **Correction to the structure-pass fact:** across 16 fresh d1c realizations
  (probe_cef sweep), exactly 1 — the seed-60000001 fresh-session realization, saved as
  d1c_black.dta — blacks out; the other 15 behave exactly like regime A (all ame SEs
  identical positive; min ame r(V) diag .0004-.0037). The parallel d1 sweep (same seeds,
  iid errors): 0/16 blackouts. Also: the DGP's y-dose pairing is sort-RNG-dependent
  (amese_cases.doh bysort tie-breaks), so "same seed" runs are only reproducible from a
  fresh session — the earlier p_d1c.log ran the same nominal command/seed after another
  wooldid call and got a CLEAN regime-A run with different ame points but identical att
  points (cell-mean invariance) — see latent issue 6.
- On d1c_black.dta the manual basic+ci reghdfe reproduces the B signature exactly:
  e(V) ≡ 0 bitwise (334/334), same 2 omitted terms, same swallowed Stata warning, margins
  r(V) = 0, all r(table) SEs missing (probe_fix PART 6).
- Raw-CRVE discriminator on the SAME data (plain regress + identical margins): AME point
  estimates IDENTICAL to the reghdfe path (-.0704, .3692, ...), ame SEs all .0324 =
  sqrt(Var_cl(phihat) = .0010517) — i.e., the ordinary regime-A rank-1 collapse — and
  positive att SEs (.399, .300, ...). **The raw plug-in CRVE on the blackout realization
  is unremarkable; the entire 20-of-20 blackout is the Stata-post-screen wipe of
  reghdfe's numerically-degenerate V.** Why this realization trips the screen and the
  other 15 don't is a property of roundoff in the degenerate directions of reghdfe's
  D·M·D product (y-dependent through the meat; the design matrix is identical between
  d1-s0 and d1c-s0, so the trigger is not the dose draw); POSTEXP (below) shows the
  screen's predicate is nonsymmetry beyond a tolerance in (1e-12, 1e-6].
- Contrast fingerprints (probe_fix3): clean d1c realization (seed 60000003): live slope
  diagonals all = .0024950659 = Var_cl(phihat) from the visible-ringer refit
  (.00249506593923219) — the exact regime-A signature. Blackout realization: the
  visible-ringer refit is ITSELF wiped (all 335 diagonals zero), so its
  "Var_cl(phihat)=0" is the wipe again; the true plug-in (from regress) is .0010517.
  verbose(1) on d1c_black (probe_fix3): reghdfe q = (2160−1)/(2160−343)·30/29 = 1.2292
  — healthy; NO reghdfe-mata zeroing message; Stata's "variance matrix is nonsymmetric
  or highly singular." fired. Dose forensics on d1c_black: min within-cell dose sd .139,
  no cohort-constant-dose cohorts — nothing special in X (and X is bitwise identical to
  the clean d1-s0 realization, since the d1/d1c generators share the same first 210 RNG
  draws), so the screen trip is driven by y through the meat's roundoff.
- Verdict: faithful relay by wooldid of an upstream (reghdfe→Stata post) numerical wipe;
  contreat-specificity of the observed pattern (binary-att-only wooldid on the same data
  posts fine) is because the binary-only model has no saturated-cell dose columns, so
  its V is far from the screen's knife edge.
- Bucket input: what the blackout HIDES is regime A (rank-1). So: ame rows — same as
  regime A's bucket input (a correct SE is possible via cross-cohort dispersion /
  bootstrap; what would be posted absent the wipe is the self-refuting rank-1 V);
  att rows — bucket (iii): positive plug-in SEs exist (regress) and are lost.
- Fix sites: same as B (post-L1099 zero-V detection + warning surfacing + regress
  fallback). Note for the MC campaign: d1c_basic_ci results must be read as a mixture:
  ~15/16 of reps are regime-A reps; ~1/16 blackout reps show up as rc-0 all-missing rows.

## D. byTime + cluster(i) (d1/d2): full-rank posted V — faithful to the bit

- probe_fix PART 5 on d1_r1: max|r(V) − J·e(V)·J′| = 0.0; max|e(V_es_ame) − r(V)| = 0.0;
  max|e(b_es_ame) − r(b)| = 0.0. Eigenvalues .0017-.0421 — full rank 8.
- Verdict: **the code computes exactly what the method says.** The byTime absorb
  (i.t#c.ringer) restores a cross-cluster moment channel (period-specific baseline dose
  slopes are shared across cohorts), so per-horizon variation reappears; but the
  own-cell slope noise is still annihilated by cohort-level CRVE (cells remain
  cluster-saturated), which is why these SEs sit below the cluster(unitid)/robust
  counterparts (struct: e.g. post3 .1479 vs .2071/.2098). Whether they are CORRECT is
  the MC's question; the audit's answer is: any defect is the estimator's, not the
  plumbing's.
- Bucket input: if MC says the SEs understate, this is bucket (iii)-by-method: the same
  correct-SE ingredients as regime A exist (cells, Jacobian, within-cell variances);
  fix sites identical (detection can't key on rank here — a structure test would need
  the cell-in-cluster saturation predicate: cluster() variable == the cohort argument
  while contreat cells are cohort x period).

## E. d4_within + cluster(i): total blackout (all-treated staggered)

- Surface (struct2.log L6690-6739; reconfirmed probe_fix3 on this session's frozen
  d4_r1.dta): wooldid rc 0; L3324 caution; e(b) 1x20 posted; e(V) ABSENT;
  e(b_es_ame)/e(V_es_ame)/e(b_es_att)/e(V_es_att) ALL ABSENT — the same f2ok predicate
  as B (task 7 confirmed for E).
- Model level (probe_fix3, manual basic-within reghdfe on d4_r1): **e(V) ≡ 0 bitwise
  (664/664 diagonals, no missing entries)** — the B signature at larger K. df_r=29,
  df_m=649, df_a=12, 13 omitted terms: ALL 11 cell dummies of cohort 30 (reghdfe's
  normalization choice in the no-never-treated design) + o.dose + the base-cell slope.
- Margins level (probe_fix4): r(V) = 0 exactly, all r(table) SEs missing — same relay
  as B. NOT strictly deterministic: extra-seed check gives s=60000011 CLEAN (evzero
  14/664, all ame SEs posted) and s=60000012 blacked out — so E is a high-frequency
  knife-edge (3 of 4 observed realizations trip; B remains 4/4).
- Raw-CRVE discriminator (plain regress on d4_r1): the raw plug-in is the REGIME-A
  picture — ame slope-block diagonals all = .0429000244 (rank-1; margins ame SEs .2071
  identical, b identical to the reghdfe path), dummy block positive. So absent the wipe,
  d4_within_ci would post regime-A rank-1 ame SEs, not a blackout. BUT: one ES_att
  horizon's raw variance is genuinely near-zero (7.86e-10 at _t=11), and on this regress
  path **margins ITSELF printed "variance matrix is nonsymmetric or highly singular."
  and blanked the entire att se row while r(V) kept positive entries** — the same Stata
  screen exists inside margins' delta-method, and the all-treated design carries genuine
  near-degeneracy in some dummy directions that d1-type designs don't.
- d5 same-seed contrast (d5_r1.dta): pure regime A — live slope diagonals all =
  .0324059659 = Var_cl(phihat); margins ame SEs .1800165712 x8 (matches the struct-pass
  displayed .180016571 exactly); no screen trip. The d4-vs-d5 difference (time-varying
  vs unit-constant dose) moves the design just enough to flip the numerical knife-edge.
- Verdict: faithful relay by wooldid; the blackout is the Stata nonsymmetry screen
  zeroing reghdfe's degenerate V (layer 3); what it hides is regime A (rank-1) plus a
  genuinely near-degenerate att horizon. Bucket input: ame — as regime A (correct SE
  possible from cross-cohort dispersion/bootstrap; the would-be posted V is the
  self-refuting rank-1); att — mostly bucket-(iii) (positive plug-in variances exist,
  lost in transit), except the near-zero horizon, which needs its own refusal.
- Fix sites: same as B (post-L1099 zero-V detection, warning surfacing, regress
  fallback), plus the margins-level screen means even a healthy e(V) does not guarantee
  posted SEs — wooldid's relay would also benefit from checking r(V) directly
  (L1771/L2186 blocks) when r(table) comes back all-missing.

## F. d3 (cohort-constant dose, all VCEs): fabricated-looking zeros with honest-but-thin
## signaling

- Surface (wooldid.log L8159-9047): 165 multicollinearity cautions (L1309-1323; one per
  treated-cell dose slope: dose constant within cohort ⇒ constant within cell ⇒ S_gh
  collinear with the cell dummy ⇒ dropped o.); then the L3318 branch: "not identified
  for 10 of 20 ... placeholder entries" — att rows keep SEs, all 10 ame rows display
  b = 0, se = missing.
- Margins-level confirmation (probe_fix4, frozen d3_r1.dta): the model has ALL 167
  dose-carrying coefficients omitted (o.) — the AME derivative is identically zero in
  the fitted model — and margins returns r(b) = EXACT zeros x8, r(V) = 0 exactly,
  se row all missing. The binary att margins is untouched (SEs .1025-.1682, matching the
  struct-pass displayed table digit-for-digit). vce(robust) reproduces the same 8/8
  missing ame SEs — the omission is X-structural (dose constant within cohort ⇒ constant
  within every cell), independent of the VCE, explaining "all VCEs" in the struct pass.
  Note the model V posts fine here (no screen trip: 166 positive dummy diagonals; the
  zeros sit on omitted columns, which Stata tolerates — POSTEXP S1).
- What e() carries: e(b) has hard 0s in the 10 ame slots (a pipeline consuming e(b) gets
  zeros indistinguishable from genuine estimates); e(V) diagonal has −1e100 sentinels for
  those slots (renders as missing SE in Stata surfaces; raw readers see −1e100);
  e(V_es_ame)/e(b_es_ame) absent (f2ok=0); e(wdidmainresults) rows carry 0/missing.
- Would a user reading only e()/the table know the ame is unidentified? The SE side —
  yes ("not identified" caution + missing SEs). The point-estimate side — NO: the 0s are
  displayed/posted without any flag that they are structural zeros of a model whose
  treatment-dose terms were all dropped; the L1313-1315 caution says only "could imply
  your results are not identified, depending on exactly what was dropped" and requires
  verbose to see what was dropped.
- Bucket input: (iv) for the ame estimand as specified — no within-cohort dose variation,
  nothing estimable; the correct response is refusal/steer, not an SE. Fix sites: the
  multicollinearity screen L1244-1332 already identifies every dropped treated-cell dose
  slope (the `catastrophe' counter); keying "all contreat cell slopes dropped ⇒ declare
  ame not estimable, blank/flag the ame rows" requires only that counter + the contreat
  flag; assembly/display sites L2820-2898, L3279-3327.

## Task 7 — F2 suppression consistency across B/C/E

One shared predicate for every family and regime: `f2ok` (L3352-3358). b and V share its
gate (L3359-3363), so blackout regimes return NEITHER e(b_es_ame) NOR e(V_es_ame) (nor
the att twins). Empirically confirmed in regimes B and E (all four ABSENT while e(b)
posts 1x20), and in C via the struct-pass surface (identical caution + VES absent;
the struct harness's "VES|1x1" line is its coercion artifact, latent issue 7).
Consumer-safety answer: consumers do NOT get "b without V" —
they get nothing from F2 in blackout regimes (safe, though the point estimates remain
only in e(b)/e(wdidmainresults)); in regime A they get b WITH a silently degenerate V
(unsafe — the actual gap).

## POSTEXP — which V property trips Stata's ereturn-post zeroing screen

Controlled experiment (probe_fix4: an eclass program `ereturn post`s hand-built 3x3 Vs,
then reads e(V) back):

| case | V property                                              | result                            |
|------|---------------------------------------------------------|-----------------------------------|
| S1   | symmetric PSD, diag (1, 1e-30, 0)                       | posted UNCHANGED                  |
| S2   | negative DIAGONAL (−1e-22), zero off-diags              | posted UNCHANGED (negative kept!) |
| S3   | asymmetry 1e-12 in one off-diagonal                     | posted unchanged                  |
| S4   | asymmetry 1e-6 in one off-diagonal                      | **ZEROED** (with the warning)     |
| S5   | symmetric, diag (1, 1e-260, 1)                          | posted unchanged                  |
| S6   | symmetric INDEFINITE, off-diag 2e-22 > diag 1e-22       | **ZEROED** (with the warning)     |
| S7   | symmetric indefinite at O(1): off-diag 2 > diag 1       | **ZEROED** (with the warning)     |
| S8   | symmetric PSD, exact-zero 2x2 block + healthy diag      | posted unchanged                  |

(S4/S6/S7's zeroing each printed exactly "warning: variance matrix is nonsymmetric or
highly singular." — the warning and the zeroing are one event.) **The trip class is
invalid covariance geometry: asymmetry beyond a tolerance in (1e-12, 1e-6], or
indefiniteness carried by off-diagonals (implied "correlations" > 1). PSD
rank-deficiency, exact-zero blocks, extreme scale, and even a pure-diagonal negative
entry all pass.** Mechanism of the blackouts, completed: reghdfe's V = D·M·D on these
designs has a (near-)null block whose diagonals are ~0/±1e-22 while its off-diagonal
roundoff is the same order — off-diagonals exceed what the diagonals allow — and/or
carries product asymmetry; Stata's screen fires, `ereturn post` stores e(V) ≡ 0, the
display shows every SE missing. (Inferred for the real cases — the pre-post V is
destroyed — but the trigger class is now bracketed by the experiment, and regress's
posted V on the same data shows the same block structure that S6 caricatures.) The same
screen exists inside margins: on the d4 regress path margins itself printed the warning
and blanked the ENTIRE r(table) se row while r(V) retained positive values (regime E).
Observed trip rates by design: A (basic, d1/d1c-typical): 0/16+; C (d1c basic): 1/16
realizations; E (d4 within): 3/4; B (byCohort): 4/4. The realization-dependence of C
and E is exactly the y-dependence of roundoff in the meat.

## Latent issues on these paths (task 8)

1. **The A/B silence-vs-caution asymmetry has a one-line cause: every detector tests
   only missingness.** L3307-3309 (e(V) caution) and L3352-3358 (F2 gate) fire on
   missing/non-positive diagonals; regime A's positive rank-1 V and regime D's
   positive understated V pass every gate silently, while B's zeroed-then-missing V
   trips them all. Same root cause (CRVE cell-saturation), opposite user experience,
   purely because one failure surfaces as "." and the other as a positive number. The
   "not identified" detector cannot fire on A as written; no detector inspects the
   structure (rank, dispersion, equality) of a present V. Fix leverage: the F2 site
   already loops over the captured r(V) diagonal (L3355-3358) — a rank/equality test
   there is a few lines.
2. **wooldid's `qui reghdfe` (L1099) swallows the only mechanistic diagnostics.** In
   every blackout the estimation itself printed "warning: variance matrix is
   nonsymmetric or highly singular" (Stata built-in), "warning: missing F statistic;
   dropped variables due to collinearity or too few clusters" (reghdfe5.mata L2706),
   and "note: dose is probably collinear with the fixed effects". None reach the user;
   the downstream caution names no cause. A post-L1099 check "e(V) exists and is
   identically zero ⇒ explain + point at VCE choice" costs one mreldif.
3. **The blackout wipe itself is upstream of wooldid** (reghdfe's numerically degenerate
   V → Stata `ereturn post` zeroes it) — wooldid is a faithful relay — but wooldid is
   the only layer that COULD notice (it sees e(V)≡0 with rc 0) and currently doesn't.
   Also actionable: wooldid already contains a plain-regress estimation path
   (`unconditionalse`/`customvce`, L1102-1108) that on the same data produces positive
   ATT SEs where the reghdfe path blacks out.
4. **AME point estimates under byCohort are normalization-dependent:** on identical data
   the reghdfe and regress parameterizations drop different collinear columns and the
   displayed ame differs (.723 vs .554 at horizon 1) while att is invariant. The
   byCohort ame estimand is thus partly an artifact of collinearity-drop order — worth
   a theorist ruling independent of the SE question.
5. **margins carries its own copy of the degeneracy screen and blanks SEs wholesale.**
   V_jj = 0 renders as se = "." in r(table); worse, when margins' delta-method V is
   near-degenerate it prints the same "nonsymmetric or highly singular" warning and
   blanks the ENTIRE se row even for quantities whose r(V) diagonal is healthy-positive
   (observed: d4 regress path, att horizons with variances .01-.03 blanked because one
   horizon sat at 7.9e-10). So a zero-variance quantity, a never-computed variance, and
   a healthy variance sharing a row with a degenerate one are all indistinguishable
   downstream; wooldid's missingness predicate conflates "margins couldn't" with
   "variance is exactly zero" with "collateral blanking".
6. **Assessment-harness nondeterminism (not wooldid.ado's results themselves):**
   amese_cases.doh builds data whose y-dose pairing depends on tie order in its bysort
   steps (duplicate rows after expand; rnormal if _n==1). What the evidence shows:
   same-seed builds reproduce bit-for-bit from a fresh session AND after plain
   reghdfe/margins/bysort activity (probe_cef's 17th-iteration d1-s0 matches
   fresh-session probe_bd exactly), but a build issued after a WOOLDID run diverged
   (p_d1c.log: identical att rows, different ame rows — the fingerprint of a
   within-cell permutation; plausibly wooldid's gtools calls (gegen/gstats/hashsort)
   perturb the tie-break state, though the precise agent inside wooldid was not
   isolated). Consequences: the prereg's "deterministic structure pass" is only
   att-deterministic across arbitrary session states; and since mc_amese.do calls
   wooldid every rep, MC per-rep seeds do not pin the ame data. Fact C should be
   re-labeled "knife-edge realization event (~1/16)".
7. **Struct-harness coercion bug:** `cap matrix V = e(V_es_ame)` (amese_structure.do)
   never fails when the e-matrix is absent — Stata coerces the undefined e() to a 1x1
   missing matrix with rc 0 — so the script's "VES|absent" branch is unreachable and
   absent matrices print as "VES|1x1|rank|0". Use `cap confirm matrix e(V_es_ame)`.
8. **e(V) = diag(se²) with fabricated zero off-diagonals** (L3280, acknowledged
   L3266-3268): cross-estimate lincom/test on e() is silently wrong in ALL regimes; F2
   is the sanctioned joint surface — which is exactly why its silent rank-1 posting in
   regime A matters most.
9. **The −1e100 sentinel** (L3313) renders correctly in Stata's display/inference
   surfaces but reaches raw e(V) consumers as a numeric "variance" with no
   machine-readable marker (regime F).

## Probe file index

- probe_bd.do/.log — regimes A/B/D wooldid e-surfaces on frozen d1_r1.dta; B manual
  reghdfe (V≡0 discovery); mkvars validation (aborted at a mata parse quirk; completed
  by probe_fix).
- probe_cef.do/.log — d1c 16-seed sweep (1 blackout, saved) + d1 16-seed sweep (0
  blackouts); aborted at the same mata quirk before E/F.
- probe_fix.do/.log — B margins detail + regress discriminator + seed-2 determinism;
  A same-file contrast; D faithfulness (all exact); CBLK dissection + regress
  discriminator (aborted at missing d1c_clean.dta).
- probe_fix2.do/.log — verbose(1) reghdfe VCE trace for B (proves reghdfe's mata did
  not zero V; q healthy) (aborted at variable collision).
- probe_fix3.do/.log — verbose(1) trace on d1c_black (screen attribution); d1c_clean
  rebuild + contrast (Var(phihat) fingerprint); dose forensics; regime E wooldid
  e-surface + model-level e(V)≡0 (aborted at a nosample/margins interaction).
- probe_fix4.do/.log — regime E margins detail + regress discriminator + d5 contrast +
  extra d4 seeds; regime F margins detail (cluster + robust); POSTEXP S1-S5.
- probe_postexp2.do/.log — POSTEXP S6-S8 (symmetric-indefinite and zero-block cases).
- mkvars.doh — reconstruction of wooldid's internal _i/_t/targetflag/null/ringer + ES
  window (mirrors wooldid.ado L543-658, L666, L1228, and the validate bump L3998-3999;
  validated by exact digit-match of the manual margins b to wooldid's displayed tables:
  B ame row .7234028/-.1309208/... == struct-pass d1_bycoh_ci ES_ame table, and CBLK att
  row -.1399923/... == struct-pass d1c_basic_ci ES_att table).
- d1_r1.dta, d1c_black.dta, d1c_clean.dta, d3_r1.dta, d4_r1.dta, d5_r1.dta — frozen
  realizations every cross-config claim was re-run against (removes the sort-RNG
  confound).
