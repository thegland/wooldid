# LENS2 on FCW-03 — realism and impact (third, independent lens)

## The question I am the tiebreaker on
Prober: the 149/150 over-coverage at cluster()==i is an ADDED cohort-composition term
injected by the nesting re-id leg, and the record's stated cause has the wrong sign.
Skeptic: PARTIAL/low — the record already retracted that cause the same day, already
states the composition mechanism as a hypothesis, and the re-id attribution is wrong.

My lens: how likely is a real applied economist to hit this, what do they lose, would
they notice.

## Established from the record BEFORE running anything (all verified by me)
1. `test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md` F3/F4 amendment 3 DOES carry the
   retraction ("A wrong analogy, retracted here so it does not propagate") and DOES state
   the mechanism as a hypothesis ("conservatism of the percentile cluster bootstrap where
   every cohort is its own resampling cluster"). The skeptic quoted it accurately.
2. BUT the retracted analogy is STILL LIVE, verbatim, in the header of
   `test/suite/mc/mc_fineclus.do` lines ~36-40: "it is the same effect the help already
   describes at 'What level to cluster at' - at or above the i level the covariance stops
   seeing the treated cells' own variation." That is the file a reader of the MC opens
   first, and it is the text the prober quoted. The amendment's own stated purpose
   ("so it does not propagate") was not carried out in that file.
3. The MC's DGP (`test/suite/mc/mc_fineclus_dgp.doh`, dgpE) sets
   `tau = mult*0.2*nu if mod(state,2)==1` OUTSIDE the seeded block. tau_s is a
   deterministic function of the cohort's size and is IDENTICAL in all 150 draws; only
   the unit RE `mu` and the AR(1) `eps` are redrawn. CONFIRMED by reading the file.
   So the MC measures the FIXED-COHORT sampling distribution.
4. The re-id attribution is refuted, but NOT by the argument the skeptic led with.
   r03b compares se(att) on the re-id leg vs the count path and finds them bit-identical
   (.0709192322 both) — but se(att) is exactly the quantity the prober's own control
   showed is regime-invariant, so that comparison cannot test the att_i claim.
   The refutation is in r03b's OTHER pair: on the COUNT path (no re-id possible at all),
   se(att) is .0709 under HET and .0026890 under HOM — a 26x composition gap on an
   estimand that never touches the re-id leg. The composition term comes from resampling
   COHORTS AS CLUSTERS, not from the re-id leg. The leg only decides att_i's SHARE of it
   (full multiplicity vs presence-only) — that is F-B/F-F, a different finding.

## What I am measuring that neither agent measured
The over-coverage is an ESTIMAND MISMATCH between the MC's truth and the bootstrap's
target, and the realism question is which target a user has and what the two numbers on
their screen look like.

- L2A.do  analytic arms, cluster(state)==i, 150+150+120 draws
    FIX    published dgpE, tau_s held fixed  -> SD_FIX  (the SD the 149/150 was judged against)
    SUP    same design, tau_s redrawn iid {.4,1.2} each draw -> SD_SUP (cohorts resampled)
    LOWFIX tau .6/1.0 (Var_s=.04) -> realistic-heterogeneity calibration
  Predicted: SD_SUP^2 = SD_FIX^2 + .16/20 exactly.
- L2B.do  25 HET + 15 LOW bootstrap draws, reps(99), each PAIRED with the analytic run on
  the SAME dataset, so boot/analytic — the gap a user sees between two runs of their own
  data — is on one row.
- STEP 0 of L2A asserts my dgpL FIX arm reproduces the published dgpE bit for bit.

## Status
- [x] record read and verified
- [ ] L2A run
- [ ] L2B run
