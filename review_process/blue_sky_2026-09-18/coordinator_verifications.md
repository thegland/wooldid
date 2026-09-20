# Coordinator-verified notes (things I checked myself, not agent claims)

Written 2026-09-17 during the probing wave.

## semielasticity / eydx / eyex: NOT a regression, and structurally invisible to the project's own oracles

Checked directly on both builds:
- pinned build ado:1476 `qui reghdfe ... , absorb(...) `vcespecifier' tol(1e-`regtol')`  -- **no `d`, no savefe, no residuals()**
- published v1.3 (main) line 1147: the same reghdfe call, **also with no `d`**
- all three ppmlhdfe call sites (ado:1462, 1468, and the rc-3200 fallback) DO pass `d`, with the
  in-code comment that ppmlhdfe requires it "for predict, which margins requires"

So the absorbed-FE omission on the reghdfe path predates the 2026 campaign and ships in the
released v1.3.

**Why this matters for how the finding is reported.** Three of the project's standard
verification instruments are structurally incapable of catching it:
 1. old-vs-new dual-build comparison -- both builds are wrong in the same way, so the diff is clean
 2. the bit-exact baseline battery -- it pins current behaviour, and current behaviour is the defect
 3. external oracles -- jwdid/csdid/did_imputation were mapped on LEVEL estimands (dydx), which the
    finding says are unaffected; the elasticity rows were never externally benchmarked
This is precisely the defect class the blue-sky review was commissioned to find: a wrong number that
every existing green instrument agrees with.

Status: the underlying finding is single-agent and NOT yet refuted. What is verified here is only
the code-archaeology half (the missing `d`, present in v1.3, contrasted with ppmlhdfe). The
measured 65% divergence and the sign-flipping pre-treatment coefficients remain to be independently
reproduced.

## Verified statically, earlier this session
- README:141 / sthlp:164 normalized-zero claim: FALSE (measured, bluesky/coord/normzero/nz.do)
- sthlp:1469 `ereturn post e(b_es_*)` handoff destroys e(): the command given is `post`, not `repost`
- wooldid.pkg pins reghdfe 5.7.3 (sthlp:46 repeats); project baselines on 6.12.3; no version check
- wooldid.pkg Distribution-Date: 20260819, stale against a build changed through 2026-09-16
- help does NOT restrict impvalues cluster geometry; sthlp:375 claims the crossed case was measured
  fine at 1e-12 (read accuracy, not IM validity); sthlp:1097 recommends impvalues as the valid fallback

## CORRECTION to my own earlier framing (impvalues), 2026-09-17

I wrote that sthlp:1097 "recommends impvalues as the valid fallback where analytic inference fails"
and used it as part of the realism case. That is WRONG and must not reach the author in that form.
sthlp:1097 sits INSIDE the contreat cell-containment paragraph, whose regime is cluster CONTAINING
the cell ("automatic when clustering on the cohort identifier i itself") - i.e. the geometry where
impvalues was independently MEASURED to work (coverage 0.935 at cluster coarser than i, 0.810 at
cluster == i). That recommendation is made in the safe regime and is not itself wrong.

The realism case rests instead on:
 - sthlp:67 and sthlp:255 - the package's headline guidance: i = onset cohorts, and "The intended
   practice is to cluster BELOW the level of i". That is the failing geometry verbatim.
 - sthlp:257 - the cluster() entry imposes exactly ONE rule on impvalues: a single cluster variable.
 - sthlp:375 - normalises the crossed case by reporting the reads "exact - measured at 1e-12 both
   for a cluster nesting i and for one crossed with it". Read accuracy != IM validity; only the
   first was ever measured.

## NEW sub-finding from the impvalues refuter: the analytic arm COARSER than i is far worse than documented

Measured, R=200: the analytic clustered arm at cluster STRICTLY COARSER than i has se/sd 0.179 and
95% coverage 0.335. sthlp:255 quantifies the at-or-above-i understatement as "roughly 0.66 to 0.75
of the reference on the diagnostic fixtures". 0.179 is far outside that band. This CONTRADICTS a
documented measured range rather than rediscovering it, so it is first-class per the coordinator
note in known_bydesign.md. Keep it separate from the impvalues finding so it does not dilute it.

## THE HEADLINE CHAIN (coordinator-verified 2026-09-18) - the remedy works on SEs and breaks points

Two findings that looked contradictory are complementary, and together they are the most
important result of this review.

1. ANALYTIC SEs UNDERCOVER at or above `i`. Known; the project's own 500-draw MC (2026-09-14)
   measured 95% coverage 0.794-0.876 at i=unit + cluster(i).
2. THE DOCUMENTED REMEDY (sthlp:255) is "a finer cluster()" - "The intended practice is to cluster
   BELOW the level of i".
3. THE REMEDY GENUINELY WORKS ON THE SE, MEASURED. The mc-fineclus-analytic cell ran ~1,500 draws
   over seven DGPs: remedy arm (i=cohort, cluster(unit)) median SE/SD 0.970, coverage 0.935;
   at-i control 0.589/0.739; cross-cutting cluster 0.999/0.950. Point estimates are EXACTLY
   invariant to cluster() - mreldif 0 on every one of 283+200+200+91 draws. The seeded positive
   control reproduced the project's own numbers (0.852/0.796 vs their 0.838/0.794), validating
   the harness.
4. BUT REACHING THAT GEOMETRY BIASES THE POINT ESTIMATE ON RAGGED PANELS. On a one-row-per-
   state-year panel the only way to get cluster() below `i` is to coarsen `i` to the onset WAVE,
   which replaces the state fixed effects with a single wave intercept. The sub-`i` permanent
   heterogeneity is then no longer absorbed.
5. `fe(state)` OR `controls(i.state)` REPAIRS IT EXACTLY - and sthlp:255's remedy sentence
   mentions neither.

COORDINATOR VERIFICATION (coord/remedy/rem.do, zero-noise DGP, true att == 1 exactly):
  balanced: i=state 1.0000000000 | i=wave 1.0000000000 | i=wave+fe(state) 1.0000000000
  ragged  : i=state 1.0000000000 | i=wave 0.9838810126 | i=wave+fe(state) 1.0000000000

MAGNITUDE - REPORT THE RANGE, NOT THE HEADLINE NUMBER. My ordinary 8% MCAR panel gives a 1.6%
bias. The refuter's 21% came from attrition landing on a TREATED state, and its -57% from a
treated state entering late; its dispersion ladder (sd(s_g) 0.5/1/2/4 -> 1.5/3.1/6.2/12.3%)
brackets my 1.6%. So: small under ordinary MCAR attrition, large when the raggedness hits a
treated unit. Do NOT lead with 21%.

This is the project's own recorded lesson applied (check the ORDER OF MAGNITUDE against the
mechanism before explaining a numeric difference).
