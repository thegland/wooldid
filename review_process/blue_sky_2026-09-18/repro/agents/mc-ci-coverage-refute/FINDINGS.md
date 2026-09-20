# SKEPTIC report — refutation attempt on the mc-ci-coverage findings

Build asserted in every probe log (`which wooldid` -> `.../scratchpad/bluesky/build/wooldid.ado`).
Every dataset is generated inside the probe do-file. No repo fixture, no prober fixture, no
`include` of anything the prober wrote. I read their do-files; I ran none of them.

My DGP throughout (deliberately unlike theirs): 60 units, `i = cg` (3 onset cohorts 2004/2006/2008
plus a never-treated group), 12 annual periods, a unit random intercept sd 0.70 that `i = cg`
does NOT absorb, AR(1)(0.4) idiosyncratic errors sd 0.55, a common linear time trend, and
WITHIN-cohort effect heterogeneity (a unit multiplier `1 + N(0,0.25)` on cohort effects
1.00 / 1.15 / 0.85). The band probe uses a relative-period-only effect so every weighting scheme
gives the same horizon truth.

---

## VERDICTS

### mcci-01 impvalues coverage — CONFIRMED (critical), with one mechanism correction
own repro: `repro_S1_impvalues.do`  (also `p1_quick.do`, `p2_mcimp.do`, `p4_mech.do`)

140-draw Monte Carlo, three arms on the SAME draws (`p2_mcimp.txt`, `an_p2.py`):

| arm | mean(b) | sd(b) | mean(se) | se/sd | cov vs att | cov vs cluster-equal truth | cov around own MC mean | df_r |
|---|---|---|---|---|---|---|---|---|
| impvalues cluster(unit) — cluster FINER than i | 0.9977 | 0.1376 | 0.0233 | 0.169 | 0.257 | 0.264 | 0.243 | 44 |
| impvalues cluster(cg) — cluster == i, textbook IM | 0.9977 | 0.1376 | 0.1092 | 0.794 | 0.971 | 0.964 | 0.964 | 2 |
| analytic cluster(unit) | 1.0145 | 0.1373 | 0.1252 | 0.912 | 0.957 | 0.950 | 0.936 | 59 |
| ES horizon post0 (impvalues) | 0.9956 | 0.1278 | 0.0286 | 0.224 | 0.314 | 0.329 | 0.336 | 44 |
| ES horizon post1 (impvalues) | 1.0023 | 0.1359 | 0.0285 | 0.209 | 0.343 | 0.336 | 0.286 | 44 |

3 binomial SE at p=.95, n=140 is 0.055, so 0.257 sits ~35 SE below nominal. The coverage number is
not an artefact of picking the wrong target: the interval misses its OWN Monte Carlo mean 76% of
the time, and mean(se)/sd(b) = 0.169 settles it without reference to any target at all.

THE MIDDLE ROW IS THE DISCRIMINATOR THE PROBER NEVER RAN. Swap only the cluster variable, keep
`i`, the data and the draws: with `cluster()` at the level of `i` the IM machinery is fine
(se/sd 0.794, coverage 0.97). So this is not "IM is broken"; it is specifically the finer-than-`i`
geometry, which is the geometry the help's "What level to cluster at" (sthlp:255) calls
"the intended practice".

EXACT MECHANISM, proved arithmetically (PART 1 of the repro). On one dataset,
`impvalues cluster(cg)` reports se 0.105582803 on 2 df, so the sd over the 3 distinct
per-cohort reads is 0.105582803*sqrt(3) = 0.18287478. If the 45 per-unit reads are those same
3 values each repeated 15 times, the reported se must be
`0.18287478*sqrt(2*15/44)/sqrt(45) = 0.022510329274`. Observed: `0.022510329`. RATIO
1.000000000000000. The 45 reads take exactly 3 distinct values.

CORRECTION to the stated mechanism. The prober wrote that the spread is "cross-cohort effect
heterogeneity, with the sampling component absent entirely". That is wrong, and the control arm
proves it: those 3 distinct reads DO carry the sampling component — the same 3 numbers, referred
to sqrt(3) and t(2), give se/sd 0.794 and 0.97 coverage. What is actually wrong is MULTIPLICITY:
each distinct read is counted once per cluster, so `ttest` divides by sqrt(G_clusters) and refers to
t(G_clusters - 1) where the honest divisor is sqrt(#distinct reads) and the honest reference is
t(#distinct reads - 1). The deflation factor is therefore ~sqrt(G_clusters/#distinct reads)
(4.69x here) and the df is inflated from 2 to 44. "Sampling component absent entirely" is the
limiting case only — reached when every cluster has an identical weight profile, where the reads
collapse to one value, `se` hits the degeneracy scrub and everything is blanked (PART 3).

Everything else in the claim reproduces: the ES horizons carry the same construction (I measured
them; the prober asserted them), there is NO caution of any kind in the log, and the geometry is
the one the help recommends. Severity critical stands: a silently wrong number, on a realistic
panel, on the inference surface the help sells as the fallback.

SCOPE I ACTUALLY TESTED: the pooled `att` row and the `att`-family event-study horizons, in binary
mode. The claim's surface also names the continuous-treatment `ame` rows (ado:2320/2375/2715); those
share the same `margins ... over(cluster)` construction by code inspection but I did not run them.

Not T4 (known_bydesign section 1): every run here is `att` alone, and the headline geometry is a
cluster nested INSIDE `i` rather than a cross-cutting cluster with a non-`att` estimand. The
reported point estimate is the cluster-equally-weighted mean (that part is T4 and I do not re-file
it); what is new is that its STANDARD ERROR is deflated by the multiplicity of the reads.

### mcci-02 simultaneous band — CONFIRMED (medium), and the open mechanism question is now answered
own repro: `repro_S2_band.do`  (MC in `p3_mcband.do` / `p3_mcband.txt`, analysis `an_p3*.py`)

250 draws, G = 12 clusters, `i = cg`, `cluster(unit)`, `esprelength(3) espostlength(2)`,
`simultaneousbands level(90)`, 5 banded horizons, `e(esband_engine) = analytic-simulated`,
`e(df_r) = 11`:

- measured SIMULTANEOUS coverage 0.788 at nominal 0.90 (3 binomial SE = 0.057, so ~6 SE low);
- POINTWISE coverage on the same runs 0.872-0.928 per horizon with se/sd 0.96-1.16, so the
  covariance and the pointwise machinery are sound and the deficit is the band's alone;
- mean reported `e(esband_crit_att)` = 2.1536; the ORACLE critical value (the empirical 0.90
  quantile of `max_h |b_h - truth_h| / se_h`) is 2.7574. The band is 78% of the width it needs.

THE DECOMPOSITION THE PROBER LEFT OPEN (400,000-draw simulations on the MC's own empirical
horizon-error correlation; coverage re-computed by rescaling each draw's own reported crit):

| reference distribution for the sup | crit | vs reported | simultaneous coverage |
|---|---|---|---|
| as wooldid reports it (normal sup on the ESTIMATED correlation) | 2.1536 | — | 0.788 |
| normal sup on the TRUE correlation | 2.2112 | +2.7% | 0.804 |
| t(11) sup, one common scale factor | 2.4950 | +15.8% | 0.856 |
| t(11) sup, per-horizon scale factors | 2.6054 | +21.0% | 0.888 |
| MC oracle | 2.7574 | +28.0% | 0.904 |

So the t-versus-z reference is the MAJORITY component of the 0.112 deficit: it accounts for 0.068
(common scale factor) to 0.100 (independent per-horizon scale factors), and the truth is between the
two, because cluster-robust standard errors across horizons share the same cluster draws and so
have correlated but not identical scale errors. Plugging the estimated correlation in as if known
contributes ~0.016 (it biases the crit DOWN 2.7%, because sample correlations from 12 clusters are
too large in magnitude, which shrinks a sup quantile). At most ~0.03 is left unattributed. The
prober's open question is answered in the direction they guessed; their "probable" tag was the right
call and the mechanism can now be called established as the dominant term rather than the whole. The claim's own numbers differ
from mine (0.752 vs 0.788 at nominal 0.90) — DGP-specific, same direction, same size class.

The self-contained repro `repro_S2_band.do` reproduces this on its own 200 draws:
simultaneous coverage 0.780 at nominal 0.90 (3 binomial SE 0.064), mean reported crit 2.1537,
oracle crit 2.8292, band = 0.761 of the width it needs, pointwise coverage 0.875-0.950.

Severity medium stands: opt-in, default off, and the deficit shrinks with G (the t/z gap on the sup quantile is
~13% of width at G=12; the corresponding univariate t(39)/z ratio at G=40 is ~3%, computed rather
than measured). Code: `wooldid_esbanda`, ado:16711-16746 —
`Z = L * rnormal(m, nd, 0, 1)` with no df adjustment, while the pointwise intervals printed on the
same table use `invttail(G-1)`. Help (sthlp:524): "a simultaneous band covers every horizon of
the family at once with that probability". It does not.

### mcci-03 esprelength moves the att — PARTIAL (low): real, but the mover is `esprelength()` alone
own repro: `repro_S3_esprelength.do`

One unchanged dataset, five runs, `e(N) = 720` and `e(df_r) = 59` in all of them:

| spec | att | se | e(escomp_refN) |
|---|---|---|---|
| (no ES) | 0.929543611 | 0.124672021 | (not returned) |
| `espostlength(2)` | 0.929543611 | 0.124672021 | 270 |
| `esprelength(2) espostlength(-999)` | 0.858508187 | 0.132884576 | 180 |
| `esprelength(2) espostlength(2)` | 0.858508187 | 0.132884576 | 180 |
| `esprelength(3) espostlength(3)` | 0.894165051 | 0.151295323 | 135 |

The behaviour is real and the stated cause (the reference-normalization pool is
`relyr < -esprelength` plus the forced first observation, so requested pre-horizons leave the pool)
is confirmed by `e(escomp_refN)` falling 270 -> 180 -> 135. There is NO caution anywhere in the
log, and the `esprelength()` option entry (sthlp:729) does not mention the effect — only the
`e(escomp_refN)` Stored Results entry (sthlp:1368-1370) does. So the doc-gap half of the claim holds.

Why PARTIAL: the claim is titled and stated as "adding `esprelength()`/`espostlength()`", which
implicates `espostlength()`. It does not move anything — `espostlength(2)` reproduces the no-ES att
and se to the last digit. And "inflates its confidence interval by up to 60%" is their dataset's
number, not a property of the option: the same option at the same setting gives +21% on mine. The
finding should read "`esprelength()` > 0 silently re-bases the headline att and widens its interval,
by an amount that depends on how much of the pre-history the request consumes".

Severity low stands.

---

## NEW FINDINGS

### skep-01 impvalues also under-covers at `i = unit` + `cluster(unit)`, where IM is legitimate
own repro: `repro_S4_imunit.do` (MC `p5_imunit.do` / `p5_imunit.txt`)

100 draws, `i = unit`, `cluster(unit)` — 45 treated clusters, each owning its own
cohort-by-relative-period coefficients, i.e. the geometry in which Ibragimov-Muller is textbook
legitimate and in which no multiplicity of reads exists (df_r = 44 = 45 - 1):

| arm | sd(b) | mean(se) | se/sd | cov vs att | cov vs cluster-equal truth | cov around own MC mean |
|---|---|---|---|---|---|---|
| impvalues (100 draws, p5) | 0.1249 | 0.0827 | 0.662 | 0.820 | 0.830 | 0.840 |
| analytic  (100 draws, p5) | 0.1273 | 0.1153 | 0.906 | 0.890 | 0.890 | 0.910 |
| impvalues (80 draws, repro_S4) | 0.1172 | 0.0823 | 0.702 | 0.850 | 0.863 | 0.863 |
| analytic  (80 draws, repro_S4) | 0.1187 | 0.1147 | 0.967 | 0.900 | 0.900 | 0.900 |

3 binomial SE is 0.065 at n=100 and 0.073 at n=80, so both coverage figures sit 4-6 SE below
nominal and both se/sd ratios are ~0.66-0.70. The finding stands on the coverage and se/sd numbers alone. The mechanism is HYPOTHESIZED, not
measured: here the 45 reads really are 45 distinct numbers (no multiplicity, unlike mcci-01), and
the natural explanation is that they are positively correlated — every treated unit's effect is
measured against the same never-treated comparison group and the same estimated period effects —
while IM's `sd/sqrt(G)` assumes independence across clusters. I did not measure that correlation.
The analytic arm's 0.890 is the documented at-or-above-`i` understatement (known_bydesign section 1);
what is new is that `impvalues` is WORSE than the path it is offered as an alternative to, in the
commonest geometry an applied user will write down. Marked `probable` because it rests on one DGP.

### skep-02 the IM degeneracy message misdiagnoses noisy data as deterministic  [SEEDED]
(The prober's mcci-01 OBSERVED field already reports this from their probeD.log; I reproduced it on
my own data and split it out as a separate UX item rather than leaving it inside a critical finding
about a silently wrong number. Credit for spotting it is theirs.)
own repro: `repro_S1_impvalues.do` PART 3 (also `p4_mech.do` block C)

`i = cg`, `cluster(reg)` with 5 regions that cross the cohorts and hold equal numbers of each, on a
panel whose idiosyncratic sd is 0.55 and unit sd 0.70. Every `se/t/p/lb/ub` is blanked and the
message reads: "the per-cluster estimates were numerically identical ... A collapsed IM standard
error usually reflects deterministic or noiseless data - please verify your input if this was not
expected." The data is not deterministic and the user has nothing to verify: the collapse is a
property of the CLUSTER GEOMETRY (identical weight profiles over the `i` x relative-period cells),
not of the data's noise. This is an observable failure rather than a silent wrong number, so it is
low — but the message sends the user to look in the wrong place, and it is the same root cause as
mcci-01 at its limit. A message naming the geometry would be right.

---

## NEGATIVE RESULTS (things I attacked and could not break)

- POINTWISE analytic intervals at `cluster()` below `i`. 250 draws x 5 horizons plus 140 draws of
  the pooled `att`: se/sd 0.91-1.16, coverage 0.872-0.957 at the nominal level. The pointwise
  machinery on the analytic path is calibrated; mcci-02 is a band-only defect and mcci-01 is an
  impvalues-only defect.
- impvalues at `cluster() == i` (3 treated cohorts, df 2): se/sd 0.794, coverage 0.964-0.971.
  Slightly conservative, as t(2) should be. The IM code path itself is correct.
- The ESTIMAND under impvalues is the cluster-equally-weighted mean (T4, already ruled); I
  measured coverage against BOTH targets and against the estimator's own MC mean, and all three
  agree, so none of the verdicts above turns on the estimand question.
- `espostlength()` does not perturb the att, the se, `e(N)` or `e(df_r)` — bit-identical to the
  no-ES run.
- `e(escomp_refN)` behaves as documented: absent with no ES window, present and shrinking once
  one is requested.

## FILES
- `p1_quick.do/.log` — first contact: mcci-03 four runs, mcci-01 five arms on one dataset.
- `p2_mcimp.do/.log/.txt`, `an_p2.py` — 140-draw 3-arm impvalues MC.
- `p3_mcband.do/.log/.txt`, `an_p3.py`, `an_p3b.py`, `an_p3c.py` — 250-draw band MC + decomposition.
- `p4_mech.do/.log` — which ES option moves the att; the i=unit and crossed-cluster boundaries.
- `p5_imunit.do/.log/.txt` — 100-draw i=unit impvalues MC.
- `repro_S1_impvalues.do`, `repro_S2_band.do`, `repro_S3_esprelength.do`, `repro_S4_imunit.do`
  — the four self-contained repros, each runnable with nothing but the pinned build.
