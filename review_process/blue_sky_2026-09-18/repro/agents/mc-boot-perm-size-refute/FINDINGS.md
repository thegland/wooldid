# SKEPTIC REPORT — refutation of the "MC size of permute()/bootstrap() p-values at few clusters" wave

All numbers below are from MY OWN generator and MY OWN draws. I read the prober's
do-files but ran none of them as evidence. Every repro named here is self-contained.

## VERDICTS

### MCBP-ESFBP — **REFUTED** (severity high -> low)
The claim: at G=10 with cluster()==i the default `att` analytic SE is 0.71-0.77 of truth
(rej .18-.24), and `esfixedbaseperiod` yields se/sd 0.93-1.01 and rej .06-.08, so the
base-period option "is the lever that fixes it".

My repro: `r01_analytic.do` (300 draws; theirs was 25 and ~49) and `r05_esfbp_altdgp.do`.

1. **The claimed fix does not hold at 300 draws.** G=10, cluster(st)==i, sharp null:
   default se/sd **0.714**, rej **.170**; + esfixedbaseperiod se/sd **0.809**, rej **.140**.
   Not 0.93-1.01, not .061-.080. A partial improvement, not a repair.
2. **On a multi-unit panel the "lever" runs BACKWARDS.** Same i=state, same cluster(st)==i,
   5 sub-units per state-year (N=500, 120 draws): default se/sd **0.632**, rej **.258**;
   + esfixedbaseperiod se/sd **0.557**, rej **.275**. A user who followed the finding's
   advice on a microdata state panel would make their standard error worse.
3. **The sd ordering in their evidence is backwards and the help predicts it.**
   sthlp:836: "pooled reference periods offer the possibility of greater precision ...
   they pool together more data". I measure sd(esfbp) > sd(default) on BOTH shapes
   (0.6356 vs 0.4742; 0.6019 vs 0.4493). Their sd(S2) < sd(S0) at n=25 was noise.
4. **The default's understatement is documented and ruled.** 0.714 is inside sthlp:251's
   "roughly 0.66 to 0.75", which also says "at that level the analytic covariance is the
   wrong object rather than a mis-scaled one" — the finding's "one of the two analytic SEs
   is simply wrong, and it is the default one" IS the help's sentence. known_bydesign §1
   records it as ruled three times.
5. **"an event-study display option, not an inference option" misreads the help.**
   sthlp:832: esfixedbaseperiod "will switch to employing a fixed reference period, meaning
   that ALL treatment effect estimates will be formed by comparing ..." — a specification
   switch, described as such, with its precision cost stated in the next paragraph.
6. Their S3==S2 sub-observation (the window is absorbed by esfixedbaseperiod) reproduces
   exactly: S2, S3 and S5 are identical to 4 dp on all 300 draws.

### MCBP-JOINTSIZE — **REFUTED** (severity medium -> none)
Measured on my own 300 draws: analytic `e(joint_horiz)` rejects a true sharp null at
**.6967** (ECDF .01/.05/.10/.20/.30/.50 = .567/.697/.750/.827/.887/.947); adversarial .2367.
- The project's record says **"65.7% at G=10 (NF-1b)"** (known_bydesign §1) and the ado's own
  SI-A1 comment (`wooldid.ado:1396`) says **"measured rejection rates 11-66% at nominal 5%"**.
  The finding quotes "11-36%" and calls .755 "2-7x that band"; its own Wilson CI contains .657.
- Worse for the claim, wooldid **prints a caution on exactly this design** (PS-64 ruling,
  ado:1435-1443) — verified firing verbatim in `r02_cautions.log:121`: "...the joint test's
  covariance is degenerate and its p-value should not be read as a calibrated test. wooldid's
  randomization inference (permute()) or bootstrap-based clustered SEs give a joint test that
  is valid under this design...". The prober ran `qui` and never saw it. Their EXPECTED
  ("a user ... is reading a statistic that fires three times in four") is that caution's text.

### MCBP-GUARD — **PARTIAL** (severity low -> low; mechanism and remedy refuted)
True: at exactly G=10 total with 6 treated, no few-cluster caution prints.
False: "the few-cluster caution keys on the TREATED cluster count" and "should key on the
total cluster count, or at least both counts". Both cautions exist. Census (`r02_cautions.log`,
bootstrap(standard), 6 treated throughout): G=8 and G=9 fire "Caution: only N clusters are
available for the bootstrap"; G=10/11/12 silent; G=10 with 4 treated fires "Caution: only 4
treated cluster(s) are available for resampling". The claimed silent band "(6-9 treated
clusters, any number of controls)" does not exist — 6-9 treated with total<10 is warned.
What survives is only that the total-cluster threshold is `< 10` and the prober's design sits
on the boundary. Also: the `< 6` treated rule is an explicit author ruling (F050, 2026-09-06)
whose own comment says the total caution "cannot see ... few TREATED clusters".

### MCBP-BAYES — see repro r03_boot.do (in progress at time of writing; final numbers in NOTES.md)
Novelty already refuted at source: ado:6491-6499 — the lines the finding cites — say
"the Bayesian bootstrap shares its [the pairs cluster bootstrap's] first two moments" and cite
MacKinnon-Webb badly sized. bootstrap(bayes)+att_i is REFUSED r(198) under default safety
with a message naming the understatement.

### MCBP-S7COVER — see repro r04_s7geom.do (in progress)
The §7 record (`git show c8c4e2c`) is M=150 at reps(199), i=state with MULTIPLE UNITS per
state (the companion arm is cluster(unit), "strictly finer than i=state"), and a DGP with
LIVE COHORT HETEROGENEITY in effects — not a sharp null. The finding says its geometry
"differs only in having two states per onset cohort instead of one". It differs on at least
four axes.
