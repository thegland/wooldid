# boot-reproduce -- FINDINGS (agent BR: bootstrap reproducibility, differential testing of the two legs)

Build: pinned `/private/tmp/.../bluesky/build/wooldid.ado` (sha c1d8979), loaded with `adopath ++`.
All probes generate their own data with the recon generator `gen_core.do`; no repo test suite used.

---

## BR-1  permsave() has no comma check; it silently writes a file whose NAME contains ", replace"
SEEDED (recon lead). CONFIRMED, certain.

`bootsave()` got a comma check on 2026-09-13 (RF-6, ado:6282-6289) precisely because
`bootsave(name, replace)` -- the obvious analogy to Stata's own `saving(f, replace)` -- was
"swallowed whole and wrote a file literally called `name, replace.dta` at rc 0".
`permsave()` (ado:5866-5881 parse check, ado:8837 writer) never got that check.

Observed (probe2.log):
  `... permute(shuffle) permreps(20) permsave(junkp, replace)` -> rc 0, and
  `file junkp, replace.dta saved`; the file is on disk as `junkp, replace.dta`.
  The same shape on bootsave: `bootsave(junkb, replace)` -> rc 198 with the RF-6 message.
Expected: permsave() should refuse the comma with the same message; the help (sthlp:495, 700)
documents permsave's name rules as bootsave's.
Repro: probe2.do (block P6).  Severity low-medium, realistic (a user copying Stata's `saving()` idiom).

---

## BR-2  The reps() resolution Caution is unreachable without cluster(), contradicting the help
SEEDED (recon lead). CONFIRMED, certain.

sthlp:534 states, unconditionally: "A rank-based p-value cannot be finer than its own grid, and
`wooldid` prints a Caution naming that grid whenever it is coarser than .05".
In the code the whole caution block (ado:6455-6478) sits INSIDE `if "`clus'" != "" {`
(the cluster-count sanity block), so no clustered run == no caution.

Observed (probe2.log):
  clustered   `att cluster(unit) bootstrap(standard) reps(12)` -> "Caution: with reps(12) the
              coefficient p-values fall on a grid of 0.167 ..." (printed).
  unclustered `att bootstrap(standard) reps(12)`               -> NO such Caution.
  unclustered `att bootstrap(standard) boottstat reps(4)`      -> NO such Caution
              (grid/floor 2/5 = .4; no p-value can reject at 5%).
  unclustered `att bootstrap(standard) bootdisplay(wald) reps(4)` -> NO such Caution.
Mitigation that limits severity: the POST-run footer still names the realized resolution.
What is lost is the pre-run warning and its actionable "use at least 39/40 replications" advice --
which is exactly what RF-2 (2026-09-13) moved this block for. permute()'s sibling caution
(permctl) is NOT conditional and fires on every run.
Repro: probe2.do (block P5).  Severity medium (doc-vs-code), realistic.

---

## NEGATIVE (clean): the two legs agree bit-for-bit on the clustered draws
SEEDED lead, attacked directly; CLEAN.

`bootstrap(standard) att cluster(c)` takes the fweight COUNT leg; adding `att_i` switches the
whole run to the physical RE-IDENTIFICATION leg (ado:378).  The code comment at ado:9340 claims
"same-seed draws equal the count path's".  Measured on bootsave() files at a common rseed:

  cluster(unit), i=unit (nesting)   : w_* cluster multiplicities BIT-IDENTICAL (max abs diff 0
                                      over 60 cluster columns x 30 reps); main_att replicate
                                      draws agree to 1.50e-15; reported att SE differs by 1 ulp.
  cluster(unit) + boottstat (reid)  : w_* identical; main_att draws 1.60e-15.
  cluster(region) (coarser, nesting): w_* identical; main_att draws 1.30e-15.
Repro: probe1.do + cmp1.do.  The claim holds; adding a second estimand does NOT move att under cluster().

## NEGATIVE (clean): Stata's own bootstrap prefix matches the count leg replication for replication
sthlp:532 claims the count leg "matches Stata's own {helpb bootstrap} prefix replication for
replication under a common seed".  Measured with the help's own recipe:
  `bootstrap attb = _b[main_att], reps(40) seed(2468) cluster(unit) idcluster(newc): wooldid ... fastnose`
  vs `wooldid ..., att cluster(unit) bootstrap(standard) reps(40) rseed(2468) bootsave()`
  -> max relative difference 2.56e-08 across all 40 replications, which is exactly the FLOAT
     storage precision of Stata's own saving() file.  The draws match.
Repro: probe4.do.  Claim upheld.

## NEGATIVE (clean): string cluster() is refused informatively on every engine
`cluster(stabb)` (a string state abbreviation) -> r(109) "string variables not allowed in option
cluster(); stabb is a string variable", identically on the analytic, count-leg, reid-leg and
permute paths.  No crash reaches the re-identification leg's gegen min()/max() predicate.

## NEGATIVE (clean): reps()/permreps() cross-engine resolution
`bootstrap(...) reps() permreps()` -> r(198); `permute(...) permreps() reps()` -> r(198);
`reps()` with neither engine -> r(198).  No silent precedence.

---

## BR-3  Without cluster(), the two legs draw DIFFERENT resamples at the same rseed(),
##       so att's reported standard error moves when a second estimand is added
SEEDED (this is exactly the recon lead's question).  CONFIRMED, certain (mechanism), and the
statistical consequence is measured to be Monte Carlo only.

Mechanism, read off the pinned build and confirmed by measurement:
 - count leg   (ado:9295) : `bsample, cluster(<cluster or row id>) idcluster(__wdbns)`
 - physical leg(ado:9238) : `bsample, cluster(<cluster or row id>) idcluster(__wdbns)`  (same form)
 - reid leg    (ado:9345) : `if "`clus'" == "" bsample`  <-- a PLAIN bsample
The reid leg's own comment (ado:9340) says "same-seed draws equal the count path's".  That is
true WITH cluster() (probe1/cmp1: the w_* multiplicity columns are bit-identical) and FALSE
without it.

Observed (probe1.log/cmp1.log, reps(30), rseed(4242), 60 units x 8 periods, no cluster()):
  `att`         -> se(att) = 0.071372   (count leg)
  `att att_i`   -> se(att) = 0.061544   (reid  leg)
  the replicate draws themselves differ by up to 1.18e-01 relative.
Is it more than Monte Carlo?  probe8.do, reps(300), two seeds per leg:
  count seed A 0.081208 | count seed B 0.075653   (same-leg spread 7.3%)
  reid  seed A 0.077783 | reid  seed B 0.078920
  cross-leg at a common seed: 4.4% (seed A) and 4.3% (seed B) -- INSIDE the same-leg
  seed-to-seed spread.  So the two legs are the same estimator; only the draw sequence differs.
Expected: a user who reruns a specification with a second estimand added, at the same rseed(),
should get the same att column.  The fix is already in the file: the physical-expansion leg at
ado:9238 draws with `bsample, cluster(<row id>) idcluster(...)` on exactly this unclustered
geometry, and probe5 shows that leg's draws coincide with the count leg's to 3 ulps.  Only the
re-identification leg at ado:9345 uses a plain `bsample'.  They do under cluster(); they do not without it.  Nothing in
the help says the seed stops being comparable across estimand lists, and the code comment says
the opposite.  The remedy already exists in the file: the physical leg uses the count leg's
bsample form on the same unclustered geometry.
Repro: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-reproduce/repro_br3.do
Self-contained run, 40 units x 8 periods, reps(200), rseed(12345):
  UNCLUSTERED  att         (count leg) se_att = 0.235921
  UNCLUSTERED  att att_i   (reid  leg) se_att = 0.231988   replicate draws differ up to 3.82e-01
  CLUSTERED    att         (count leg) se_att = 0.169566
  CLUSTERED    att att_i   (reid  leg) se_att = 0.169566   replicate draws differ by 2.61e-14
Severity medium-low, realistic, certain about the divergence, and measured to be MC in effect.

## BR-4  The bootstrap resamples the `smp' universe, not e(sample): phantom clusters are drawn,
##       e(boot_df) counts them, and the documented match with Stata's own prefix fails
SELF-FOUND.  CONFIRMED, certain about the facts.

`wooldid_bootrun` builds its resampling universe as `smp' = (user `if') & nonmissing cluster
(ado:8915-8925), then resamples THAT.  Rows and whole clusters that the estimation drops --
missing outcome, missing control, singleton, outside the event-study window -- are still in the
resampling universe.  Stata's own `bootstrap` prefix, which sthlp:532 says wooldid's count leg
"matches ... replication for replication under a common seed", resamples e(sample).

Observed (probe13.log; 60 units x 8 periods, 15% of outcomes missing at random, and unit 5's
outcome missing entirely -- an ordinary administrative-panel shape):
  plain analytic run     : e(N) = 397, e(N_clust) = 59   (unit 5 is NOT in the estimation)
  bootstrap(standard)    : e(N) = 397, e(boot_df) = 59
                           -> the default bootdf(clusters) is G-1, so the engine used G = 60
  bootsave()             : 60 w_ columns, INCLUDING w_5 for the unit with no estimable data;
                           its mean multiplicity over the replications is 1.2, i.e. it is drawn
  Stata's bootstrap prefix on the same data/seed: e(N_clust) = 59
  wooldid's draws vs the prefix's: MAXRELDIF = 1.07e-01, ALL 40 replications differ.
  On the SAME design with no missing outcomes (probe4.log) the two agree to 2.6e-08, which is
  the float precision of Stata's own saving() file -- i.e. they match exactly.
Expected: (a) the help's replication-for-replication claim should hold, or should be qualified;
(b) `bootdf(clusters)` should be G-1 over clusters that contribute to the estimation, not over
clusters that merely have a nonmissing cluster id; (c) a w_ column for a cluster with no
estimable data is a phantom in the saved draw file.  Each replication draws G = 60 clusters of
which a random number are phantoms, so the effective number of contributing clusters varies
replication to replication instead of being fixed -- a small distortion here, larger the more
units lack outcome data.
probe14.do ISOLATES THE CAUSE with a two-arm test on the same data and seed:
  ARM A, the all-missing-y unit left in the data : e(N_clust)=59, e(boot_df)=59 (so G=60),
         wooldid vs Stata's prefix MAXRELDIF = 1.0678e-01  (they disagree on every replication)
  ARM B, the same unit physically dropped        : e(N_clust)=59, e(boot_df)=58 (so G=59),
         wooldid vs Stata's prefix MAXRELDIF = 2.6110e-08  (float precision -- they MATCH)
  ARM C, a cluster variable with missing values  : rows with a missing cluster ARE excluded
         correctly (e(N_clust)=56), but the missing-outcome unit is still counted (boot_df=56
         implies G=57), so the phantom is specifically the estimation-dropped cluster.
Note that in ARM A the SAME e() carries e(N_clust)=59 and e(boot_df)=59, i.e. the posted t
reference has one degree of freedom more than the run has clusters.  The same inflated count
feeds wooldid_bootctl's "at least 2 clusters" refusal and its "fewer than 10 clusters" and
"fewer than 6 treated clusters" cautions, which are computed on the same `smp' universe --
so a run with 9 contributing clusters and 2 all-missing ones would not get the few-cluster
caution at all.


### BR-4 ESCALATION -- the same estimation sample, two different inferences (repro_br4.do)
14 units x 8 periods, 5 treated 2003 / 4 treated 2005 / 5 never; FIVE units (three treated, two
never) have NO usable outcome at all.  ARM A leaves them in the data with a missing outcome;
ARM B drops them.  BOTH ARMS ESTIMATE THE IDENTICAL SAMPLE: e(N) = 72 and e(N_clust) = 9 in each.

  ARM A (phantoms present) : e(boot_df) = 13   -> the engine used G = 14
                             NEITHER few-cluster caution is printed
  ARM B (phantoms dropped) : e(boot_df) =  8   -> G = 9  (correct)
                             "Caution: only 9 clusters are available for the bootstrap ..."
                             "Caution: only 4 treated cluster(s) are available for resampling ..."

So on identical estimation samples the posted t reference distribution is t(13) instead of t(8).
t(13,.975) = 2.160 against t(8,.975) = 2.306: 6.8% narrower.  Precisely what that moves, per
sthlp:1340: the Wald family (e(boot_ci_normal), e(boot_wald_mainresults)) and everything that
reads e(V) with e(df_r) -- a hand `test' or `lincom', esttab, coefplot, the honestdid handoff --
and, under `bootdisplay(wald)', the displayed interval and p-value themselves.  The default
percentile columns are rank-based and do not move.  Both of the few-cluster cautions
that sthlp:532 promises ("A caution fires when fewer than 6 treated clusters are available for
resampling, beside the existing caution at fewer than 10 clusters in total") are suppressed.
Presence of an unusable unit in the data file, which changes nothing about what is estimated,
changes the inference and silences the warnings that the design is fragile.

Severity HIGH: this is an inference path that does not deliver what it claims (too-narrow
intervals, too-small p-values) plus two documented safety cautions going silent, on the most
ordinary shape of applied data there is.  Certain.
Primary repro: repro_br4.do (self-contained, no generator).  Corroborating: probe13.do, probe14.do.

## OPEN / IN PROGRESS
- (closed) probe8: the unclustered leg gap is Monte Carlo, not a different estimator. See BR-3.
- (closed) probe11c: e(gteffects) cell columns agree between the count and reid legs to
  1.94e-14 (nesting regime, 144 cells) and 3.99e-15 (keep-i regime); the subgroup pairwise
  contrast row agrees to 1.70e-15.  wooldid_bootcells' new-i -> old-i remap is correct.
- (closed) probe9: `att_i` draws are IDENTICAL whether att_i is asked alone or with att;
  att draws agree to 6.6e-16 whether att is asked alone or with all four estimands;
  an `if'-restricted run's count and reid legs agree to 1.3e-15.
- (closed) probe3b: a FRESH Stata session reproduces probe3's numbers bit for bit.
- (closed) probe13: bootsave()/permsave() on a pre-existing file without `replace' both stop
  the run with r(602) before any replication.

## CAVEAT for verifiers
cmplib.do uses `confirm variable`, which accepts unique ABBREVIATIONS.  In probe9.log the line
"CMP|(i) att_i alone vs att+att_i|main_att draws MAXRELDIF=6.0138e-02" is SPURIOUS: the
att_i-only file has no main_att column, and `confirm variable q1main_att` matched
`q1main_att_i`.  Every other comparison in cmp1.log/cmp3.log compares exactly-named twins.
