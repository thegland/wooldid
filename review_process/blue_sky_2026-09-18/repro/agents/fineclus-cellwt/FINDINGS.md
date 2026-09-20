# fineclus-cellwt — FINDINGS (gate archaeology on the refusals T1/T2/T3 removed)

All paths absolute. Build: .../bluesky/build/wooldid.ado (v1.4, sha c1d8979).

---

## GATE INVENTORY (what was removed, and what stands in its place)

| old site | old rc | refused | pinned build |
|---|---|---|---|
| ado:5990 permctl | 459 | `permute()` with cluster not nesting `i` | **STILL LIVE at ado:6078**, conditioned on `"`clus'" != "`ivar'" & "`onc'" == ""` — lifted only under `permonsetcohorts` |
| ado:6353 bootctl | 198 | `bootstrap(standard)` + cell-weighted estimand + **no** `cluster()` | GONE |
| ado:6364 bootctl | 459 | `bootstrap(standard)` + cell-weighted + cluster not nesting `i` | GONE |
| ado:6365 bootctl | 459 | `bootstrap(standard)` + `boottstat` + `cluster()` not nesting `i` | GONE |

Replacement (ado:8987-8998 and 9341-9350), verbatim:
```
local wdb_nest 0
if "`reid'" != "" & "`clus'" != "" { ... assert cmin==cmax by(i) ... wdb_nest = (_rc==0) }
if "`clus'" == "" bsample
else bsample, cluster(`clusv') idcluster(__wdbns)
if `wdb_nest' qui gegen double __wdbnewi = group(__wdbns `i')
else          qui gen   double __wdbnewi = `i'
```
ado:6451 states the policy: "Nothing is screened in its place, deliberately."

---

## F-A (CONFIRMED, deterministic).  Under the keep-i regime the estimand weights are
## blind to cluster multiplicity; only cluster PRESENCE survives.

repro: /private/tmp/.../probes/fineclus-cellwt/probe1b_mechanism.do  (log alongside)

The two regimes were built BY HAND on one 240-row panel — industry 1's rows appended
once with the original unit ids (what `__wdbnewi = i` produces) and once with fresh
unit ids (what `group(__wdbns, i)` produces).  Measured:

```
BASE         att=  1.0934572399  att_i=  1.1815039002  att_it=  1.0934572399
KEEPI_DUP1   att=  1.0154662125  att_i=  1.1437911577  att_it=  1.0557390439
REID_DUP1    att=  1.0154662125  att_i=  1.1026004128  att_it=  1.0154662125
KEEPI_DUP5   att=  1.0039783501  att_i=  1.1896036901  att_it=  1.0993223580
REID_DUP5    att=  1.0039783501  att_i=  1.0745092129  att_it=  1.0039783501
```
`att` is BIT-IDENTICAL between the two regimes on both duplications, so the stage-1
regression is invariant to the re-identification: the ONLY thing the branch changes is
the estimand aggregation weight.  `att_i` moves 1.1438 -> 1.1026 on one duplicated
cluster (3.6%) and 1.1896 -> 1.0745 on five (10.7%).  Under keep-i a twice-drawn
cluster's units count ONCE; under re-id they count twice.

This is exactly the collapse the design note at ado:9334-9337 names as the reason the
nesting regime must NOT keep `i`.  It is unavoidable and correct when the cluster is
genuinely finer than `i`.  It is NOT correct when the cluster is coarser than `i` for
all but one unit — see F-B.

---

## F-B (headline; measurement running).  One straddling unit silently flips the whole
## sample from the multiplicity bootstrap to the presence bootstrap.

The predicate is ALL-OR-NOTHING over the sample: `assert cmin==cmax if smp` fails if a
SINGLE unit carries two cluster codes.  A firm that reclassifies its industry code
mid-panel, a person who moves state, a hospital that changes referral region — each
flips a 40-cluster coarse-cluster bootstrap into the keep-i regime for every unit in
the sample.  Old behaviour: r459, refused.  New behaviour: silent, and by F-A the
estimand's cluster-multiplicity variance is gone.

repro: probe1_straddle.do (arms A = nesting, B = one cell changed, C = the straddler
dropped).  RESULTS PENDING.

---

## F-C (CONFIRMED by source archaeology; runtime check pending).  The reps() resolution
## caution is unreachable without cluster() — and T2 just made that the third regime.

`git show 813c414` (RF-2, 2026-09-13) moved the `if `reps' > 0 { ... }` resolution
caution from an unconditional site (then :6011) to inside `if "`clus'" != "" {`
(pinned build :6459 opens, :6501 closes; the caution sits at :6473-6478).  So
`bootstrap(standard) reps(10)` with no `cluster()` prints nothing, while the same run
with `cluster(unit)` prints "Caution: with reps(10) the coefficient p-values fall on a
grid of 0.200 ...".  ado:6489-6490 shows the author knew the unclustered leg is not
covered by the few-treated-cluster caution; there is no such note for this one.

---

## F-D (CONFIRMED by source).  The bootstrap regime has no record anywhere.

`wdb_nest` is `c_local`-ed back to the shell (ado:8998) and used only to decide whether
to print the A10 covariate caution (ado:407).  It is not in `e()`, not displayed, and
not in `bootsave()`.  T1's counterpart IS recorded — `e(perm_irebuilt)` (ado:195),
kept expressly as "the record" when the author declined a display line (2026-09-16).
So a user cannot tell which of two different variance concepts produced their SE.

---

## F-E (CONFIRMED).  The mechanism of the 149/150 over-coverage — and the reason the
## record gives is the wrong one.

repro: probe3_overcov.do   (p3.csv; p3_partial_reps199.csv is the reps(199) prefix)

`test/suite/mc/mc_fineclus.do` records P5 = att_i covered 149/150 = .9933 at
cluster(state) == i, and explains it as "the same effect the help already describes
under 'What level to cluster at': at or above the i level the covariance stops seeing
the treated cells' own variation.  It lands conservative here."  That cannot be the
mechanism: the effect it names is the documented ANALYTIC understatement, measured at
SE/SD 0.66-0.75 (sthlp:255) — it makes intervals too NARROW.  Over-coverage is the
opposite sign.

The actual mechanism is an ADDED term.  At cluster() == i the re-id leg resamples
cohorts and re-identifies each drawn copy, so the bootstrap distribution of att_i
carries a cohort-COMPOSITION component Var_s(tau_s)/S_treated.  The MC that produced
149/150 holds tau_s FIXED — in its DGP `tau = mult*0.2*nu` is a deterministic function
of the state index, redrawn never — so the sampling distribution it measures contains
no such component.  The bootstrap is estimating the variance of a
superpopulation-of-cohorts estimand while the MC measures a fixed-cohorts one.

Measured, on a re-implementation of that DGP with the noise scaled to 2% (so the true
sampling SD of att_i is ~0.003) and cluster(state) == i, reps(99):

```
                        b(att_i)              se(att_i)
 HET  tau = .4 / 1.2   0.801175892637377      0.0868634699692049
 HOM  tau = .8 for all 0.801175892637377      0.0019218211999355
```
Same seed, same noise draws, same att_i truth (0.8), same cluster count: the ONLY
difference is Var_s(tau_s) = 0.16 versus 0.  The bootstrap SE is 45x larger in HET.
Predicted composition SE sqrt(0.16/20) = 0.0894; measured 0.0869.

The keep-i arm on the same data (cluster(unit), strictly finer) gives 0.0244 in HET
and 0.0022 in HOM, i.e. it carries 28% of the composition SE.  NOTE WHAT THAT MEANS
FOR P4: keep-i is NOT composition-free on this design.  Treated states hold 2 or 6
units, and a 2-unit state goes entirely undrawn with probability (1-2/160)^160 = e^-2
= .135, so states DO vanish from replications and some composition variance survives.
P4's 142/150 = .9467 is therefore not evidence that the keep-i rule is correct; it is
evidence that on this design a presence term worth 28% of the composition SE happened
to land the total near the fixed-cohort sampling SD (ratio 0.963).  A reader must not
cite P4 as a clean bill of health for the keep-i branch.

The two arms are not two implementations of one variance; they are two different
variance concepts, and which one a run gets is decided silently by the nesting
predicate.

WHAT THIS MEANS FOR THE READER OF A wooldid TABLE.  Under bootstrap(standard) at
cluster() == i, att_i's interval covers the SUPERPOPULATION att_i (over cohorts drawn
like these) and not the att_i of the cohorts in the data.  Neither the table nor the
help says which.  A user who reads it as the second gets .9933 coverage at a nominal
.95, and one who reads a cluster-finer-than-i run the same way gets a different object
again.

HOW F-E BEARS ON F-B AND F-F, stated so nobody over-reads either.  Whether the
composition term BELONGS in att_i's standard error is a question about the estimand -
fixed cohorts or a superpopulation of them - and that is a ruling for the author, not
something this probe settles.  The severity of F-B and F-F does NOT rest on "the
nesting SE is the correct one".  It rests on this: two datasets that differ in one
discarded cell, with identical estimation samples, identical point estimates and
bit-identical att standard errors, get att_i standard errors 23% apart, because a
predicate that tracks neither question silently picks between the two variance
concepts - and no e(), no display line and no help sentence says which was picked.

### F-E, the quantitative closure (probe3, full-noise arms)

DGP as re-implemented above, 40 states, T=6, reps(99), 5 bootstrap draws per arm,
Monte Carlo sampling SD from 70 fresh draws.  HET and HOM are PAIRED: the same noise
draws, and b(att_i) is bit-identical between them, so only Var_s(tau_s) differs.

```
  MC sampling SD of att_i (the truth coverage was judged against)      0.13328
  NEST cluster(state)==i  boot SE  HET  RMS 0.15501   ratio to SD  1.163
  NEST cluster(state)==i  boot SE  HOM  RMS 0.12167   ratio to SD  0.913
  FINE cluster(unit)      boot SE  HET  RMS 0.12830   ratio to SD  0.963
  implied composition SE, paired HET-HOM: .0903 .1343 .0675 .0655 .1053
                                          (predicted sqrt(.16/20) = .0894)
```
Those ratios predict BOTH published numbers:
  nominal .95 at ratio 1.163 -> 2*Phi(1.96*1.163)-1 = .977   observed P5 149/150 = .9933
  nominal .95 at ratio 0.963 -> 2*Phi(1.96*0.963)-1 = .941   observed P4 142/150 = .9467
(149 against Binomial(150,.977) is +1.3 sd; 142 against Binomial(150,.941) is -0.4 sd.)

The HOM arm is the one that carries the effect the record NAMED: with the composition
term switched off, the cluster()==i bootstrap SE is 0.913 of the sampling SD, i.e. it
DOES understate, in the documented direction.  P5 over-covers because the composition
term is larger than that understatement, not because of it.

---

## F-F (CONFIRMED, and sharper than F-B).  The PHANTOM STRADDLER: a row the estimator
## discards decides which variance the bootstrap computes.

repro: probe1d_phantom.do  (p1d.csv)

wooldid_bootrun decides the regime on `smp' (ado:8910-8916), which is ONLY the user's
`if' (==1) plus a non-missing cluster value.  It does not screen rows the estimator
drops - a missing outcome, a zero weight, a period outside the usable window.  The
predicate (ado:8991-8993) then runs `if `smp'' on that wider row set.

Arms A and D differ in exactly one thing: the industry code recorded on unit 1's t=3
row, whose OUTCOME IS MISSING IN BOTH ARMS and which therefore never enters
estimation.  Measured (reps(199), 10 industries, 30 units, T=6):

```
 arm nest  e(N)          b(att_i)              se(att_i)            se(att)
  A   1    179   2.755071041835752   0.4348901884075443   0.4338581068233801
  D   0    179   2.755071041835752   0.3308192001858818   0.4338581068233803
  A   1    179   2.755071041835752   0.4724939276687015   0.4715096894625723
  D   0    179   2.755071041835752   0.3655332376168671   0.4715096894625724
  (second dataset, same pattern: .43498/.47448 -> .33027/.36696)
```
e(N) is identical, b(att_i) is identical to the last digit, and se(att) is identical
to FIFTEEN significant figures - so the resamples themselves are the same draws, the
same fits, the same everything.  The only thing that changed is `__wdbnewi', and
se(att_i) falls 22.6-24.0%.

ado:8980-8986 argues the `if `smp'' is load-bearing because evaluating the predicate
on the raw frame "instead" could misclassify a coarser cluster and "silently take the
keep-i rule and collapse att_i's composition variance".  `smp' is not the estimation
sample, so that exact misclassification is still reachable - it just needs the
straddling row to be in the design sample and out of the estimation sample, which a
missing outcome does.  sthlp:431 already documents, for permute(), that rows "wooldid
later removes by its deterministic markout, e.g. missing outcomes or zero weights,
still count toward the design"; nothing says it for the bootstrap regime.

---

## F-C (CONFIRMED, runtime + source).  REGRESSION: the reps() resolution caution is
## unreachable without cluster(), and T2 just made that a recommended configuration.

repro: probe2a_caution.do  (probe2a_caution.log)

Three paired runs on one dataset, reps(10) throughout:
```
 RUN1  no cluster(), att att_i att_it   rc 0   NO caution printed
 RUN2  + cluster(unit), same otherwise  rc 0   "Caution: with reps(10) the coefficient
                                                p-values fall on a grid of 0.200 ..."
 RUN3  no cluster(), att only           rc 0   NO caution printed
```
Cause, from `git show 813c414` (RF-2, 2026-09-13): the `if `reps' > 0 { ... }` block
was moved from an unconditional site (then :6011) into the body of
`if "`clus'" != "" {` (pinned build: opens :6459, closes :6501, caution at :6473-6478).
RUN1 is also the configuration old ado:6353 refused with r198 until 2026-09-16, so the
combination "newly permitted AND silently uncautioned" is one day old.
Severity medium: a missing caution, not a wrong number - but it is the caution that
tells a user a printed p of 0 is the resolution floor, and the unclustered leg is the
one where `bootdf(clusters)` is refused precisely because the resampling unit is the
observation (ado:6328-6331), so the authors have thought about this leg elsewhere.

---

## F-D (CONFIRMED, source).  Nothing anywhere records which bootstrap regime ran.

`wdb_nest` is computed in wooldid_bootrun (ado:8987-8998), `c_local`-ed back, and used
at ado:407 only to decide whether to print the A10 covariate caution.  It is not in
`e()`, not in the displayed footer, not in `bootsave()`.  Its T1 counterpart IS
recorded - `e(perm_irebuilt)` (ado:195, sthlp:1352) - and was kept expressly as "the
record" when the author declined a display line on 2026-09-16 ("Not to be reopened").
Given F-B and F-F, a user cannot tell from any output which of two variance concepts,
differing by 23% here, produced the standard error they are reading.
Severity low-medium on its own; it is what makes F-B and F-F silent.

---

## F-G (CONFIRMED; scope note, not a defect).  The permute() gate T1 was written to
## open is still closed for the configuration T1's own motivation names.

repro: probe4_permgate.do  (probe4_permgate.log)

12 states x 4 counties, T=8; states 1-4 adopt at t=4, 5-8 at t=6, 9-12 never.  i=state
(the help's "coarse i"), cluster(county) (the help's "finer cluster"):
```
 P1  permute(shuffle), cluster(county)                      r459 "requires the
       cluster() variable to nest i"
 P2  the same + permonsetcohorts                            r459 "exactly one level of
       state per distinct onset value ... 384 in-sample observation(s) sit in onset
       groups that span more than one level of state"
 P4  bootstrap(standard), SAME configuration                rc 0
```
`review_process/feature_tickets_2026-09-15.md` opens T1 with "a user who follows both
pieces of advice - coarse i, finer cluster - cannot use permute() at all."  After T1
that user still cannot, whenever two cohorts share an onset year, which is the ordinary
staggered-adoption shape.  permonsetcohorts requires a BIJECTION with the onset
partition, so the remedy is to redefine `i` as the onset grouping - which changes the
estimand `att_i` averages over.  The help is not wrong (sthlp:68 and :431 both state
the bijection); the gap is that the ticket's motivating case is not covered and nothing
says so.  Severity low-medium, documentation/scope.
(P3, which would have shown the remedy working when i IS the onset grouping, errored
r411 on my own bad fixture - cohid coded 0 for never-treated, which wooldid requires to
be a positive integer.  Untested, not a finding.)

---

## CLEAN — probed hard, came back clean

### C-1.  The unclustered cell-weighted bootstrap (the removed r198 gate) gives an SE
### of the right size on an ordinary single-row-cell panel.

repro: probe2b_noclus_se.do  (p2b.csv).  Readout PRE-REGISTERED in the file's footer
before the first number landed: ratio in [0.85, 1.15] = clean, < 0.85 = understated,
> 1.15 = the F-E composition mechanism on cells.

40 units, T = 6, one row per (i,t) so every treatment cell is a singleton, iid errors,
bootstrap(standard) reps(299) with NO cluster().  Monte Carlo sampling SD from 120
fresh datasets; bootstrap SE from the first draws of the same sequence.

```
  att     MC SD = 0.23453 (n=120)   boot RMS se = 0.25436 (n=7)   ratio = 1.0845
  atti    MC SD = 0.22390 (n=120)   boot RMS se = 0.22790 (n=7)   ratio = 1.0179
  attit   MC SD = 0.23453 (n=120)   boot RMS se = 0.23396 (n=7)   ratio = 0.9976
```
All three land inside the pre-registered clean band.  What would have shown up: the
mechanism of F-A applied to cells - a cell drawn twice is two identical rows whose
cell mean is unchanged, so the equal-cell weighting counts it once and only cell
PRESENCE survives - predicts a ratio near 0.76 if the presence term were the whole
story.  It is not: cell presence and the regression's response to row multiplicity
evidently compose to about the right total here.  The removed r198 guard's stated
reason ("observation-level resampling cannot re-identify duplicated treatment cells")
is mechanically true and inferentially harmless on this shape.
CAVEAT, stated because it bounds the claim: this is the singleton-cell shape with iid
errors and no cluster().  It does not license the keep-i branch under cluster(), which
F-B/F-F measure separately and which is a different object.

### C-2.  The point estimates never move.
Across every arm of every probe - the one-cell straddle, the phantom straddle, HET vs
HOM - b(att), b(att_i) and b(att_it) are bit-identical between the regimes.  The
regime switch touches only the replication frame, never the observed fit.  A defect
that moved the reported estimate would have shown here and did not.

### C-3.  The observation-weighted att standard error is untouched by the regime flip.
se(att) ratio B/A = 1.0029 over six paired runs in probe1c, and bit-identical to 15
significant figures in probe1d.  This is the predicted signature (att carries cluster
multiplicity through duplicate ROWS, which both regimes produce) and it is also the
control that rules out "the two arms simply drew different resamples".

### C-4.  The three removed refusals all genuinely run.
Every configuration the three retired gates used to refuse returns rc 0 on the pinned
build: cell-weighted + no cluster() (probe2a RUN1, old r198), cell-weighted + cluster
not nesting i (probe1c arm B and probe1d arm D, old r459), and the same configuration
under bootstrap at i=state/cluster(county) (probe4 P4, old r459).  No crash, no
uninformative error, no failed replications.

---

## UNTESTED (relevant things I did not get to)

- **boottstat under the keep-i regime.** T3's own constraint is that any scheme must
  keep duplicate draws distinguishable at the cluster level or the E[m^2]/E[m] = 2-1/G
  (~sqrt(2)) per-replication SE inflation returns.  Under keep-i the subinstr still
  puts the rep VCE on `__wdbns`, so it should be intact, but I did not measure
  E[se_r^2]/se_obs^2 on a finer-than-i fixture.  This is the one leg of the three
  retired gates whose *replacement mechanism* I did not exercise directly.
- **The `if` / weight route into the phantom straddle.** F-F used a missing outcome.
  `smp` also ignores zero weights and out-of-window periods; a `wooldid ... if x==1`
  whose excluded rows carry a different cluster code is the same defect through a
  different door, and the `if` route interacts with BS-37's exact-one rule.
- **permute() P3**: whether permonsetcohorts + cluster finer than i actually completes
  when `i` IS a valid onset grouping (my fixture coded never-treated as 0, which
  wooldid rejects r411).  The remedy T1 shipped is therefore untested here.
- **att_itime** — I probed att_i and att_it; att_itime is the third cell-weighted
  estimand and takes the same `wdb_reid` leg, but I never ran it.
- **Whether a caution or e() marker would be cheap.** I did not price a fix; ado:6451's
  "nothing is screened in its place, deliberately" is a policy about refusals, not
  obviously about disclosure.
- **bootstrap(bayes)** under a cluster not nesting i: `wdb_reid` is keyed on
  `bootstrap == "standard"`, so the bayes leg never takes the re-id path at all and the
  cell-weighted gate there (ado:6437-6444) still fires.  Not probed.
