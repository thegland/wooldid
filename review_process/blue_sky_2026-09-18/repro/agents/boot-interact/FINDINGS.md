# boot-interact — FINDINGS
Agent: bootstrap() x everything else. Method: interaction matrix + covariance provenance tracing.
Build: pinned `.../bluesky/build/wooldid.ado` (v1.4, sha c1d8979), loaded with `adopath ++`.
All probe do-files and .log files are in this directory.

---

## BI-1 (HIGH) — `pretrendpower()` under `bootstrap()` prints an uncorrected k-dimensional
## Wald test on a covariance estimated from `reps` draws, referred to F(k, e(df_r))

**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-interact/p12b_hotelling.do`
(run: `stata-mp -e do p12b_hotelling.do`; log `p12b_hotelling.log`, completes in ~9 min)
(supporting: `p10_ptpboundary.do`, `p08_lowest.do`)

**TWO DEFECTS, TWO FIXES.** (i) the joint pre-trend p is referred to the WRONG reference
distribution — a df correction repairs it; (ii) the GLS slope, its CI and both
minimum-detectable-trend columns are different NUMBERS because S^-1 is a garbage estimate
when n is close to k — no df correction repairs those, only a guard keyed on
`e(boot_completed)` relative to the block dimension, or a refusal. The `bootdf()` table below
is the evidence that they are separate: the slope and its se are bit-identical across every
`bootdf()` setting while the p moves by eight orders of magnitude.

**THE CONTROL.** Under the ANALYTIC path the three quantities agree EXACTLY:
`e(roth_prep_att)` = `e(joint_horiz)` = hand `test` p = .01003909. F(k, e(df_r)) is the right
reference when `e(V_es_*)` is not a sample covariance. The bootstrap path broke a reference
that was correct.

**Claim.** `wooldid_ptp` (ado:10545) computes the joint pre-trend test as
`Ftail(keff, dfr, W/keff)` with `dfr = e(df_r)`. After `bootstrap()` the matrix W is built
from is `e(V_es_att)`, a SAMPLE covariance from n = `e(boot_completed)` replications, while
`e(df_r)` is the bootstrap reference df for a ONE-dimensional t (G-1 = 39 here). The
k-dimensional quadratic form needs Hotelling's reference,
`(n-k)/((n-1)k) * W ~ F(k, n-k)`; wooldid's version inflates the F statistic by
`(n-1)/(n-k)` — 6.0 at n=7,k=6; 3.5 at n=8; 1.56 at n=15; 1.31 at n=60,k=15.
The same covariance drives the GLS slope, its `level()` CI and both minimum-detectable-trend
columns, so those move too.

**Observed** (40 states x 20 periods, cluster(state), esfixedbaseperiod,
esprelength(6) espostlength(4) pretrendpower(80); k = 6 pre coefficients):

| run | KEFF | roth slope | slope CI | joint pre-trend p | mdt_wald |
|---|---|---|---|---|---|
| analytic | 6/6 | -.06623 | (-.1220, -.0105) | **.01004** | .1101 |
| reps(7)  | 6/6 | **+.08134** | (+.0667, +.0960) | **2.36e-12** | **.02885** |
| reps(8)  | 6/6 | +.10420 | (+.0542, +.1542) | 2.703e-06 | — |
| reps(10) | 6/6 | +.11512 | (+.0649, +.1654) | 1.788e-06 | — |
| reps(15) | 6/6 | -.04258 | (-.1167, +.0315) | .3654 | — |
| reps(25) | 6/6 | -.05956 | (-.1249, +.0058) | .2981 | — |
| reps(50) | 6/6 | -.04690 | (-.1076, +.0138) | .3299 | — |

At reps(7) the reported pre-trend slope has the WRONG SIGN, its 95% CI excludes zero on the
wrong side, the joint pre-trend p says 2e-12 where the settled value is ~0.3, and `mdt_wald`
says the design detects a trend 3.8x smaller than it actually can — a false reassurance in
exactly the column the diagnostic exists to provide.

**The only caution printed at reps(7)** is the generic p-value-GRID caution
(ado:6477, "with reps(7) the coefficient p-values fall on a grid of ..."), which is about
Monte-Carlo resolution of the percentile p and explicitly blesses low reps ("reps() has no
minimum - small values are a legitimate way to time a specification"). Nothing mentions the
Roth block, the covariance rank or the replication count behind it.

**The Roth rank guard is the wrong guard.** `wooldid_ptp` screens with
`diag0cnt(invsym(S))` on the k x k pre block. That fires only while completed <= k
(measured: reps(3) -> KEFF 2, reps(4) -> KEFF 3, reps(6) -> KEFF 5, each with the documented
rank-deficiency caution). At completed = k+1 it goes silent and stays silent, which is
precisely where the Hotelling factor is largest.

**Documentation contradiction.** sthlp:532 says of the bootstrap joint tests:
"Neither is the Wald test on e(V) that a hand `test` reproduces; that one is a normal-theory
test on the bootstrap covariance and **is not printed**." `pretrendpower()` prints exactly
that statistic, as "joint pre-trend test: F(6, 39) = ..., p = ...", under the bootstrap
header. The same entry also says "the existing rank guard applies to the bootstrap block as
it does to the analytic one, skipping with a caution if it is not invertible" — true, and not
the relevant guard.

**Internal contradiction in one run.** With `jointtests pretrendpower(80)` at reps(8) the run
prints two bootstrap pre-trend p-values that differ by five orders of magnitude:
`e(joint_horiz)` = .1111 (the rank test on the draws, immune) and `e(roth_prep_att)` =
2.703e-06 (the uncorrected Wald). At reps(400) they agree (.4838 vs .4928).

**Expected.** Either the Hotelling/`F(k, n-k)` reference, or a guard keyed on
`e(boot_completed)` relative to the dimension of the tested block, or a refusal — the same
treatment `jointtests` already gets. The help's own sentence says this statistic is not printed.

**Mechanism CONFIRMED numerically** (p12_hotelling.log). The hand `test` p and
`e(roth_prep_att)` are identical to every printed digit at every reps — they are one
statistic. Under the ANALYTIC path all three agree exactly
(roth_prep = joint_horiz = hand test = .01003909), which is the control: there the
F(k, e(df_r)) reference is the right one. Under bootstrap they diverge:

| reps n | W | wooldid roth_prep = hand test | Hotelling p, F(k, n-k) | wooldid jointtests (rank) | (n-1)/(n-k) |
|---|---|---|---|---|---|
| analytic | 19.81 | .01003909 | — | .01003909 | — |
| 7  | 159.73 | **2.362e-12** | .3482 | .1250 | 6.00 |
| 8  | 55.30  | **2.703e-06** | .3006 | .1111 | 3.50 |
| 10 | 57.47  | **1.788e-06** | .0912 | .0909 | 2.25 |
| 15 | 6.756  | .3654 | .6418 | .3125 | 1.56 |
| 25 | 7.559  | .2981 | .4553 | .2692 | 1.26 |
| 50 | 7.163  | .3299 | .3937 | .3333 | 1.14 |
| 200 | 5.494 | .4941 | .5013 | .4527 | 1.03 |

The Hotelling-corrected value tracks wooldid's own (immune) rank test; the printed value
does not.

**`bootdf()` proves the point** (p14_ptpvariants.log, reps(8), identical data and seed).
The Roth slope and its standard error are IDENTICAL across every `bootdf()` setting
(slope .10420368, se .02469802) — the covariance is the same object. Only the reference
distribution moves, and with it the printed joint pre-trend p:

| variant | e(df_r) | reference | roth_prep |
|---|---|---|---|
| default `bootdf(clusters)` | 39 | F(6, 39) | 2.703e-06 |
| `bootdf(normal)` | . | chi2(6) | **4.023e-10** |
| `bootdf(3)` | 3 | F(6, 3) | .04797 |
| `boottstat` | 39 | F(6, 39) | 2.703e-06 |
| (Hotelling, n=8) | — | F(6, 2) | .3006 |

So a user option documented as choosing the reference distribution for the COEFFICIENT
statistics silently moves `pretrendpower()`'s joint pre-trend p by eight orders of magnitude,
and `bootdf(normal)` — the default whenever `cluster()` is absent — is the worst case.

**Ordinary-reps version** (p09_calib.log): `esprelength(16)`, 11 realized pre coefficients,
`reps(60)`: `e(roth_keff_att)` = 11 of 11 (no rank caution — and note the ANALYTIC run on the
same data has keff 9 of 11 and DOES print one), W = 10.542, wooldid prints p = .4984, the
Hotelling reference gives .6431, inflation (60-1)/(60-11) = 1.20. **Zero caution lines are
printed in that run** (2/60 = .033 < .05). A 20-period event study with reps(60) is ordinary
practice, so the defect is not confined to toy replication counts.

**realistic:** true
**severity:** high  **confidence:** certain (mechanism), probable (that Hotelling is the exact
correction rather than an approximation — b_obs is not independent of S)
**surface:** `bootstrap()` x `pretrendpower()`; `wooldid_ptp` ado:10342-10560

---

## BI-2 (MEDIUM) — the posted bootstrap `e(V)` is singular whenever completed reps <= the
## number of reported coefficients, and nothing records it; the help invites a hand `test`

**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-interact/p07_testrank.do`
(supporting: `p01b_rank.do`, `p03_conseq.do`)

**READ THE CALIBRATION FIRST.** The arithmetic below is GENERIC to bootstrap + Wald in Stata
(measured: Stata's own `bootstrap` prefix does the same thing, p13_calib2.log). wooldid's
specific exposure is three things and nothing else: `reps()` has no floor while `permreps()`
has one of 20; wooldid's own tests are guarded while the hand route is not; and sthlp:1340
actively promotes the hand route. This is a documentation and asymmetry finding, not a
computation bug.

**Claim.** `wooldid_bootwire` reposts `e(V)` with `ereturn repost V=`, deliberately bypassing
Stata's positive-definiteness screen (ado:10037-10048; the comment says so and says downstream
consumers "must screen for a non-positive-finite diagonal first"). `rank(e(V)) =
min(e(boot_completed) - 1, k)`. Nothing in `e()` or in any caution reports the rank; the only
low-reps caution is keyed on `2/reps > .05` (reps < 40), a different and unrelated condition.

**Observed** (k = colsof(e(b))):
- 40x12 panel, k = 12: reps(8) rank 7, min eig -6.5e-17; reps(19) rank 12; reps(21)/40/200 rank 12.
- 40x20 panel, k = 14: reps(4) rank 3; reps(8) rank 7; reps(12) rank 11; reps(30)/400 rank 14.
- `e(V_es_att)` inherits it: 12 horizons, reps(8) -> rank 7, reps(12) -> rank 11.

Hand `test` over the six pre-period coefficients on the SAME runs:
analytic .01004 | reps(4) .8568 (df 3, drop 1) | **reps(8) 2.703e-06 (df 6, drop 0)** |
reps(12) 6.988e-05 (df 6, drop 0) | reps(30) .2782 | reps(400) .4928.
At reps(4) Stata detects and drops constraints; at reps(8) and reps(12) it does not
(`drop = 0`, full df, no message) and returns the anti-conservative value. wooldid's own
guarded joint test on the same coefficients says .1111 and .0769 there.

`honestdid` on the rank-7-of-12 bootstrap `e(V_es_att)` runs at rc 0 (p03_conseq.log).

**The realistic, zero-caution case** (`p09_calib.log`; trimmed self-contained twin
`p09b_zerocaution.do`, which was still running when the budget closed): `esprelength(8) espostlength(8)
subgroups(grp)` on the same 40-state panel gives `colsof(e(b))` = 72 and 64 event-study
coefficients — an ordinary long event study with three subgroups. With `reps(50)`:
`rank(e(V))` = **49 of 72**, min eigenvalue -5.9e-16, and **not one caution is printed**
(2/50 = .04 < .05, so the ado:6477 grid caution does not fire; grep of the run block returns
0 lines containing "Caution"). A hand `test` over the 64 event-study coefficients returns
**p = 9.948e-33 on df 49** with `drop = 1`: Stata silently dropped 15 constraints and
reported a p-value anyway. (Footnote, not a comparator: the analytic run on the same 64 coefficients gives
p = 4.681e-17 on df 36, i.e. 28 constraints dropped — the analytic `e(V)` is `diag(se^2)`
plus spliced event-study blocks with zero cross-family off-diagonals by the PS-62 ruling, so
it is degenerate for a different, already-ruled reason. Neither side is a clean reference for
a 64-constraint hand test; that is itself the point.)

**Expected.** sthlp:1340 says "e(df_r) is aligned with the reference distribution actually
reported, so a `test` or `lincom` run by hand after a bootstrap fit reproduces the
bootdisplay(wald) p-values under every bootdf() setting", and names esttab/coefplot/honestdid
as readers of e(V). That is true for a one-dimensional t and silently false for a
k-dimensional Wald once reps is small relative to k. A recorded `e(boot_V_rank)`, or a
caution when `e(boot_completed) <= colsof(e(b))`, would close it.

**Calibration note (honest), MEASURED** (p13_calib2.do). `bootstrap, reps(8) cluster(state):
regress y i.t` on a 40x20 panel posts `e(V)` of rank 7 of 21 and `test` over six time dummies
returns **p = 2.883e-81, df 6, drop 0, no warning** (reps(15) .001774, reps(50) .026803,
reps(400) .093996; analytic cluster-robust .162645). So the arithmetic is entirely generic to
bootstrap + Wald in Stata and is NOT a wooldid defect.
What is wooldid-specific is only: (a) `reps()` has no floor and low reps are explicitly
blessed, while `permreps()` has a hard floor of 20 on the same argument; (b) wooldid's own
tests are guarded against exactly this while the hand route is not; and (c) sthlp:1340
actively promotes the hand route without qualifying it to one dimension. The difference from
Stata's prefix is that Stata never PRINTS a multi-dimensional Wald on its bootstrap
covariance — wooldid does, under `pretrendpower()`. That is BI-1, and it is the real defect;
BI-2 is the documentation half of it.

**realistic:** true  **severity:** low-medium (the arithmetic is Stata's; the doc invitation and the
permreps()/reps() floor asymmetry are wooldid's)  **confidence:** certain
**surface:** `bootstrap()` `reps()`; `wooldid_bootwire` ado:10037-10048

---

## BI-3 (see p09) — the ordinary-reps version of BI-1 (pending)
`esprelength(16)` (15 pre coefficients) with `reps(60)`: 2/60 = .033 < .05 so NO low-reps
caution fires at all, and the inflation factor is (60-1)/(60-15) = 1.31.

---

## PROVENANCE BASELINE — what is CORRECT (negative results)
- The seeded headline lead ("pretrendpower runs Roth's formula on a covariance that is NOT
  the one the displayed SEs come from") is **guarded in code and confirmed clean by run**:
  `wooldid_bootwire` swaps `e(V_es_*)` to the bootstrap block and preserves the analytic as
  hidden `e(V_es_analytic_*)`; `wooldid_ptp` (ado:10427) gates on the presence of that twin
  and skips with a caution when the swap failed. Measured after every bootstrap run:
  `e(wdid_vcvstate)` = "bootstrap", `e(V_analytic)` and `e(V_es_analytic_att)` both present,
  `mreldif(e(V_es_att), e(V_es_analytic_att))` = .017-.153 (genuinely swapped),
  `e(boot_ves_failed)` empty. At reps(1) the swap correctly FAILS and the caution fires.
- `e(df_r)` IS reposted to the bootstrap df (ado:10225-10246) before ptp reads it at ado:449,
  with `e(df_r_analytic)` kept. Provenance of the scalar is right; BI-1 is about its REUSE.
- The bootstrap joint tests are a genuine statistic-bootstrap on a consistent metric:
  W_r from the replication's OWN analytic `e(V_es_*)` (`wooldid_bootjwr`, ado:15238), W_obs
  from the OBSERVED analytic twin (`wooldid_bootjwobs`, ado:15311), both behind a scale-free
  eigen guard. At reps(400) the hand Wald and the rank test agree (.4928 vs .4838).
- `wooldid_bootves` (ado:14582) rebuilds `e(V_es_*)` with the analytic matrix's own row and
  column stripes and refuses (`wdb_vesok=0`) on any unmapped name, so there is no silent
  misalignment between `e(b_es_*)` and `e(V_es_*)` for ptp's positional indexing.
- `reps(1)`: completed 1, `e(boot_N_V)` 0, every reported se/t/p/CI missing, the V_es swap
  declines with a named caution. Correct observable failure, not a defect.

---

## BI-4 (LOW-MEDIUM) — `bootstrap()` + `timetrends` + `cluster()` at the level of `i`:
## the table is rescued, `e(V)` is not, and the help says the opposite

**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-interact/p15_timetrends.do`
(first seen as cell 19 of `p06a_matrix.do`)

**Claim (the failure is LOUD, not silent).** On the design sthlp:1059 documents (cohort trend terms with `cluster()` on `i`),
the estimation command returns an all-zero covariance and wooldid posts the estimate WITHOUT
a `V`. Under `bootstrap()` the replications still vary, so wooldid computes and PRINTS a
complete bootstrap table — se, t_wald, p and CI on every row. But `ereturn repost V =` then
fails **rc 503**, so `e(V)` never becomes the bootstrap covariance, `e(wdid_vcvstate)` reads
`absent`, and `e(V_analytic)` is never created. Every downstream reader — `test`, `lincom`,
`esttab`, `coefplot`, `honestdid` — fails (loudly, rc 321) on a run whose own output is fully
populated with valid bootstrap inference.

**Observed** (p15_timetrends.log; 40 states x 12 periods, cluster(state), esfixedbaseperiod,
esprelength(3) espostlength(3), timetrends, reps(60)): `e(wdid_vcvstate)` = "absent",
`e(boot_v_posted)` = 0, `confirm matrix e(V)` -> rc 111, `e(V_boot)` -> rc 0,
`e(V_analytic)` -> rc 111.  The printed table is fully populated with bootstrap quantities
(at reps(20) on the same design: `main att` = -1.5685, se .74513, t -2.11, p 0.000).  A hand
`test main_att = 0` and `lincom main_att` both fail **rc 321** - a LOUD, observable failure,
not a silent wrong number.  Caution printed: "the bootstrap covariance could not be posted to
e(V) (rc 503); it remains available as e(V_boot). Reported tables carry the bootstrap
standard errors either way." 

**Expected.** The run-time caution is honest and names `e(V_boot)`, so the code discloses it.
sthlp:1059 does not: it states flatly that on this design "every standard error, test
statistic, p-value and interval in the run is missing" and that `e(wdid_vcvstate)` reads
`absent` "so honestdid, pretrendpower(), the joint tests and simultaneousbands have no
input" — with no carve-out for `bootstrap()`, under which the first clause is false and the
second is true for a different reason (the repost failed, not the analytic fit).

**realistic:** true  **severity:** low  **confidence:** certain
**surface:** `bootstrap()` x `timetrends` x `cluster(i)`; `wooldid_bootwire` ado:10048

---

## BI-5 (LOW) — the bootstrap simultaneous band's critical value is the max of the
## completed draws at low reps, and only the p-grid caution mentions reps

**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-interact/p03_conseq.do` (the reps(2000) arm was
killed for slot budget; reps(8)/reps(12) and the analytic arm completed)

**Observed** (40 states x 20 periods, esprelength(6) espostlength(6), simultaneousbands):
analytic simulated band critical value 2.6504; `reps(8)` 2.0712 with `e(esband_n_att)` = 8;
`reps(12)` 2.4471. A 95% sup-t quantile taken over 8 draws IS the maximum of 8 draws, so the
band is structurally too narrow. On the other design (p09, 11 pre coefficients) analytic
2.6920 vs `reps(60)` 2.7353 — by reps(60) the band is fine.
`e(esband_n_att)` records the draw count, which is the disclosure; no caution is keyed on it.

**realistic:** true (bounded to low reps)  **severity:** low  **confidence:** probable
(the reps(2000) reference was not obtained)
**surface:** `bootstrap()` x `simultaneousbands`; `wooldid_esband` ado:10879

---

## INTERACTION MATRIX (p06a_matrix.do / p06b_matrix.do, reps(20), cluster(state))
Metric per cell: `e(wdid_vcvstate)`, `e(boot_method)`, and sqrt(e(V)[1,1]) against
sqrt(e(V_analytic)[1,1]).

**CORRECT** (vcvstate = bootstrap, analytic twin preserved, SE actually moved):
bare (1.229) | `att_i` (1.398) | `att_it` (1.229) | `att_itime` (1.398) | all four estimands |
`subgroups()` | `jointtests` | `adversarial` | `simultaneousbands` | `gteffects` |
`controls()` (1.290) | `interactivecontrols()` (1.125) | `anticipation()` (1.138) |
`esrelativeto(-2)` (1.247) | `boottstat` | `bootdisplay(bc)` | `bootdisplay(bca)` |
`bootdf(normal)` | `bootdf(7)` | `level(90)` | `ties` | `histogramcohorteffects` |
`summarystats` | `nofastcollapse` | `emptycellsoverride` | `cleanmatrices` |
`pretrendpower()` | `safety(off)` | `fe()`

**REFUSED, informatively, matching the documented incompatibility list:**
`impvalues` | `fastnose` | `unconditionalse` | `customvce()` | `permute()` |
`subgroups()`+`suppresssgprimaryeffects` | `oldsyntax` | `bootalpha()` under
`bootstrap(standard)` | `permonsetcohorts` without `permute()` | `ccc_absorb()` with
unit-level `i` (r459, documented)

**SILENTLY IGNORED:** none found.
**WRONG:** `timetrends` at `cluster(i)` — see BI-4.
**NOT TESTED (my own option-syntax errors, not wooldid behaviour):** `poisson` (my outcome
has negative values), `contreat()`, `oostreatedcontrols`, `coarsecohortcontrols()`,
`regtol()`, `dots()` (duplicate option in my harness), weights (aw/fw/pw), `if`,
`bootsave()`, multi-way `cluster()`.
