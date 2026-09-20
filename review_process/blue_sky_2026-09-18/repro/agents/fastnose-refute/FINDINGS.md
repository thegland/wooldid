# fastnose-refute — FINDINGS

Pinned build asserted in every log (`which wooldid` -> scratchpad/bluesky/build, v1.4).
Every repro below generates its own panel and can be run alone.

---

## VERDICT: FASTNOSE-03 — **CONFIRMED**, severity **high** (as filed; see the caveat below)

**Own repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fastnose-refute/r1_fn03.do`
(escalation: `r3_escalate.do`, controls: `r4_controls.do`, mechanism: `r6_omarker.do`)

My panel, not theirs: 105 units, T = 2010..2019, onsets 2013/2016/2019 + never-treated,
dose varying WITHIN cohort, subgroup = a late-adopter region (`region==2` <=> the 2019
cohort), `esprelength(3) espostlength(2)`.

### Measured (r1_fn03.log)
`contreatcontrols(basic)`:
```
BASIC_ANA_bES_ame_s1_COLS: post0 post1                 numpre_ame_s1 = 0
BASIC_FST_bES_ame_s1_COLS: pre3 pre2 pre1 post0 post1  numpre_ame_s1 = 3
                            values:  0  0  0  -.2629225468  .2017838828
BASIC_ANA_bES_ame_s2_COLS: c1   (1x1, value MISSING)   numpre_ame_s2 = .
BASIC_FST_bES_ame_s2_COLS: pre3 pre2 pre1 post0        numpre_ame_s2 = 3
                            values:  0  0  0  0
```
`contreatcontrols(untreatedZero)`: same shape (analytic s1 `post0 post1`, s2 `post0`;
fastnose s1 `pre3..post1`, s2 `pre3..post0`, all extra entries exactly 0).

Control, same spec **without** `subgroups()`: both modes give `pre3 pre2 pre1 post0 post1`
and `numpre_ame` = 3 — identical. The divergence is subgroup-specific, as claimed.

`e(es_ndegen_ame_s#)` is ABSENT in both modes, so nothing flags the zeros. In the analytic
run four Cautions name the drops verbatim ("the continuous-treatment relative period
effect(s) pre3 pre2 pre1 post0 are not estimable for subgroup 2 (margins reports them as
omitted)", "no continuous-treatment relative period effect is estimable for subgroup 2, so
its continuous-treatment event study is skipped entirely"). **Under `fastnose` not one of
those Cautions fires.** The printed subgroup-2 table then reads `pretreatment | ame | 0`
and `ES_ame | pre3 pre2 pre1 post0 | 0 0 0 0`.

### Mechanism — PROVEN, not inferred (r6_omarker.log, no wooldid involved)
`wooldid_sglbl` (ado:12575) drops a horizon iff its margins colname carries the `o.`
prefix. On a period-specific dose-slope design with three terms omitted for collinearity:
```
margins , dydx(x) over(t) df(500)  ->  r(table) colnames: 1bno.t 2o.t 3o.t 4.t 5.t 6.t   (3 markers)
margins , dydx(x) over(t) nose     ->  r(table) colnames: 1bn.t  2.t  3.t  4.t 5.t 6.t   (0 markers)
```
`nose` strips the marker. ado:1981 forces `sgmarginsdof` to `df(df_r)` for the BINARY
subgroup-by-horizon call (BS-26); the CONTINUOUS twin at ado:3290 still threads
`marginsdof`, which ado:1973 has set to `nose`. Confirmed at the wooldid level by r4:
BINARY + subgroups + ES + bootstrap runs clean (B3 rc 0), the identical contreat spec dies
(B4 rc 459).

### Severity: high stands, and the root cause reaches further (see SKEP-01/SKEP-02)
`high` is right for the claim AS FILED: the user typed `fastnose`, the help tells them to
re-run without it, and no inference is attached to the fabricated rows. I considered
raising it to critical and decided against it for one reason I checked rather than assumed
— on the POOLED (no-`subgroups()`) path the analytic run ALSO reports
`pre3 pre2 pre1 = 0` (r1: `NOSG_ANA_bES_ame` = 0 0 0 -.1752816978 .2017838828), flagged by
a Caution and `e(es_ndegen)`. Publishing a structural zero as an ame horizon is documented
design on the pooled path (sthlp:1325, sthlp:903). What is NOT design is the subgroup
row-set divergence, and the total suppression of every drop Caution and of
`e(es_ndegen_ame_s#)` alongside it. The same root cause does escalate on two surfaces the
prober never reached — `permute()`'s default and `bootstrap()` — filed separately below so
the severity is not double-counted.

### Corrections to the prober's write-up
- Their extra row was a post-treatment horizon; on my geometry the extra rows are the
  pre-treatment ones (plus, for subgroup 2, a `post0` of exactly 0 where the analytic run
  reports no ame event study at all). Both are the same defect; "a POST-treatment horizon"
  is not essential to it.
- sthlp:903 is a partial hedge they did not cite: the dose-slope structural-zero screen
  "is skipped under {it:fastnose}". It covers the *aggregate* ame screen, not the
  subgroup-by-horizon row set, so sthlp:292's "fast mode therefore reports the same
  event-study rows as the full run" is still contradicted. Worth quoting both.

---

## VERDICT: FASTNOSE-01 — **CONFIRMED**, severity **medium** (unchanged)

**Own repro:** `.../fastnose-refute/r2_fn01_fn02.do`  (Part 1)

My own binary panel (96 units, T=2001..2012, 4 treated cohorts + never-treated),
`set seed 4242` before each of three arms, `permreps(40)`, `cluster(cohort)`.

```
ARM A  permute + permanalyticse           rc 0  vcvstate=analytic-blocks  e(V) rc 0  test rc 0   perm_analyticse=1
ARM B  permute + permanalyticse + fastnose rc 0  vcvstate=absent          e(V) rc 111 test rc 321 lincom rc 321  perm_analyticse=1
ARM C  permute (fast default)              rc 0  vcvstate=absent          e(V) rc 111             perm_analyticse=.
```
`e(wdid_vcespec)` in ARM B is the permanalyticse string verbatim ("se_analytic = the
analytic cluster(cohort) standard error, shown for reference only"), and the footer prints
"The standard-error column is the ANALYTIC standard error (se_analytic, permanalyticse)
... e(V) and e(V_es_*) hold the analytic covariance" while `e(V)` does not exist. The
printed `se_analytic` column is `.` on every row. ARM B vs ARM C: `e(perm_p)` in %21x and
`e(wdidmainresults)` are identical line for line — `permanalyticse` is a complete no-op.

### Two corrections, both in the prober's favour
1. They cite only code. **sthlp:1234** (the engine-compatibility table) states: "{cmd:fastnose}
   refused / it {ul:is} the default. Fast mode is what makes the permutation affordable;
   {cmd:permanalyticse} opts back in." So the help says `permanalyticse` opts back in and
   the code lets a typed `fastnose` override it — a documented-vs-actual divergence, not
   only an internal inconsistency.
2. The run does **not** print the fast-mode header. Line 318 of the log is
   "Wooldid Estimation for Outcome: y; Standard Errors: randomization inference (...)",
   and the "SEs and related materials are not estimated in fast mode." line is absent.
   Nothing in the output says fast mode is in force.

Severity stays medium: every number is a correct permutation quantity and the missing
`e(V)` fails loudly (r321), so the damage is the false footer, the false
`e(wdid_vcespec)` and the machine-readable `e(perm_analyticse)=1`.

---

## VERDICT: FASTNOSE-02 — **REFUTED**, severity **none**

**Own repro:** `.../fastnose-refute/r2_fn01_fn02.do`  (Part 2)

The behaviour reproduces exactly: `e(gteffects)` is 32x15 in both modes and the two
matrices are byte-identical in `%21x`, printed under
"Standard Errors: not computed (fast mode) / SEs and related materials are not estimated
in fast mode."

But the claim's charge — "sthlp:286 says ... and the gteffects paragraph never mentions
fastnose", i.e. a documentation gap — is **false**. sthlp:1149-1150, inside the
`gteffects` paragraph:

> "Not available with {cmd:contreat()} ... **Fully available under {it:fastnose} (the
> regression's own VCE is intact even though the {it:margins} stage skips variances),
> making {it:gteffects} the one inference surface that survives fast mode.**"

That is the exact sentence the finding says is missing, it names the exact mechanism the
finding "discovered" (the stage-1 VCE), and the fastnose paragraph's "SEs and related
materials are not estimated in fast mode." line is itself documented at sthlp:286 as the
fast-mode marker. There is no divergence to close.

---

# NEW FINDINGS (mine; same root cause as FN-03 but surfaces the prober never reached)

## SKEP-01 [high] `permute()`'s fast-mode DEFAULT publishes a whole subgroup ame event study the analytic run skips entirely, as zeros carrying permutation p-values and shift intervals
**Repro:** `.../fastnose-refute/r3_escalate.do` (ARM P), controls in `r4_controls.do` (ARM V)

`wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) fe(unit) cluster(cohort)
subgroups(region) esprelength(3) espostlength(2) permute(shuffle) permreps(20)` — no
`fastnose` typed — returns rc 0 and prints, for subgroup 2:
```
                  estimate  se_analytic  t_perm      p_perm   lb_95ci  ub_95ci  relyr
      main:ame           0            .       .          .4  -.1923960217  .0001773005   .
pretreatment:ame         0            .       .           1             0            0   .
   ES_ame:pre3           0            .       .           1             0            0  -3
   ES_ame:pre2           0            .       .           1             0            0  -2
   ES_ame:pre1           0            .       .           1             0            0  -1
  ES_ame:post0           0            .       .           1  -.1923960217            0   0
```
The analytic run on the identical data (r4 ARM V) reports **no `ES_ame` row at all** for
subgroup 2 and says why ("no continuous-treatment relative period effect is estimable for
subgroup 2, so its continuous-treatment event study is skipped entirely"). Subgroup 1 is
hit too: `ES_ame:pre3/pre2/pre1 = 0` with `p_perm` 1, 1, .95 where the analytic run reports
only `post0` and `post1`.

Why high: (a) `permute()` takes the fast path BY DEFAULT (ado:5448), so the user never
types `fastnose` and is never told to re-run; (b) the zeros carry inference, so a subgroup
pre-trend check reads as a flawless null (three pre-horizons at exactly 0, p = 1) that is a
pure artifact of the regression's omitted terms; (c) every Caution the analytic run prints
about the drops is suppressed, and `e(es_ndegen_ame_s#)` is absent, so nothing anywhere
flags them.

NOT the complaint: on the POOLED path the analytic run reports `pre3 pre2 pre1 = 0` as
well, flagged by `e(es_ndegen)` — that is documented design (sthlp:1325). The defect is the
SUBGROUP row set, where the analytic run drops the rows and says so and the permute()
default publishes them with no flag of any kind.

Sharpening: subgroup 1's `ES_ame:pre1` comes back `p_perm` = .952 with shift interval
[-.1859513891, 0]. The observed statistic there is a structural zero while the
replications mix structural zeros with genuine estimates from reassignments that made the
horizon estimable — so that p is not a permutation p-value of anything. "Inference
attached" is worse than "value 0".

Root cause: shared with FASTNOSE-03 (ado:3290 threads `marginsdof` = `nose`; ado:2912's
binary twin was given `df(df_r)` by BS-26 and this one was not).

## SKEP-02 [high] `bootstrap()` and `permute() permanalyticse` abort r(459) on contreat() + subgroups() + an event-study window
**Repro:** `.../fastnose-refute/r4_controls.do` (B1-B5), also `r3_escalate.do` (ARM PA/BT)

```
Error: the rep-0 coefficient vector does not match the observed run's column set;
bootstrap accumulation would be misaligned. This is unexpected - please report it.
r(459)
```
(the `permute() permanalyticse` arm prints the same with "permutation accumulation".)

Control matrix on one panel, `bootstrap(standard) reps(20)`:
```
B1 contreat + ES, NO subgroups          rc 0
B2 contreat + subgroups, NO ES window   rc 0
B3 BINARY  + subgroups + ES             rc 0      <- BS-26's fix works on the binary call
B4 contreat + subgroups + ES            rc 459
B5 B4 with robust instead of cluster()  rc 459
```
So the abort is exactly the three-way interaction, and not cluster-dependent. The
observed run is analytic and drops the non-estimable horizons; the rep-0 comparator runs
`fastnose` (ado:322 for bootstrap, ado:126 for permute) and keeps them, so the column sets
cannot align. This is the failure the BS-26 comment at ado:1976-1978 says it fixed — fixed
for the binary call at ado:2912, never for the continuous twin at ado:3290. Continuous
treatment + heterogeneity by subgroup + an event study + a bootstrap is a mainstream
applied specification and it cannot be run at all.

This is not an edge case in the data either. The structural zeros exist because `dose` is
0 in untreated periods — which is exactly what sthlp:897 instructs ("Must be defined for
all observations (set to 0 if not naturally defined within untreated cohort-periods)"). So
the abort fires on the help's own recommended input shape. The "This is unexpected -
please report it" wording confirms it is an unmet internal assertion, not a design refusal.

Root cause: shared with FASTNOSE-03.

## SKEP-03 [low] Where a subgroup's continuous event study is "skipped entirely", `e(b_es_ame_s#)` is still posted as a 1x1 matrix of MISSING under the junk colname `c1`, with `e(numpre_ame_s#)` = .
**Repro:** `.../fastnose-refute/r1_fn03.do` (ARM basic / ANALYTIC, subgroup 2)

```
BASIC_ANA_bES_ame_s2_COLS: c1
  __000002[1,1]   c1
  r1               .
BASIC_ANA_numpre_ame_s2 = .
```
This is the ANALYTIC run, right after "no continuous-treatment relative period effect is
estimable for subgroup 2, so its continuous-treatment event study is skipped entirely."
`e(b_es_ame_s#)` is a documented honestdid/pretrends handoff object; here it exists but
holds a single missing value under a column name that is not a horizon label, and
`e(numpre_*)` is missing rather than 0. Withholding the matrix (as the dose-slope screen
does elsewhere, sthlp:1325) would be the consistent behaviour. Low: the value is missing,
so a consumer is unlikely to be silently misled.
