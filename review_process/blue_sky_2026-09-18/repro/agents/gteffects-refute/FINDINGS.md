# SKEPTIC pass — verdicts on GTE-1..GTE-5, plus own findings

Build: pinned `scratchpad/bluesky/build` (sha c1d8979). All evidence from do-files in this
directory; none of the original prober's scripts was run.

---

## VERDICTS

### GTE-3 — CONFIRMED (severity medium, unchanged)
Own repro: `r02_drops.do` (blocks P1A/P1B/P1C), log `r02_drops.log`.

30 firms, 2000-2010, 12 treated in 2004 (4 of them with a true effect of 4, the other 8 of 1,
so the TRUE ATT over all 12 is 2.0), 18 never-treated.

| arm | att | N | N_dropped | treated i in e(gteffects) |
|---|---|---|---|---|
| P1A all firms observed 2000-2010 | 1.995807 | 330 | 0 | 12 |
| P1B firms 1-4 enter the register only in 2004 (== their onset) | **1.004200** | 286 | 28 | **8** |
| P1C same data, i = onset cohort instead of firm | 1.672831 | 314 | 0 | 1 |

`grep -c Caution` over the P1B block: **0**. The only trace is `N = 286 (28 obs dropped from
initial sample)`, the same line an `if`-drop or a markout produces. Mechanism as claimed
(ado:12020 `gegen mint = min(t) if ifvar==1, by(i)`; ado:12055-12057 `nopre = max(ttre<=mint)`,
`replace ifvar = 0 if nopre==1` — unit-level, no message).

CORRECTIONS to the prober's write-up:
1. The "41% swing between two legitimate choices of i" is **not evidence of a defect**. With a
   coarse i, `min(t)` is the COHORT's minimum, so a cohort with any early-entering member
   survives and the late entrants' post-onset rows are retained inside its cells (P1C: 0 dropped).
   That is a different estimand by construction, not a bug. Drop that leg of the argument.
2. The reported number is not *wrong* given identification — a treated unit with no pre-period
   cannot contribute a DiD contrast, and csdid/did_imputation drop it too. The defect is the
   **silence** and the **help wording**, nothing more. Severity medium is nevertheless right,
   on the project's own internal precedent: the arithmetically identical drop performed twelve
   lines earlier for `anticipation()` (ado:12039-12047) DOES print a named Caution.

### GTE-1 — PARTIAL (severity low, unchanged)
Own repro: `r01_low4.do` block A, log `r01_low4.log`.

Multiplier held CONSTANT within each (i,t) cell (`pop = 1 + 0.5*coh + 0.2*(t-1999)`, i = coh):

- no multiplier: per-cell `w_att_it` mass min = max = 1.00000000; per-cohort `w_att_i` 0.1667-0.3333
  (= 1/6 and 1/3, exactly one unit of mass per cohort).
- `customweightmultiplier(pop)`: per-cell `w_att_it` mass 2.50000000 - 4.00000000; per-cohort
  `w_att_i` 0.4167-1.3333.
- MECHANISM PINNED: max over post rows of |cell mass of `w_att_it` − the cell's own value of `pop`|
  = **2.2204e-15**. The cell mass is exactly the cell's multiplier value, as
  `w*v/W_it` with `W_it = total(esamp*w)` predicts.

CORRECTION: the code matches the `customweightmultiplier()` paragraph (sthlp:1107) *literally* —
"the weights normally used to produce the requested type of treatment effect (ATT, ATT_i, etc.)
will be multiplied by this variable" is precisely `w*v/W_it`. The only real gap is an
**asymmetry inside the help**: the `att` paragraph (sthlp:217) carries the qualifier
"multiplied by the customweightmultiplier() variable if one is specified", while the `att_it`,
`att_i` and `att_itime` paragraphs (sthlp:221/227/234) say "simple, unweighted mean" with no such
qualifier. That is a missing clause in three paragraphs, not a computation defect, and the
prober's own honest caveat already says so. PARTIAL.

### GTE-2 — CONFIRMED as a measurement, but it is a NEGATIVE result, not a defect (severity none)
Own repro: `r03_degen.do` blocks Q1/Q2, log `r03_degen.log`.

Under `fastnose`: `e(b_es_att)` is a real **1 x 6** matrix, `e(numpre_att) = 3`,
`e(V_es_att)` is **absent** (`e(wdid_vcvstate) = absent`; `matrix z = e(V_es_att)` yields a
1x1 missing scalar, not a matrix). Analytic control on the same data: `e(V_es_att)` is 6x6,
`e(wdid_vcvstate) = analytic-blocks`. The run prints
`Caution: no horizon of the att event study has an identified variance, so e(b_es_att) posts
without e(V_es_att)` — so the Caution sthlp:1325 promises DOES fire.

**METHOD WARNING for anyone re-running this**: `cap matrix Z = e(V_es_att)` returns `_rc == 0`
even when `e(V_es_att)` does not exist, because Stata falls back to evaluating `e(name)` as a
missing SCALAR and builds a 1x1 matrix. Test `rowsof()`, or read `ereturn list`. My first pass
got the opposite answer from exactly this trap.

CORRECTION / ADDITION: the prober did not state the sharpest version of this. On the pinned
build, sthlp:1325's list reads `(fastnose, typically timetrends, and any run whose estimation
returns an all-zero covariance)` — `impvalues` is **already off it**, `fastnose` is still on.
sthlp:284 independently says every `impvalues` margins call runs `nose` so `e(b_es_)`,
`e(V_es_)` and `e(numpre_)` are not returned. So ruling item 11 of 2026-09-13 was **half built,
and the half that was built was the correct half**: removing `impvalues` was right, and removing
`fastnose` would have made the help WRONG. The loose end is the unbuilt half of a ruling that
should not be built. Nothing in the code or the help is currently divergent, so this is a
process note with severity **none**, not a finding.

### GTE-5 — CONFIRMED (severity low, and see NEW-2 which supersedes its importance)
Own repro: `r01_low4.do` block D and `r03_degen.do` block R.

i = cohort (3 cohorts), `timetrends esprelength(2) espostlength(2) gteffects`, 13 cells:

| VCE | se missing | se exactly 0 | se < 1e-8 | min se |
|---|---|---|---|---|
| iid (default) | 0 | 0 | 0 | 4.16e-01 |
| robust | 0 | 0 | 0 | 3.72e-01 |
| cluster(unit) — below i | 0 | 0 | 0 | 1.62e-01 |
| cluster(coh) — AT i | 0 | 0 | **5** | **6.74e-15** |

So "under timetrends the coefficient-level variances are typically not identified and the se
column is missing" describes none of the four. CONFIRMED as a help-scope defect.

CORRECTION: the prober stopped one step short. In the arm where the sentence is supposed to
apply (timetrends AT the cluster level of i) the se column is not missing EITHER — it is
populated with values of order 1e-15 and p-values of order 1e-29. See NEW-2.
(At 7 clusters the same spec did go fully missing: 54/54.)

### GTE-4 — CONFIRMED (severity low, unchanged)
Own repro: `r01_low4.do` block B.

Same data, one estimand at a time, 9-row table, count of non-missing weight cells:

| requested | w_att | w_att_it | w_att_i | w_att_itime |
|---|---|---|---|---|
| att | 9 | 0 | 0 | 0 |
| att_it | 0 | 9 | 0 | 0 |
| **att_i** | 0 | **9** | 9 | 0 |
| **att_itime** | 0 | **9** | 0 | 9 |

The help's parenthetical ("missing for estimands not requested") is inaccurate for `w_att_it`.
CORRECTION worth adding to the write-up: the column is not stray — sthlp:1325 says att_it, att_i
and att_itime "share a single event-study horizon internally", and ado:1643's comment says
`w_att_it` is produced "whenever a non att weight is requested, since they are needed for the
event horizon estimates". So it is load-bearing and gating it on the option would be the wrong
fix; the one-line help edit is the right one.

---

## NEW FINDINGS (found in this pass)

### NEW-1 [high] esfixedbaseperiod silently deletes any treated unit that does not OBSERVE the exact reference relative period
Repro: `repro_esf_refperiod_drop.do` (self-contained) / `r02_drops.do` blocks P2A-P2G.

30 firms, 2000-2010, 12 treated in 2006. Firms 1-4 have a true effect of 4, firms 5-12 of 1 —
**true ATT over all 12 treated = 2.0**. Firms 1-4 have the outcome missing in 2005 (relyr = -1)
and in NO other year; they still have relyr -6..-2 and 0..4.

| arm | att | N | N_dropped | treated i in e(gteffects) |
|---|---|---|---|---|
| A complete outcome, esfixedbaseperiod | 1.997455 | 330 | 0 | 12 |
| **B one missing year at relyr -1, esfixedbaseperiod** | **1.016014** | 286 | **44** | **8** |
| C identical data, pooled reference branch | 1.994924 | 326 | 4 | 12 |
| D identical data, esfixedbaseperiod esrelativeto(-2) | 2.150968 | 326 | 4 | 12 |
| (r02 P2F/P2G) esfixedbaseperiod + timetrends, hole at -1 OR at -2 | 0.502213 | 286 | 44 | — |

`grep -c Caution` over the whole repro log: **0**.

Mechanism (ado:855-861): `refperiodobserved = max(_t == refperiod) by i`, then
`replace ifvar = 0 if refperiodobserved == 0 & ttre != .` — the whole unit, post-treatment
observations included. Under `timetrends` (ado:868-875) the unit must observe BOTH the reference
period and the one before it, so a hole at either -1 or -2 deletes it.

Why this outranks GTE-3: (a) the deleted units are fully identified — five pre-periods and five
post-periods each — and the pooled branch on the same bytes estimates them correctly; (b) the
trigger is ONE missing outcome cell in ONE specific year, the commonest imperfection in panel
data, and units that miss a particular year are exactly the units that are likely to differ;
(c) the rule is nowhere in the help — sthlp:722-727 (esfixedbaseperiod, esrelativeto) is silent
and sthlp:73 says only "cohorts treated before the sample begins are dropped"; (d)
`esfixedbaseperiod` with the default `esrelativeto(-1)` is not an exotic option — sthlp:1325
makes it a REQUIREMENT for the honestdid handoff, so it is the configuration a careful user is
steered into; (e) `e(N_dropped) = 44` rather than 4 is the only hint, and it is indistinguishable
from an `if`-drop.

Econometrically the drop is unavoidable once the base period is fixed at -1 — the defect is that
it is silent, undocumented, and moves the reported ATT by 49% of its true value here.

### NEW-2 [high] gteffects reports numerically-zero cell variances as REAL standard errors: se ~ 1e-15, p ~ 1e-45, 95% CI of width ~1e-14
Repro: `repro_gtse_degenerate.do` (self-contained) / `r03_degen.do` block R / `r04_statepanel.do` S4.

The gteffects help states, verbatim: *"degenerate zero variances are always reported as missing,
never zero"*. The guard that implements it is `wooldid.ado:13496`, inside mata `wooldid_gtx`:

    se = (v < . & v > 0) ? sqrt(v) : .

EXACT zero only, with no relative or absolute floor. When the cluster-robust VCE is rank-deficient
in the direction of a cell coefficient, the diagonal residue is ~1e-30, not 0, and it passes.

40 states, 20 years, three adoption waves, 10 never-treated; clean balanced data; i = onset
cohort; `cluster(cohid)` (i.e. cluster AT the level of i). 39 gteffects rows; 9 of them:

    i=4 relyr=-2  b=  1.09932691  se=1.035e-15  p=1.839e-45  ci=[ 1.09932691 ,  1.09932691]
    i=4 relyr=-1  b= -0.96059515  se=1.748e-15  p=1.329e-44  ci=[-0.96059515 , -0.96059515]
    i=4 relyr= 0  b=  1.65373779  se=2.685e-15  p=9.440e-45  ci=[ 1.65373779 ,  1.65373779]
    ... every cell of cohort 4, pre-treatment rows included

Controls on the same data: `cluster(state)` (below i) min se 3.21e-01, 0 tiny; `robust` min se
7.47e-01, 0 tiny. No Caution prints in any arm (`grep -c Caution` = 0), `e(es_ndegen_att)` is not
posted, and `e(wdid_vcvstate)` reads the ordinary `analytic-blocks`.

Reach (r03/r04): fires with 3 clusters (5/13 rows, min p 5.45e-30), with 7 clusters
(7/56 rows, **min p 9.5e-88**) and with 4 clusters on the 40-state panel (9/39 rows,
min p 1.84e-45), with and without `timetrends`. Does NOT fire under iid, robust, cluster below i,
a cross-cutting cluster, or with a unit-level i + cluster(unit).

Why it matters: the pre-treatment rows degenerate too, so a user reading pre-trends off
`e(gteffects)` — which the help sells as "directly comparable to the group-time tables of csdid
or xthdidregress" and as "the one inference surface that survives fast mode" — sees
p = 1e-41 on a placebo cell and concludes a catastrophic pre-trend violation. This is an
inference surface delivering the exact opposite of what the help promises for it, on a
configuration (cluster at the level of i) that the package documents, permits, and computes for.
A relative floor of the kind already used elsewhere in the package (the IM guard is
`se < 1e-10*|b|`, per the project's own record) would catch every row above: the worst
se/|b| ratio here is 1.9e-14.

#### NEW-2, reach map (`r05_map.do` / `r05_map.log`)
K onset cohorts + never-treated, 30 units each, 24 periods, `esprelength(3) espostlength(3)`:

| K | cluster AT i: tiny se / rows | min p | cluster BELOW i (unit) | 3-group cross-cutting cluster |
|---|---|---|---|---|
| 3 | 21 / 66 | 5.97e-44 | 0 | 0 |
| 6 | 18 / 123 | 2.65e-87 | 0 | 0 |
| 12 | 12 / 210 | 9.32e-167 | 0 | 0 |
| 20 | 4 / 270 | 8.36e-258 | 0 | **4 / 270, min p 1.66e-27** |

Two things this settles. First, it does not wash out with cluster count: 21 clusters still
produce degenerate cells, with a printed p of 8e-258. Second, **cluster-at-i is not the only
route** — at K=20 a 3-group CROSS-CUTTING cluster variable (a plausible `cluster(region)` with
three regions) produces the same thing. The general trigger is a cluster-robust VCE that is
rank-deficient in the direction of a cell coefficient, which is guaranteed whenever the cluster
count is small relative to the number of treated cohort-by-period coefficients.

---

## WHAT I ATTACKED AND FOUND CLEAN (negative results)
- **Boundary behaviour of both silent drops.** Pooled branch with EVERY treated unit entering at
  its own onset: `r(459)`. esfixedbaseperiod with EVERY treated unit missing at relyr -1:
  `r(459)`. Both informative, no garbage. (`r02_drops.do` P3A/P3B.)
- **The `mint` computation is NOT taken on the raw user `if`.** `ifstatement` is
  `` `ifvar' == 1 `` (ado:709) and `ifvar` already folds the markout of missing y / weights /
  controls (ado:676, 704), so a unit whose early outcomes are missing correctly counts its first
  USABLE period as `mint`. The "missing early outcome creates a phantom pre-period" route does
  not exist.
- **gteffects cell SEs are clean under iid, robust, cluster below i, and unit-level i +
  cluster(unit)** on 3-, 4-, 7-, 13- and 21-cluster designs (r03/r04/r05): 0 degenerate rows in
  every one of those 14 arms.
- **The `att_it` weighted-mean identity survives `customweightmultiplier()`** — the cell masses
  are exactly the cell's multiplier value to 2.2e-15, so `e(gteffects)` still exposes the
  reweighting truthfully.
- **The fastnose Caution promised at sthlp:1325 does fire** (verbatim in `r01_low4.log`), for
  both the att and att_it event studies.

---

## TWO NOTES THAT PRE-EMPT AN "ALREADY RULED" READING

**NEW-2 is the cell-level face of a ruled aggregate result, not a rediscovery of it.**
`ado:1767` is `cap matrix gtsed = vecdiag(e(V))` — the gteffects cell SEs are read straight off
reghdfe's stage-1 covariance, the same object that then feeds `margins` for the aggregates. The
project has RULED the aggregate consequence (analytic SEs understated to 0.66-0.75 of truth when
`cluster()` is at or above `i`; sthlp:251, se_coverage_2026-09-14). An aggregate averages many
cells and lands at 0.7x, which looks plausible. An INDIVIDUAL cell lands at 1e-15, which does not,
and which the gteffects help explicitly promises will be shown as missing. One surface was ruled;
the other contradicts a verbatim promise in its own paragraph.

**NEW-1's remedy already exists on the same code line.** `wooldid_anticesf` is invoked at
`ado:862` on this exact drop site — it attributes reference-period drops when `anticipation()`
caused them, with a named Caution. The author has therefore already decided this drop is worth
naming; the finding is that it goes unnamed when the cause is a hole in the data rather than the
`anticipation()` shift. The fix is a count and a sentence, not a change of behaviour.
Honest caveat on the size of the swing: the 2.00 -> 1.02 movement is DGP-engineered (the dropped
firms were given an effect of 4 against the others' 1). The defect is the silence; the magnitude
in any real application is whatever the dropped units' heterogeneity happens to be, and it is
zero only if missingness at the reference period is unrelated to treatment effects.
