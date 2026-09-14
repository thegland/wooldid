# The standard-error coverage MC that Part 8 C.6 said nobody had run — 2026-09-14

**What C.6 left open.** Thirty draws at G = 40 showed wooldid's analytic standard errors sitting
20-40% below the sampling standard deviation on two arms, while csdid and did_imputation — which
return the SAME point estimates from the same data under the same `cluster(i)` — covered
correctly in the same draws. Because the points agree, the benchmarks isolate the variance path
exactly. Thirty draws is a lead, not a finding, and C.6 wrote down the shape needed before anyone
could rule. This is that measurement at **500 draws**, run on the author's instruction of
2026-09-14 ("run it at G=40 only").

**Status of the other half.** The author had already closed the *help* half of C.6 on 2026-09-13
by removing the sentence at `wooldid.sthlp:67` that steered users toward the `i` geometry in
question ("I am skeptical of it"). This measures the statistical half, which that ruling did not
touch. The result validates the removal.

## The result

500 draws, G = 40, both DGPs, four wooldid arms, both external packages as controls in every cell.
No estimation errors on any arm. Full table in `RESULTS_G40.txt`; medians across horizons:

| arm | hom SE/SD | hom cov95 | stag SE/SD | stag cov95 |
|---|---|---|---|---|
| wooldid `esfixedbaseperiod`, **i = unit id** | 0.949 | 0.934 | **0.742** | **0.838** |
| wooldid `esfixedbaseperiod`, i = cohort id | 1.019 | 0.953 | 0.972 | 0.934 |
| wooldid pooled `esprelength(0)`, **i = unit id** | **0.799** | **0.876** | **0.661** | **0.794** |
| wooldid pooled `esprelength(0)`, i = cohort id | 1.010 | 0.950 | 1.006 | 0.946 |
| csdid `long2` *(control)* | 0.961 | 0.939 | 0.927 | 0.926 |
| did_imputation *(control)* | 0.963 | 0.939 | 0.943 | 0.936 |

**C.6's own flag rule, written before any of this ran:** an arm fails if its mean reported SE over
empirical SD is below 0.90, or its 95% coverage below 0.93, **while both controls in the same
draws sit inside those bands**.

**On `hom` the rule fires cleanly.** The pooled `esprelength(0)` arm at `i` = unit id reports
standard errors at 0.799 of the sampling SD and covers at 0.876 against a nominal 0.95, while
csdid and did_imputation cover at 0.939 in the same draws. C.6's S2 is **confirmed at 500 draws.**

**On `stag` the rule does not formally apply, and the reason is a technicality worth stating
rather than hiding.** csdid's coverage there is 0.926 — four thousandths below the 0.93 gate — so
"both controls inside band" fails and the rule abstains. But the magnitudes are not close: the two
wooldid `i` = unit arms cover at 0.838 and 0.794 against controls at 0.926 and 0.936. C.6's S1
predicted 0.838 on this exact arm and the 500-draw run returns **0.838**. Read the gate as
abstaining, not as exonerating.

**The decisive detail is in the per-horizon table.** On `stag` the empirical SD of wooldid's
estimates and csdid's are *identical to four decimals, horizon by horizon* — 0.1563, 0.1645,
0.1489, 0.2339, 0.1607, 0.1576 — because they are the same point estimates. What differs is the
reported standard error: 0.1219 against 0.1530 at the first pre horizon, 0.0845 against 0.1389 at
the fourth. Same estimates, same sampling variability, different reported uncertainty. There is no
DGP explanation available for that.

**`i` = the cohort id is clean everywhere** — 0.95 to 1.02 on the ratio and 0.934 to 0.953 on
coverage, across both DGPs and both estimators. The defect is specific to `i` = the unit id, which
is exactly what the removed help sentence used to recommend.

## What this does NOT establish

One cluster count (G = 40), one panel length, `cluster(i)` only, two DGP families. C.6 asked for
G in {20, 40, 200}; the author scoped this run to G = 40. So the result confirms the effect exists
and is large at G = 40; it says nothing about whether it shrinks with more clusters. Two horizons
on the staggered `esfixedbaseperiod` arm are much worse than their neighbours (ratios 0.567 and
0.591 against ~0.75), which is unexplained and is a thread, not a conclusion.

## Reproducing

Raw per-draw output is not tracked — it is ~6.5MB and the driver is fully seeded, so it
regenerates exactly. Four shards, 125 draws each:

    for s in 1 2 3 4; do
      <stata> -e do se_cov_mc.do <path-to>/wooldid.ado $s 125 700000 <outdir> 40 hom+stag &
    done
    python3 analyze_secov.py <outdir>
