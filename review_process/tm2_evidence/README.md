# T-M2 mandated evidence bank (2026-08-24/25)

Defect T-M2 (review_process/mini_assessment_2026-08-23.md section 2b): in wooldid.ado's
continuous-treatment subgroup joint-pre-trend block, the two sg-continuous sites read
`r(p)[1,2]` after the `margins ... contrast(overjoint)` call (`local joint_Choriz_s<sg>`
and `local joint_Choriz_i_s<sg>`). For continuous estimands r(p) is 1x1 with the real
joint pre-trend p at [1,1]; Stata's lenient expression-context subscripting makes the
[1,2] read silently return missing, so e(joint_Choriz_s*)/e(joint_Choriz_i_s*) store "."
where significant p-values exist. (The adjacent binary twins read [1,2] correctly
because binary r(p) is 1x2.)

**Author's mandated protocol**: full wooldid console output plus the relevant matrices
(r(p) and r(b) at the two sg-continuous sites, captured MID-RUN) from at least 3
different example dataset runs, including data rich enough to rule out input degeneracy.

## Build provenance

Instrumented and stock copies were taken from the repo's `wooldidinstall/wooldid.ado`
BEFORE commit `7054878` landed (copy verified byte-identical to
`git show f4aeb60:wooldidinstall/wooldid.ado`; the ado's last-touching commit at that
vintage is `201a830`). This is deliberate: `7054878` (AMESE R1/R3/R4) adds screens that
suppress `joint_Choriz*` on cohort-clustered contreat configs, making the defect sites
unreachable on the pinning config at the new HEAD. Defect sites at this vintage:
L2688 (`joint_Choriz_s<sg>`, site1) and L2699 (`joint_Choriz_i_s<sg>`, site2).

Instrumentation (see `instrumented.diff`) is display-only: immediately before each
`local ... = r(p)[1,2]` assignment it prints a `TM2DUMP site<N> sg=<sg>` marker with
`rowsof/colsof(r(p))`, `el(r(p),1,1)`, the literal `r(p)[1,2]` read, then
`matrix list r(p)` and `matrix list r(b)`; the else-branches print an ABSENT marker.

## Identity check (instrumented build changes nothing)

Full e() dumps (%21x hex, every scalar/macro/matrix cell, plus rc and varlist, the
battery `dumpres` format) were compared per dataset:

- `cont_sg_stock.csv` vs `cont_sg_instr.csv`: **byte-identical**.
- `rich_stock.csv` vs `rich_instr.csv`: **byte-identical**.

So the TM2DUMP lines are display-only; instrumented e() results equal stock exactly.

## Runs and results

All runs Stata MP (batch `-e`), one instance at a time; rc=0 for every delivered run.
"[1,1]" below is the TM2DUMP-captured `el(r(p),1,1)`; "[1,2] read" is what the code
actually reads; "e() stored" is the post-run `e(joint_Choriz_s<sg>)` /
`e(joint_Choriz_i_s<sg>)` value (TM2ERET lines at the end of each log).

### Dataset 1: battery `cont_sg` config, exact (frozen test/data/dta_cont.dta)

`wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i)
esprelength(2) espostlength(2) subgroups(sg) jointtests` with `sg = mod(i,2)+1`.
Default estimand is att-only (wooldid.ado L4008), so only site1 executes here.

| site | sg | r(p) dims | [1,1] | [1,2] read | e() stored |
|------|----|-----------|-------------------|---|---|
| 1 | 1 | 1x1 | 6.209134328e-06 | . | e(joint_Choriz_s1) = "." |
| 1 | 2 | 1x1 | 4.108584953e-07 | . | e(joint_Choriz_s2) = "." |

The main (non-sg) continuous site, which reads [1,1], returned e(joint_Choriz) =
.0030804770378025 in the same run — the per-shape index convention the siblings follow.
Console prints "p-value for joint test that sg=1/2 pre-treatment ame horizon estimates
are 0: ." while the binary att twins print real values (0.1341 / 0.4122).

Supplementary leg (same data + `att att_it` prepended, otherwise identical): site2 also
executes; r(p) 1x1 at both sites, sg1 [1,1] = 6.209134328e-06 and sg2 [1,1] =
4.108584953e-07 at site1 AND site2; e(joint_Choriz_i_s1) = e(joint_Choriz_i_s2) = ".".

### Dataset 2: tidose-style rich panel (test/suite/preps.json `tidose_panel`, seed 447102)

50 units, onsets 5/8/10 + never-treated, time-invariant dose; `wooldid y i t ttre,
contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2)
subgroups(sg) jointtests`, sg = mod(i,3) in {0,1,2}. Att-only default: site1 executes.

| site | sg | r(p) dims | [1,1] | [1,2] read | e() stored |
|------|----|-----------|---|---|---|
| 1 | 0 | 1x1 | . | . | e(joint_Choriz_s0) = "." |
| 1 | 1 | 1x1 | . | . | e(joint_Choriz_s1) = "." |
| 1 | 2 | 1x1 | . | . | e(joint_Choriz_s2) = "." |

This leg confirms the 1x1 SHAPE on a third dataset but is itself input-degenerate: the
time-invariant unit-level dose makes the per-sg continuous pre-treatment margins come
back all-omitted (r(b) all zeros, [1,1] itself missing; main e(joint_Choriz) also ".") —
the known battery dupcol/tidose behavior. It cannot distinguish the [1,2] bug from a
genuinely missing p, which is exactly why the protocol demands Dataset 3.

### Dataset 3: purpose-built rich clean case (rules out input degeneracy)

15 cohorts (12 staggered onsets t=4..15 + 3 never-treated, one per subgroup) x 8 units
per cohort x 16 periods = 1,920 obs; dose varies across units within every cohort-period
cell (`exp(0.6*runiform() + 0.15*mod(i+unit,6))`); y has a real post dose effect
(0.6*dose) and a real sg-differentiated dose-interacted pre-trend
((0.08+0.06*sg)*dose*(t-ttre) pre-treatment), noise sd 0.35, seed 20260824.
`wooldid y i t ttre, att att_it contreat(dose) contreatcontrols(basic) cluster(i)
esprelength(2) espostlength(2) subgroups(sg) jointtests`, sg = mod(i,3)+1 in {1,2,3}.
Both sites execute (att + att_it).

| site | sg | r(p) dims | [1,1] | [1,2] read | e() stored |
|------|----|-----------|-----------------|---|---|
| 1 | 1 | 1x1 | .001191340272 | . | e(joint_Choriz_s1) = "." |
| 2 | 1 | 1x1 | .001191340272 | . | e(joint_Choriz_i_s1) = "." |
| 1 | 2 | 1x1 | .003952649944 | . | e(joint_Choriz_s2) = "." |
| 2 | 2 | 1x1 | .003952649944 | . | e(joint_Choriz_i_s2) = "." |
| 1 | 3 | 1x1 | .01215965453 | . | e(joint_Choriz_s3) = "." |
| 2 | 3 | 1x1 | .01215965453 | . | e(joint_Choriz_i_s3) = "." |

All six [1,1] p-values are genuinely significant (<0.05, two sgs <0.01) — the joint
tests have something real to find, yet e() stores "." at every site. The main sites
(reading [1,1]) store e(joint_Choriz) = e(joint_Choriz_i) = .0091446674654929 in the
same run. The STOCK run on this dataset (`tm2_rich_stock.log`, the before picture)
prints "." for all six per-sg ame/ame_it joint-test console lines while the main ame
joint test (0.0092) and the all-SGs-equal tests (0.0000) print real values; its e()
dump is byte-identical to the instrumented run's.

Design notes: the dataset must not be cohort-degenerate for these margins — a
time-invariant unit dose with i = unit id reproduces the tidose all-omitted degeneracy,
and unit-level i (one cohort per unit) under cluster(i) leaves the margins V wholly
unidentified. The delivered design follows the dta_cont / demo_vesame_pathology.do
identification pattern (i = cohort, units nested, dose varying within cohort-period
cells). Two earlier generator iterations with those pathologies were discarded.

## File inventory

- `README.md` — this file.
- `instrumented.diff` — `diff -u wooldid_stock.ado wooldid_instr.ado` (display-only dumps).
- Logs (complete Stata consoles incl. TM2DUMP mid-run matrix dumps and TM2ERET e() records):
  - `tm2_cont_sg_stock.log`, `tm2_cont_sg_instrumented.log` — Dataset 1, stock/instrumented.
  - `tm2_cont_sg_attit_instrumented.log` — Dataset 1 supplementary (att att_it; site2 coverage).
  - `tm2_tidose_instrumented.log` — Dataset 2.
  - `tm2_rich_instrumented.log`, `tm2_rich_stock.log` — Dataset 3, instrumented/stock.
- e() identity dumps (battery dumpres format, %21x): `cont_sg_stock.csv`,
  `cont_sg_instr.csv`, `cont_sg_attit_instr.csv`, `tidose_instr.csv`, `rich_instr.csv`,
  `rich_stock.csv`.
- Run scripts (self-contained, absolute paths to the tmp builds):
  `tm2_cont_sg_stock.do`, `tm2_cont_sg_instrumented.do`,
  `tm2_cont_sg_attit_instrumented.do`, `tm2_tidose_instrumented.do`,
  `tm2_rich_instrumented.do`, `tm2_rich_stock.do`.

The instrumented/stock ado copies themselves live in the (ephemeral) job tmp dir
`/Users/hegland/.claude/jobs/640b9294/tmp/tm2/`; `instrumented.diff` against the
`f4aeb60`-vintage ado blob reconstructs the instrumented build exactly.

Nothing in this directory is tracked by git (deliberately; do not add).
