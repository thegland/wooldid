## 6. Conclusions, one per surface

### 6.1 Per-element permutation p-values (`e(perm_p)`) - CALIBRATED

All six reported inference rows - `main_att`, `pretreatment_att` and the four event-study
horizons - sit at nominal within Monte-Carlo error at both .05 and .10 on n = 84 sharp-null
draws. The largest deviation at .05 is `main_att` at .0952 (+1.90 MC se) and the smallest is
`ES_att_pre1` at .0000 (-2.10, the conservative direction); at .10 every element is inside
1.7 MC se. The distributional checks, which use the whole p rather than one tail, are
cleaner: the chi-square on 10 equal-probability bins runs 4.8 to 19.1 against a .01 critical
value of 21.67 for every element, and every element's mean p is within 1.9 MC se of its exact
theoretical value .5050. (A printed "MC se 0.0000" wherever a share comes out exactly 0 or 1
- `ES_att_pre1` at .05 here, the tau = 3 coverages below - is the boundary artifact of
`sqrt(phat(1-phat)/n)`, not a zero-variance cell; the deviation column, which uses the se
under the nominal value, is the one to read there.) **Verdict: no evidence of miscalibration
on any element.**

### 6.2 One-sided shares (`e(perm_p_lo)`, `e(perm_p_hi)`) - CALIBRATED

Eleven of the twelve one-sided 5% cells are inside 1.5 MC se of nominal. The exception,
`ES_att_post1`'s left-tail share at .1190 (+2.90 MC se), is the largest deviation anywhere in
the sharp-null design; the same horizon's right-tail share is .0357 and its two-sided p is
.0357, so this is a one-sided fluctuation and not a shift of the horizon's null distribution
(its mean p is within 1.9 MC se and its chi-square is 5.5). With roughly thirty 5%-level
cells reported, one at |z| ~ 2.9 is what this many cells produce.

One directional pattern is visible in the table and should be named rather than left for a
reader to spot: **every** element's mean `p_lo` falls below its theoretical .5025 (.4496,
.4606, .4743, .4760, .4734, .4397; -0.8 to -2.0 MC se each) and every element's left-tail 5%
share exceeds its right-tail share - the observed statistic ranks low relative to its draws
slightly more often than it ranks high. The six elements are statistics on the same rep and
are strongly correlated within it, so this is roughly two independent observations, not six,
and noise is the plain reading. It is a pattern to re-check at larger S, not a finding at
this n.
**Verdict: calibrated; the flagged cell is listed in section 6.9 as the brief requires, but
is not read as a defect.**

### 6.3 The three alternative conventions - EXACT AGREEMENT WITH THE PACKAGE; SIZES AT NOMINAL

The package was re-run with `permabsp`, with `permrawp` and with `permside(left)` on the same
`rseed()` on four D1 reps, and the p it printed for `main_att` matched, to the last bit
(max abs diff 0.000e+00), the value the driver computed from `e(permdist_b)` and
`e(perm_comparator_b)` - 4 of 4 for each of the three. The two free identities hold on all 84
reps: the recomputed equal-tailed p equals `e(perm_p)` exactly, and the recomputed left-tail p
equals `e(perm_p_lo)` exactly. **`permside(left)` *is* `e(perm_p_lo)`** - a user who wants the
one-sided number does not need a second run.

Sizes at .05 (n = 84, MC se under nominal .0238): `permside(left)` .0714 (+0.90),
equal-tailed .0952 (+1.90), `permabsp` = `permrawp` .1071 (+2.40). At .10: .1310, .1310,
.1429 (all within 1.4 MC se). `permabsp` and `permrawp` are identical **by construction** at
199 replications and at every other conventional count (section 4), so their agreement here
is arithmetic, not evidence, and this design cannot separate them. The .1071 cell is 9
rejections where 4.2 are expected (one-sided binomial p = .024), which is unremarkable across
thirty cells but is the one worth re-checking at larger S.
**Verdict: the conventions are implemented as documented and none is miscalibrated.**

### 6.4 The fixed-metric joint pre-trend test `e(joint_horiz)` - CALIBRATED

.0476 against nominal .05 (-0.10 MC se) and .0833 against .10 (-0.51 MC se) on n = 84, with
every value on the 1/200 lattice, which also confirms the posted number is the permutation p
computed on all 199 complete draws rather than an analytic value surviving under a permutation
header. This is the cleanest surface in the report. **Verdict: calibrated.**

### 6.5 The max-|studentized| step-down `e(perm_stepdown)` - FWER CONTROLLED

| | at .05 | at .10 |
|---|---|---|
| uncorrected: any horizon's `p_raw` <= alpha | .1190 | .2381 |
| FWER: any horizon's `p_stepdown` <= alpha | .0952 (+1.90 MC se) | .1310 (+0.95) |
| FWER: any horizon's `p_singlestep` <= alpha | .0952 | .1310 |
| per-horizon `p_raw`, post0 / post1 | .0714 / .0595 | .1190 / .1429 |

Two facts. First, **the correction removes the multiplicity**: with two post horizons the
uncorrected any-rejection rate is .119 at .05 and .238 at .10 - about twice nominal, as two
roughly independent 5% tests should be - and the step-down brings both back inside 2 MC se of
nominal. Second, `min p_stepdown` and `min p_singlestep` reject on exactly the same reps,
which is not a bug: the step-down's first step *is* the single-step max over all horizons, so
with H = 2 they can differ in which horizon is credited but never in whether anything is
rejected. The per-horizon `p_raw` values were also checked against the package's own
`permabsp` p for the same horizon on the verification reps and matched exactly (4 of 4 per
horizon), confirming the studentization by one constant is rank-preserving as intended.
**Verdict: FWER is controlled and the correction is doing the work it claims.**

A note on reading this cell: at an interim n = 60 it stood at .1167 (+2.37 MC se) and settled
to .0952 (+1.90) by n = 84. That is the size of the noise here, and the reason section 7
insists on the n-limits.

### 6.6 The shift interval `e(perm_ci)` under the sharp null - EXACTLY THE DUAL OF THE TEST

| column | coverage of 0 | 1 - (the .05 size from 6.1) |
|---|---|---|
| `main_att` | .9048 | 1 - .0952 |
| `ES_att_post0` | .9286 | 1 - .0714 |
| `ES_att_post1` | .9643 | 1 - .0357 |

Coverage equals one minus the equal-tailed size on every column - and the identity was checked
rep by rep, not just share by share: `ll <= 0 <= ul` agreed with `p > .05` on **84 of 84 reps
for all three columns**. That is what the construction promises - the interval is the estimate minus the equal-tailed quantiles of the
null draws, so it contains 0 exactly when the two-sided p exceeds .05 - and it is a useful
internal-consistency check that the interval and the p-value are computed off the same
quantiles. **No interval was blank on any of the 84 reps** (`e(perm_kernel)` row `offcentre`
never fired): the null band never excluded 0 in this design. **Verdict: correct by
construction and correct in fact.**

### 6.7 The shift interval under a planted effect - CALIBRATED AT tau = 1, OVER-COVERING AT tau = 3, AND THE WIDTH GROWS EXACTLY AS F035 SAYS

| tau | `main_att` | `ES_att_post0` | `ES_att_post1` | n |
|---|---|---|---|---|
| 1 (D2) coverage | .9250 (-0.73) | .9750 (+0.73) | .9000 (-1.45) | 40 |
| 3 (D3) coverage | 1.0000 (+0.92) | 1.0000 (+0.92) | 1.0000 (+0.92) | 16 |

and the paired per-rep width ratios against the same rep's tau = 0 interval:

| ratio | `main_att` | `ES_att_post0` | `ES_att_post1` |
|---|---|---|---|
| tau = 1 / tau = 0 (mean) | 1.209 | 1.079 | 1.075 |
| tau = 3 / tau = 0 (mean) | 2.380 | 1.620 | 1.744 |

At tau = 1 coverage is at nominal on all three columns; no interval was blank. At tau = 3 the
interval covers on **every** rep and is 1.6 to 2.4 times as wide as the same rep's tau = 0
interval - the documented behaviour, and the mechanism is visible directly in
`e(perm_sd)`: the null-distribution SD for `main_att` goes .425 (tau = 0) -> .495 (tau = 1)
-> .904 (tau = 3). The reassigned onsets move real post-treatment outcomes into the draws, the
null distribution widens with the true effect, and the shift interval inherits that width.
**Verdict: the help text's claim that the interval "is calibrated only near the null" is
confirmed quantitatively - it is calibrated at tau = 1 (an effect of about 1 residual SD) and
conservative by a factor of ~2.4 in width at tau = 3.** No cell here is a defect, and the tau = 3 over-coverage is
**not** in section 6.9's list: at n = 16 the deviation is only +0.92 MC se. The paired width
ratio, a mean of a per-rep ratio, is the informative statistic for that arm - not the
coverage share.

### 6.8 `permtstat` - THE t-SCALE DOES NOT NORMALISE UNDER `cluster(i)`, AND THIS RUN CANNOT TEST THE WIDTH CLAIM

At tau = 1, n = 24: coverage .9167 / .8750 / .8750 (-0.75, -1.69, -1.69 MC se) on
`main_att` / post0 / post1, no blank intervals, and the analytic SE was available on every rep
for every column (0 of 24 missing or non-positive). The `permtstat` interval is **wider** than
the default interval on the same draws - paired mean ratios 1.101, 1.196, 1.111 - and the SD
of the null **t**-ratios is 1.43, 1.67, 1.55, not ~1. That is the `cluster(i)` caution the
package already prints (`permtstat` with `cluster()` at the `i` level: with one observation
per unit-period the treated cells are absorbed by their own cohort-period terms and the
cluster-robust SE cannot see their variance): studentizing by an SE that does not track the
draw's dispersion inflates the spread of the null t rather than standardising it. This is
consistent with that caution **and** with the generic downward bias of a cluster-robust SE at
G = 12 clusters, which on its own puts the SD of a t-ratio above 1; this run does not separate
the two mechanisms.

**A limitation to record, not a finding:** this run has no `permtstat` arm at tau = 0, so it
**cannot** test the help's claim that `permtstat` "removes that growth for the post-treatment
rows". The only clean test is `permtstat` width at tau = 1 against `permtstat` width at
tau = 0 on the same draws, and D1 is the default arm. Adding a fifth phase (`permtstat`,
tau = 0) is the single highest-value addition to the next run of this driver, and it costs one
extra run per rep.

### 6.9 Surfaces off nominal by more than 2 MC standard errors

Four cells, all from the sharp-null design's 5% tail; none from any coverage table.

| surface | value | nominal | deviation | n |
|---|---|---|---|---|
| `ES_att_pre1` equal-tailed p @ .05 | .0000 | .05 | -2.10 MC se (conservative) | 84 |
| `ES_att_post1` `e(perm_p_lo)` @ .05 | .1190 | .05 | +2.90 MC se | 84 |
| `main_att` under `permabsp` @ .05 | .1071 | .05 | +2.40 MC se | 84 |
| `main_att` under `permrawp` @ .05 | .1071 | .05 | +2.40 MC se (the same rejection set as `permabsp`) | 84 |

Around thirty 5%-level cells are reported, so under perfect calibration one or two at
|z| ~ 2.4 are expected; the last two rows are one event counted twice. None of the four
survives the global checks on the same element (mean p and the 10-bin chi-square are both
inside their bands everywhere), and the flagged elements are not the same across statistics.
**This is a list to re-check at larger S, not a list of defects.** If any one of them is to
be chased, it is `main_att` under `permabsp` - the one cell where the two-sided default
(.0952) and the absolute-value convention (.1071) both lean high on the same coefficient.

## 7. What this run can and cannot detect

The nominal levels are attained exactly on the 199-replication lattice (section 4), so there
is no lattice slack and every table above is a straight test of a share against a known
constant. What limits the run is `n`, not the design. At n = 84 the two-MC-se band around a
nominal 5% is [0.3%, 9.7%] and around a nominal 10% is [3.5%, 16.5%]; at n = 40 the band
around 95% coverage is [88.2%, 100%]; the tau = 3 arm (n = 16) and the `permtstat` arm
(n = 24) are pilots whose coverage shares are barely informative - their **paired width
ratios**, which are means of a per-rep ratio rather than tail shares, are what those two arms
actually establish. The run therefore detects **gross** miscalibration - a surface at twice or
half its nominal size - and not a one-percentage-point tilt. The two global checks (the
chi-square on 10 equal-probability bins, and the mean of the p against its exact theoretical
value .5050 equal-tailed / .5025 one-sided) use the whole p-value distribution rather than one
tail and carry most of the calibration weight at this `n`; the 5%/10% tail shares are the
headline but the noisiest cells in the report.

For a follow-up run, in priority order: (a) raise S on the sharp-null design to 300+ so the
5% cells have MC se ~1.3pp; (b) add a `permtstat` tau = 0 phase so section 6.8's width claim
can be tested; (c) raise the tau = 3 arm above n = 16.

## 8. Files

| file | what |
|---|---|
| `test/suite/mc/mc_ri2.do` | the driver (four phases, wall-budgeted, one CSV row per rep and arm) |
| `test/suite/mc/run_mcri2.sh` | shard launcher with the instance-cap wait |
| `test/suite/mc/analyze_mcri2.py` | the analyzer that produced every table above |
| `test/suite/mc/out2/permutation/mc_ri2_s1.csv`, `mc_ri2_s2.csv` | the raw per-rep records |
| `test/suite/mc/out2/permutation/lattice_exact.txt` | the exact lattice enumeration of section 4 |
| `test/suite/mc/out2/permutation/smoke/` | the plumbing smoke pass (`permreps(39)`, one rep per phase, all four identities checked) |
| `test/suite/mc/out2/permutation/F058_perm_prose.md`, `F058_perm_conclusions.md`, `assemble_f058_perm.sh` | the pieces this record was assembled from |

Every CSV row carries its `dgpseed` and `rseed`, so any rep is reproducible on its own:
rebuild the R1 panel at that `dgpseed` and re-run the call in section 2 with that `rseed()`.
Nothing under `wooldidinstall/` was modified. One read-only `git log --oneline -3` was run
for orientation in the first minute of the session; no other `git` command was issued.
