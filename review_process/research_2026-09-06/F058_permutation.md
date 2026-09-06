# F058 - extended calibration simulations: PERMUTATION

**Status: research record. READ-ONLY on the package - no file under `wooldidinstall/` was
touched. One read-only `git log --oneline -3` was run for orientation in the first minute;
no other `git` command was issued.** Author ruling 2026-09-06 on F058: run the extended
calibration simulations. This is the permutation part (`mc_ri2.do`); the bootstrap part is a
separate record.

Package tree: branch `bugfixes`, HEAD `01a7f7a` (the final tree for these surfaces, after
F064 `permside()`/one-sided shares and F065 the max-|studentized| step-down).

*Provenance note.* `wooldidinstall/wooldid.ado` was written to by another agent at 17:13,
after both of my shards had finished (17:11:54). This does not contaminate the run: each
shard `do`s the `.ado` **once**, at shard start, so an edit landing mid-run cannot reach a
running shard. Shard 1 loaded the file at 16:09:53 and shard 2 at 16:10:18, 25 s apart, so
every one of the 176 rows below comes from one consistent snapshot of the code. Nothing in
this work modified any file under `wooldidinstall/`.

## 1. What was measured

Everything the `permute()` header can put on screen, on two designs.

| surface | e() name | design |
|---|---|---|
| per-element permutation p (pooled att, pretreatment att, each ES horizon) | `e(perm_p)` | sharp null |
| one-sided shares | `e(perm_p_lo)`, `e(perm_p_hi)` | sharp null |
| the three alternative p conventions | `permabsp`, `permrawp`, `permside(left)` | sharp null |
| fixed-metric joint pre-trend test | `e(joint_horiz)` (`jointtests`) | sharp null |
| max-\|studentized\| step-down over the post horizons | `e(perm_stepdown)` rows `p_raw`, `p_singlestep`, `p_stepdown` | sharp null |
| shift interval | `e(perm_ci)` rows `ll`, `ul` | planted effect |
| shift interval on the t-ratio scale | `e(perm_ci)` under `permtstat` | planted effect |

## 2. Designs

**Shared DGP** - `mc_ri.do`'s R1, which is the configuration the shipped RI size battery
already uses:

```
12 clusters i, 10 periods t = 1..10
onset labels {5 x3, 7 x3, . x6} assigned by a uniform shuffle
y_it = a_i + 0.1 t + e_it,  a_i ~ N(0,1) iid,  e_it ~ N(0,1) iid
call: wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) jointtests
                          permute(shuffle) permreps(199) rseed(<per-rep>)
```

The reported window therefore carries six inference rows - `main_att`,
`pretreatment_att`, `ES_att_pre2`, `ES_att_pre1`, `ES_att_post0`, `ES_att_post1` - so the
per-horizon sizes, the joint pre-trend test (two pre horizons, the minimum the analytic
test needs) and the step-down (two post horizons) all exist on one run.

**Design 1 - sharp null (phase D1).** `tau = 0`: the outcome does not depend on the
assignment at all, so the sharp null holds and every permutation p-value is being asked for
its size.

**Design 2 - planted constant effect (phases D2, D3, D4).** `y := y + tau*1{t >= ttre}`
with `tau = 1` (D2, default arm; D4, `permtstat` arm) and `tau = 3` (D3, default arm).
The estimand is `pooled_reference_period`, so with a constant effect the truth is `tau` on
`main_att` and on every post horizon, and 0 on the pre-treatment rows.

**Pairing.** The DGP seed and the `rseed()` are functions of `(shard, rep)` only, so rep
`r` of D1, D2, D3 and D4 shares one draw of `(a_i, e_it, onset labels)` and one set of 199
permutation assignments. Width growth is therefore a **paired per-rep ratio**, not a
difference of two independently noisy means, and D1 supplies the `tau = 0` width for free.

## 3. How the alternative p conventions were obtained

`permabsp`, `permrawp` and `permside()` are alternative **readings of one accumulator**:
`wooldid_permtail` (the `absp` convention code, 0/1/2/3/4) computes all of them from the
same replication vector, and none of them changes the draws. The driver therefore

1. runs the default arm once per rep and stores `e(permdist_b)` (the replication
   coefficients) and `e(perm_comparator_b)` (the observed statistic that is ranked);
2. recomputes, in Mata, `wooldid_permtail`'s four conventions on the `main_att` column -
   equal-tailed `min(1, 2 min(plo, phi))`, `permabsp` `(1+ma)/(n+1)`, `permrawp` `ma/n`,
   `permside(left)` `plo`, with `plo = (1+m0)/(n+1)`, `phi = (1+m1)/(n+1)`,
   `ma = #{|s_r| >= |s_obs|}`;
3. **re-runs the package** with `permabsp`, with `permrawp` and with `permside(left)` on
   the *same* `rseed()` on every 15th rep, and asserts exact equality against the values
   computed in step 2.

Step 3 is the licence for step 2: three genuine package runs cost the same as three reps of
the headline size measurement, so the identity is checked on a subset and the sizes are
then reported on the full sample. Section 5.2 reports the agreement. Two further identities
are checked for free on every rep: the driver's recomputed **equal-tailed** p must equal
`e(perm_p)`, and its recomputed left-tail p must equal `e(perm_p_lo)`.

## 4. What "calibrated" means here, exactly

At `permreps(199)` the permutation p-value lives on a 200-point lattice and **the nominal
levels .05 and .10 are attained exactly**, so there is no lattice slack to hide behind.
Write `R` for the rank of the observed statistic among the 200 values (199 draws plus the
observed); under exchangeability `R` is uniform on 1..200.

| convention | p | rejects at .05 iff | size |
|---|---|---|---|
| equal-tailed (default) | `2 min(R, 201-R)/200` | `R <= 5` or `R >= 196` | 10/200 = .0500 |
| `permabsp` | `(1+ma)/200` | `ma <= 9` | 10/200 = .0500 |
| `permrawp` | `ma/199` | `ma <= 9` (since 10/199 = .0503 > .05) | 10/200 = .0500 |
| `permside(left)` | `R/200` | `R <= 10` | 10/200 = .0500 |

and at .10 all four give 20/200 = .1000. (Enumerated, not simulated:
`out2/permutation/lattice_exact.txt`.)

Note the by-product, which bears on F052: **`permrawp` and `permabsp` have the same
rejection set here, and generally.** `permrawp` rejects iff `ma <= alpha*n`, `permabsp` iff
`ma <= alpha*(n+1) - 1 = alpha*n + alpha - 1`; the two coincide whenever no integer lies in
`(alpha*n + alpha - 1, alpha*n]`, which holds at every conventional replication count for
which `alpha*(n+1)` is an integer - 99, 199, 999, 1999 at both `alpha = .05` and `.10`. The
two options therefore differ in the **printed number** (`permrawp` can print an exact 0,
`permabsp` cannot go below `1/(n+1)`) and in the Monte-Carlo se / Clopper-Pearson interval
attached to it - not in what a 5% or 10% rule decides. This simulation consequently cannot
separate their sizes and is not evidence that it could.

The same argument applies to `e(joint_horiz)`: its metric `Sigma` is built from the 199
draws *plus* the comparator, so it is a symmetric function of all 200 and `W_obs` is
exchangeable with the `W_r`; `p = (1 + #{W_r >= W_obs})/200` then rejects at .05 iff
`W_obs` is among the top 10 of 200. Nominal exactly.

Consequently the yardstick used below is the **MC standard error under the nominal value**,
`sqrt(alpha(1-alpha)/n)`, and a cell is flagged when the observed share is more than 2 of
those away from nominal. (The tables also print `sqrt(phat(1-phat)/n)`, the conventional
reported MC se.)

## 4b. Run parameters - what S was chosen, and why not 300

The brief's target was S = 300 per design. That was not attainable inside the ~90-minute
runtime cap on this machine and S was cut accordingly, as the brief allows.

| item | value |
|---|---|
| shards | 2 (the cap on concurrent shards of my own), seeds `900000 + shard*100000 + rep` |
| replications per run | `permreps(199)` throughout |
| per-shard wall budget, by phase | D1 1900 s, D2 800 s, D3 300 s, D4 600 s (3600 s = 60 min per shard) |
| observed seconds per `wooldid` call | ~35-48 s default arm, ~70-95 s under `permtstat` (machine shared with 3-5 other Stata jobs) |
| verification cadence | every 15th D1 rep also runs the three package convention arms |
| launch -> both shards done | 16:09:53 -> 17:11:54 = 62 min wall (0 min lost to the cap wait) |
| realised n | D1 sharp null 84 reps (42 per shard, verified from the CSVs); D1 verification arms 4 reps x 3 package runs; D2 tau=1 default 40; D3 tau=3 default 16; D4 tau=1 permtstat 24 |
| realised phase elapsed | shard 1: D1 1906 s, D2 818 s, D3 330 s, D4 644 s; shard 2: 1913 / 812 / 331 / 640 s |

At ~40 s per default call, S = 300 for the sharp-null design **alone** is 300 x 1.2 x 40 s
= 4 hours of Stata on 2 shards, i.e. ~2 hours of wall time before the planted-effect
designs start. The budgets above were set to spend the available hour in priority order:
the sharp-null sizes (the headline) first, then the tau = 1 coverage, then the tau = 3
width growth, then the `permtstat` arm as a pilot. Each phase carries a wall-clock guard
and each rep is appended to the CSV with its own open/write/close, so a shard stopped at
its budget still leaves every completed rep on disk. The realised n per cell is printed in
every table below and must be read with the tables.

**The instance-cap counter.** The brief's literal probe

```
pgrep -fl stata-mp | grep -c " -e do "
```

also matches other agents' *shell pipelines*: their own monitoring commands contain both
`stata-mp` and ` -e do `, so the probe counts them as Stata jobs. At 16:06 it returned 12
while `ps` showed four real Stata processes (PIDs 2865, 4559, 4782, 98497 - three other
agents' shards plus one smoke run); the literal rule would have blocked every launch
indefinitely. The launcher therefore counts processes whose ARGV[0] is the Stata binary:

```
ps -eo args= | grep -c "^/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do"
```

which is what the cap means. Both shards were launched only with that count below 6 (3 and
4 respectively), and never more than 2 of my own ran at once. Flagging this because it is a
deliberate reinterpretation of an explicit instruction, not a silent deviation.

**One log file, two shards.** Stata's `-e do` names its log from the do-file arguments, and
both shards were given the same arguments apart from the shard number, so both wrote to
`permutation 1900 800 300 600 900000 200 15 199.log` and their lines interleave there. The
CSVs are the authoritative record; the log is usable only for the eight
`MCRI2_PHASE_DONE` lines quoted above, which all survived.

## 5. Results


rows read: 176  from 2 shard files

- phase D1 arm absp: 4 rows, 4 with rc=0
- phase D1 arm left: 4 rows, 4 with rc=0
- phase D1 arm main: 84 rows, 84 with rc=0
- phase D1 arm rawp: 4 rows, 4 with rc=0
- phase D2 arm main: 40 rows, 40 with rc=0
- phase D3 arm main: 16 rows, 16 with rc=0
- phase D4 arm tstat: 24 rows, 24 with rc=0

### 5.1 Integrity of the run

- D1 reps with rc=0: 84
- column-name lookups that fell back to position: 0
- reps with e(perm_completed) != 199: 0; with e(perm_failed)>0: 0
- usable draws per rep for main_att: min 199 max 199
- recomputed equal-tailed p vs e(perm_p) main_att: 0 mismatches of 84
- recomputed left-tail p vs e(perm_p_lo) main_att: 0 mismatches of 84

#### lattice check (a p off the 1/200 lattice means that column lost draws)

| element | reps on lattice | n | off-lattice examples |
|---|---|---|---|
| main_att | 84 | 84 | [] |
| pretreatment_att | 84 | 84 | [] |
| ES_att_pre2 | 84 | 84 | [] |
| ES_att_pre1 | 84 | 84 | [] |
| ES_att_post0 | 84 | 84 | [] |
| ES_att_post1 | 84 | 84 | [] |
| e(joint_horiz) | 84 | 84 | [] |

### 5.2 Package verification of the alternative p-value conventions

Every VERIFY-th D1 rep was re-run through the package with permabsp / permrawp /
permside(left) on the SAME rseed (identical draws). The package p is compared with
the value the driver computed from e(permdist_b) + e(perm_comparator_b).

| arm | n verified | exact matches | max abs diff |
|---|---|---|---|
| absp | 4 | 4 | 0.000e+00 |
| rawp | 4 | 4 | 0.000e+00 |
| left | 4 | 4 | 0.000e+00 |

### 5.3 D1 - size under the sharp null (tau = 0)

Nominal size is attained EXACTLY on the 199-replication grid: with the observed
statistic ranked R uniform on 1..200, the equal-tailed p is 2*min(R,201-R)/200, so
P(p<=.05)=10/200=.05 and P(p<=.10)=20/200=.10; the same holds for permabsp,
permrawp and permside(). Any departure below is a real calibration finding.


#### equal-tailed p, alpha = 0.05

| element | size | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| main_att @0.05 | 0.0952 | 0.0320 | 84 | +1.90 |  |
| pretreatment_att @0.05 | 0.0476 | 0.0232 | 84 | -0.10 |  |
| ES_att_pre2 @0.05 | 0.0595 | 0.0258 | 84 | +0.40 |  |
| ES_att_pre1 @0.05 | 0.0000 | 0.0000 | 84 | -2.10 |  |
| ES_att_post0 @0.05 | 0.0714 | 0.0281 | 84 | +0.90 |  |
| ES_att_post1 @0.05 | 0.0357 | 0.0202 | 84 | -0.60 |  |

#### equal-tailed p, alpha = 0.10

| element | size | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| main_att @0.10 | 0.1310 | 0.0368 | 84 | +0.95 |  |
| pretreatment_att @0.10 | 0.0595 | 0.0258 | 84 | -1.24 |  |
| ES_att_pre2 @0.10 | 0.1310 | 0.0368 | 84 | +0.95 |  |
| ES_att_pre1 @0.10 | 0.1190 | 0.0353 | 84 | +0.58 |  |
| ES_att_post0 @0.10 | 0.1190 | 0.0353 | 84 | +0.58 |  |
| ES_att_post1 @0.10 | 0.1548 | 0.0395 | 84 | +1.67 |  |

#### one-sided shares e(perm_p_lo) / e(perm_p_hi), alpha = .05

| element / side | size | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| main_att p_lo @.05 | 0.0714 | 0.0281 | 84 | +0.90 |  |
| main_att p_hi @.05 | 0.0595 | 0.0258 | 84 | +0.40 |  |
| pretreatment_att p_lo @.05 | 0.0357 | 0.0202 | 84 | -0.60 |  |
| pretreatment_att p_hi @.05 | 0.0238 | 0.0166 | 84 | -1.10 |  |
| ES_att_pre2 p_lo @.05 | 0.0833 | 0.0302 | 84 | +1.40 |  |
| ES_att_pre2 p_hi @.05 | 0.0476 | 0.0232 | 84 | -0.10 |  |
| ES_att_pre1 p_lo @.05 | 0.0833 | 0.0302 | 84 | +1.40 |  |
| ES_att_pre1 p_hi @.05 | 0.0357 | 0.0202 | 84 | -0.60 |  |
| ES_att_post0 p_lo @.05 | 0.0714 | 0.0281 | 84 | +0.90 |  |
| ES_att_post0 p_hi @.05 | 0.0476 | 0.0232 | 84 | -0.10 |  |
| ES_att_post1 p_lo @.05 | 0.1190 | 0.0353 | 84 | +2.90 |  |
| ES_att_post1 p_hi @.05 | 0.0357 | 0.0202 | 84 | -0.60 |  |

#### alternative conventions, main_att (computed on every D1 rep)

| convention | size | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| equal-tailed (default) @0.05 | 0.0952 | 0.0320 | 84 | +1.90 |  |
| permabsp @0.05 | 0.1071 | 0.0337 | 84 | +2.40 |  |
| permrawp @0.05 | 0.1071 | 0.0337 | 84 | +2.40 |  |
| permside(left) @0.05 | 0.0714 | 0.0281 | 84 | +0.90 |  |
| equal-tailed (default) @0.10 | 0.1310 | 0.0368 | 84 | +0.95 |  |
| permabsp @0.10 | 0.1429 | 0.0382 | 84 | +1.31 |  |
| permrawp @0.10 | 0.1429 | 0.0382 | 84 | +1.31 |  |
| permside(left) @0.10 | 0.1310 | 0.0368 | 84 | +0.95 |  |

#### mean p-value (global check; theory: equal-tailed .5050, one-sided .5025, SD .2887)

| element | mean equal-tailed p | MC se | dev (MC se) | mean p_lo | n |
|---|---|---|---|---|---|
| main_att | 0.5173 | 0.0315 | +0.39 | 0.4496 | 84 |
| pretreatment_att | 0.5433 | 0.0315 | +1.22 | 0.4606 | 84 |
| ES_att_pre2 | 0.4748 | 0.0315 | -0.96 | 0.4743 | 84 |
| ES_att_pre1 | 0.5063 | 0.0315 | +0.04 | 0.4760 | 84 |
| ES_att_post0 | 0.5170 | 0.0315 | +0.38 | 0.4734 | 84 |
| ES_att_post1 | 0.4470 | 0.0315 | -1.84 | 0.4397 | 84 |

#### fixed-metric joint pre-trend test e(joint_horiz) (jointtests)

| test | size | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| joint_horiz @0.05 | 0.0476 | 0.0232 | 84 | -0.10 |  |
| joint_horiz @0.10 | 0.0833 | 0.0302 | 84 | -0.51 |  |

### 5.4 D1 - multiplicity over the post-treatment horizons (e(perm_stepdown))

Two post horizons (post0, post1). FWER = share of reps rejecting ANY horizon.

| quantity | share | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| FWER step-down min(p_stepdown)<=0.05 | 0.0952 | 0.0320 | 84 | +1.90 |  |
| FWER single-step min(p_singlestep)<=0.05 | 0.0952 | 0.0320 | 84 | +1.90 |  |
| UNCORRECTED min(p_raw)<=0.05 | 0.1190 | 0.0353 | 84 |  | no nominal: this is the multiplicity inflation being corrected |
| per-horizon p_raw post0<=0.05 | 0.0714 | 0.0281 | 84 | +0.90 |  |
| per-horizon p_raw post1<=0.05 | 0.0595 | 0.0258 | 84 | +0.40 |  |
| per-horizon p_stepdown post0<=0.05 | 0.0595 | 0.0258 | 84 |  | per-horizon adjusted |
| per-horizon p_stepdown post1<=0.05 | 0.0357 | 0.0202 | 84 |  | per-horizon adjusted |

| quantity | share | MC se | n | dev (MC se) |  |
|---|---|---|---|---|---|
| FWER step-down min(p_stepdown)<=0.10 | 0.1310 | 0.0368 | 84 | +0.95 |  |
| FWER single-step min(p_singlestep)<=0.10 | 0.1310 | 0.0368 | 84 | +0.95 |  |
| UNCORRECTED min(p_raw)<=0.10 | 0.2381 | 0.0465 | 84 |  | no nominal: this is the multiplicity inflation being corrected |
| per-horizon p_raw post0<=0.10 | 0.1190 | 0.0353 | 84 | +0.58 |  |
| per-horizon p_raw post1<=0.10 | 0.1429 | 0.0382 | 84 | +1.31 |  |
| per-horizon p_stepdown post0<=0.10 | 0.0714 | 0.0281 | 84 |  | per-horizon adjusted |
| per-horizon p_stepdown post1<=0.10 | 0.0714 | 0.0281 | 84 |  | per-horizon adjusted |

Cross-check: p_raw (step-down, |t| share) vs the package permabsp p for the same
horizon on the verification reps - studentizing by one constant preserves ranks,
so these must agree exactly.

| horizon | n verified | exact matches | max abs diff |
|---|---|---|---|
| ES_att_post0 | 4 | 4 | 0.000e+00 |
| ES_att_post1 | 4 | 4 | 0.000e+00 |

#### uniformity of the equal-tailed p on 10 equal-probability bins (chi-square, df=9)

| element | chi2 | n | p>.01 crit 21.67? |
|---|---|---|---|
| main_att | 8.38 | 84 | ok |
| pretreatment_att | 19.10 | 84 | ok |
| ES_att_pre2 | 14.57 | 84 | ok |
| ES_att_pre1 | 5.52 | 84 | ok |
| ES_att_post0 | 4.81 | 84 | ok |
| ES_att_post1 | 5.52 | 84 | ok |

### 5.5 Shift-interval coverage e(perm_ci) at level(95)


#### tau=0 (D1)  (n reps = 84, true effect on every post horizon = 0.0)

| column | coverage of tau | MC se | n non-blank | dev (MC se) | mean width; blank share |
|---|---|---|---|---|---|
| tau=0 (D1) main_att | 0.9048 | 0.0320 | 84 | -1.90 | width 1.664; blank 0/84 = 0.000 |
| tau=0 (D1) ES_att_post0 | 0.9286 | 0.0281 | 84 | -0.90 | width 2.301; blank 0/84 = 0.000 |
| tau=0 (D1) ES_att_post1 | 0.9643 | 0.0202 | 84 | +0.60 | width 2.326; blank 0/84 = 0.000 |

#### tau=1 (D2, default)  (n reps = 40, true effect on every post horizon = 1.0)

| column | coverage of tau | MC se | n non-blank | dev (MC se) | mean width; blank share |
|---|---|---|---|---|---|
| tau=1 (D2, default) main_att | 0.9250 | 0.0416 | 40 | -0.73 | width 1.939; blank 0/40 = 0.000 |
| tau=1 (D2, default) ES_att_post0 | 0.9750 | 0.0247 | 40 | +0.73 | width 2.462; blank 0/40 = 0.000 |
| tau=1 (D2, default) ES_att_post1 | 0.9000 | 0.0474 | 40 | -1.45 | width 2.454; blank 0/40 = 0.000 |

#### tau=3 (D3, default)  (n reps = 16, true effect on every post horizon = 3.0)

| column | coverage of tau | MC se | n non-blank | dev (MC se) | mean width; blank share |
|---|---|---|---|---|---|
| tau=3 (D3, default) main_att | 1.0000 | 0.0000 | 16 | +0.92 | width 3.561; blank 0/16 = 0.000 |
| tau=3 (D3, default) ES_att_post0 | 1.0000 | 0.0000 | 16 | +0.92 | width 3.567; blank 0/16 = 0.000 |
| tau=3 (D3, default) ES_att_post1 | 1.0000 | 0.0000 | 16 | +0.92 | width 3.973; blank 0/16 = 0.000 |

#### tau=1 (D4, permtstat)  (n reps = 24, true effect on every post horizon = 1.0)

| column | coverage of tau | MC se | n non-blank | dev (MC se) | mean width; blank share |
|---|---|---|---|---|---|
| tau=1 (D4, permtstat) main_att | 0.9167 | 0.0564 | 24 | -0.75 | width 2.025; blank 0/24 = 0.000 |
| tau=1 (D4, permtstat) ES_att_post0 | 0.8750 | 0.0675 | 24 | -1.69 | width 2.813; blank 0/24 = 0.000 |
| tau=1 (D4, permtstat) ES_att_post1 | 0.8750 | 0.0675 | 24 | -1.69 | width 2.766; blank 0/24 = 0.000 |

#### paired width growth (same DGP seed and rseed across tau; ratio per rep)

| comparison | column | n paired | mean ratio | median ratio |
|---|---|---|---|---|
| tau=1 / tau=0 | main | 40 | 1.209 | 1.211 |
| tau=1 / tau=0 | post0 | 40 | 1.079 | 1.080 |
| tau=1 / tau=0 | post1 | 40 | 1.075 | 1.096 |
| tau=3 / tau=0 | main | 16 | 2.380 | 2.563 |
| tau=3 / tau=0 | post0 | 16 | 1.620 | 1.578 |
| tau=3 / tau=0 | post1 | 16 | 1.744 | 1.758 |
| tau=1 permtstat / tau=0 | main | 24 | 1.275 | 1.223 |
| tau=1 permtstat / tau=0 | post0 | 24 | 1.249 | 1.239 |
| tau=1 permtstat / tau=0 | post1 | 24 | 1.150 | 1.124 |
| tau=1 permtstat / tau=1 default | main | 24 | 1.101 | 1.057 |
| tau=1 permtstat / tau=1 default | post0 | 24 | 1.196 | 1.188 |
| tau=1 permtstat / tau=1 default | post1 | 24 | 1.111 | 1.114 |

#### scale of the permutation null distribution e(perm_sd)

| phase / arm | mean sd main_att | mean sd post0 | mean sd post1 | n |
|---|---|---|---|---|
| D1 / main | 0.425 | 0.586 | 0.594 | 84 |
| D2 / main | 0.495 | 0.623 | 0.634 | 40 |
| D3 / main | 0.904 | 0.918 | 1.008 | 16 |
| D4 / tstat | 1.426 | 1.672 | 1.552 | 24 |

D4's row is on the t-ratio scale (permtstat); the others are in outcome units.

#### permtstat health under cluster(i)

- main: analytic SE missing or non-positive on 0 of 24 reps
- post0: analytic SE missing or non-positive on 0 of 24 reps
- post1: analytic SE missing or non-positive on 0 of 24 reps

### 5.6 Pre-treatment p-values under a planted post-treatment effect

| phase (tau) | mean p ES_att_pre2 | mean p ES_att_pre1 | mean p pretreatment_att | n |
|---|---|---|---|---|
| D1 (tau=0) | 0.475 | 0.506 | 0.543 | 84 |
| D2 (tau=1) | 0.445 | 0.517 | 0.548 | 40 |
| D3 (tau=3) | 0.562 | 0.669 | 0.750 | 16 |

### 5.7 Surfaces off nominal by more than 2 MC standard errors

- ES_att_pre1 @0.05: 0.0000 vs nominal 0.050 (-2.10 MC se, n=84)
- ES_att_post1 p_lo @.05: 0.1190 vs nominal 0.050 (+2.90 MC se, n=84)
- permabsp @0.05: 0.1071 vs nominal 0.050 (+2.40 MC se, n=84)
- permrawp @0.05: 0.1071 vs nominal 0.050 (+2.40 MC se, n=84)

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
