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
