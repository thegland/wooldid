# SKEPTIC report -- refutation pass over the "boot-reproduce" findings

All evidence below is from **my own** self-contained do-files in this directory,
on **my own** generated data. I read the original prober's do-files but ran none
of them. Build asserted pinned at `probe_br1.log:59`
(`.../scratchpad/bluesky/build/wooldid.ado`).

| id | verdict | severity (theirs -> mine) |
|----|---------|---------------------------|
| BR-4 | CONFIRMED (one sentence of the EXPECTED corrected) | high -> high |
| BR-2 | CONFIRMED | medium -> medium |
| BR-3 | CONFIRMED (one line-number mis-citation) | low -> low |
| BR-1 | CONFIRMED | low -> low |

Found by me while refuting (seeded = false):

| id | what | severity |
|----|------|----------|
| SK-1 | the RF-6 comma refusal prints a remedy that is itself wrong | low |
| SK-2 | `permute()` shares the same universe defect, where it reaches the randomization design itself | medium |
| SK-3 | sort-order reproducibility of the two legs -- clustered arms CLEAN, unclustered arms below | low |

Negative results (attacked, clean): the documented `e(boot_rngstate)` replay, the
same-`rseed()` repeat, the `bootsave()` dataset characteristics, and the
order-invariance of both clustered legs.

---

## BR-4 -- CONFIRMED. My repro: `probe_br4.do` / `probe_br4.log`

Two arms of the same 9-unit panel. ARM A leaves 6 extra units in the file whose
outcome is missing on **every** row; ARM B deletes them. The estimation sample is
bit-identical across the arms.

```
              e(N)   e(N_clust)  e(boot_df)  e(df_r)   cautions printed
ARM A (in)     90        9          14         14      NONE
ARM B (out)    90        9           8          8      BOTH
att            +1.5e1e5ebcbec84X-001 in both arms (bit-identical)
```

The headline is sharper than the original statement: **one and the same `ereturn`
posts `e(N_clust)=9` and `e(boot_df)=14`.** wooldid already knows the right cluster
count; the bootstrap uses a different one. `e(df_r)` is set to that inflated value
at `ado:10233` ("align e(df_r) with the reference distribution actually used, so
that test and lincom run by hand ... reproduce the reported p-values"), so the
inflated count reaches every hand `test`/`lincom`, `esttab`, `coefplot` and the
honestdid handoff, not only the printed Wald family. t(14,.975)=2.145 against
t(8,.975)=2.306: intervals 7.0% too narrow.

The phantom clusters are not merely counted, they are **drawn**. `probe_br4.log`
BLOCK 2 (30 real units + 10 all-missing-outcome units, `bootsave()`, reps 200):
`e(N_clust)=30`, `e(boot_df)=39`, **40** `w_` columns, phantom unit 35 has mean
multiplicity 1.02 and max 4, and the total drawn mass over the 30 *real* clusters
has mean 29.9, **SD 2.60, min 22, max 36**. A correct cluster bootstrap holds that
total at exactly 30. What runs is an m-out-of-n bootstrap with random
m ~ Binomial(G+P, G/(G+P)).

**Correction to their EXPECTED.** They wrote "The default percentile columns are
rank-based and do not move." That is wrong, and the w_ evidence above is why: the
replication draws themselves come from a different universe, so the percentile
p-value and interval move too. Only the *df* path is percentile-immune. The
systematic size of the move is modest -- E[G/m] = 1 + Var(m)/G^2, about +4.4% on
variance (+2.2% on the SE) at the 9-real/6-phantom geometry -- and at reps(40) it
is swamped by Monte Carlo noise (my arm A se_att 0.170106 vs arm B 0.188565 is a
noise-dominated difference, NOT evidence of the systematic term; I say so rather
than claim the 10%).

The help's own words are on the defect's side: the refusal beside these cautions
is documented as "bootstrap() requires at least two clusters **in the estimation
sample**", and it too counts the user-if universe.

## BR-2 -- CONFIRMED. My repro: `probe_br2.do` / `probe_br2.log`

`sthlp:534`: "wooldid prints a Caution naming that grid whenever it is coarser
than .05" -- unconditional. Measured:

```
A  att cluster(unit) bootstrap(standard) reps(12)        -> caution printed (log:90)
B  att             bootstrap(standard) reps(12)          -> SILENT pre-run
C  att             bootstrap(standard) boottstat reps(4) -> SILENT pre-run
D  att             bootstrap(standard) bootdisplay(wald) reps(6) -> SILENT pre-run
E  att cluster(unit) permute(shuffle) permreps(20)       -> sibling caution printed
```

The string `Caution: with reps(` occurs exactly once in the log, on arm A. Cause is
structural: the whole block at `ado:6455-6478` is inside `if "`clus'" != "" {`
opened at `ado:6456`. Mitigation, which is real: the post-run footer still names
the realized resolution ("The rank-based coefficient p-values have resolution
2/12 = 0.167"). What is lost is RF-2's entire point -- warn *before* the
replications are paid for -- and the "use at least 40 replications" guidance.
Medium stands.

## BR-3 -- CONFIRMED. My repro: `probe_br3.do` / `probe_br3.log`

Compared `e(bootdist)` (the replication coefficient draws), column `main_att`,
between `att` and `att att_i` at one `rseed(20260917)`, reps(120), 40x8 panel:

```
UNCLUSTERED  mreldif(att column) = 0.4505   se_att 0.182663 vs 0.158003
             rep1 count 0.87470193 / reid 0.49623578  (etc.)
CLUSTERED    mreldif(att column) = 1.58e-15 se_att agrees to 2 ulps
```

Mechanism confirmed and it is the author's own: the count path at `ado:9301` notes
"map-frame draws == full-data draws under the cluster()+idcluster() form, measured;
**plain bsample diverges - A17**", and the unclustered branch of the
re-identification leg at `ado:9345` is exactly `if "`clus'" == "" bsample`. So the
divergence was measured and documented in this file, and the one branch that takes
the diverging form is the one whose comment (at `ado:9322`, not `ado:9340` as the
finding cites) asserts "same-seed draws equal the count path's". Same estimator,
Monte Carlo difference only. Low stands.

## BR-1 -- CONFIRMED. My repro: `probe_br1.do` / `probe_br1.log`

```
bootsave(skboot, replace) -> r(198), RF-6 message
permsave(skperm, replace) -> rc 0, "file skperm, replace.dta saved"
confirm file "skperm, replace.dta" -> rc 0      (21 rows: 20 reps + the rep-0 row)
confirm file "skperm.dta"          -> rc 601
```

`ado:5866` even says "permsave(): the name rules are bootsave()'s". They are not.

## SK-1 (NEW, mine) -- the RF-6 comma refusal's own printed remedy is wrong

`probe_br2.do` / `probe_br2.log`. The refusal reads
`... so write bootsave(skboot  replace) replace.` -- built by
`subinstr("`wdbfn'", ",", " ", .)` at `ado:6287`, which turns
`skboot, replace` into the two-token string `skboot  replace` rather than into
the file name. Following that instruction literally:

```
wooldid ..., bootsave(skb2  replace) replace   -> SK1_RC=0
confirm file "skb2.dta"           -> 601
confirm file "skb2 replace.dta"   -> 0
```

So the fix that exists to stop a silently mis-named file tells the user to write
a different silently mis-named file. Same family as BR-1 and in the same commit.

---

## BR-4 realism, measured (`probe_br4c.do` / `probe_br4c.log`)

The all-missing-outcome route is only one way out of `e(sample)`. I tested the two
routes an applied economist hits far more often. The same 9 estimation units in
every route; the 6 extras differ only in *how* they fail to be estimated.

```
ROUTE Z   extras absent                     e(N)=90 e(N_clust)=9 e(boot_df)=8  BOTH cautions
ROUTE S   extras present as SINGLETON units e(N)=90 e(N_clust)=9 e(boot_df)=14 NO cautions
          (one period each; reghdfe drops singletons)
          att bit-identical to Z:  +1.91fcdc8648b4aX-001

ROUTE Zc  extras absent, controls(x)        e(N)=90 e(N_clust)=9 e(boot_df)=8  BOTH cautions
ROUTE C   extras present, controls(x) with  e(N)=90 e(N_clust)=9 e(boot_df)=14 NO cautions
          x missing on every extra row
          att bit-identical to Zc: +1.6e03f12af8532X-001
```

Both realistic routes fire, and in both the extras are invisible: `e(N)`,
`e(N_clust)` and the point estimate are bit-identical to the arm where the extras
are not in the file at all, and nothing in the output says six extra resampling
clusters were invented. A one-period unit in an unbalanced administrative panel,
or a covariate not collected for one state, is not an edge case. This is what
makes BR-4 realistic rather than a curiosity about all-missing outcomes.

## Mechanism, pinned in the source

`ado:682` is where the estimation sample diverges from the bootstrap frame:
```
qui replace `ifvar' = 0 if `i' >= . | `t' >= . | `y' >= . | `w' >= . | `w' <= 0
```
and `ado:704` adds the same for `controls()`, `ccc_absorb()`, `contreat()`,
`interactivecontrols()` and a non-`i` `cluster()`. The estimation runs on
`ifstatement` = `ifvar == 1`. But what is forwarded to the engines is the RAW user
`if`: `ado:615` `c_local wdp_if `"`if'"'`, passed as `ifr()` at `ado:170`
(permrun), `ado:394` (bootrun) and `ado:4927`/`ado:4941` (the two ctl programs),
and rebuilt there as `smp`/`insmp` from `((`ifr') == 1)` alone
(`ado:8917-8925`, `ado:5957`).

The code comment that governs this is `ado:663-670`, and it asserts the opposite:

> "all six forwarded copies of the expression share it (dropcount below,
> wooldid_validate's ifassert, permctl's insmp, permrun's and bootrun's smp), **so
> the estimation sample, the permutation design and the bootstrap frame agree**."

They share the exact-one rule for the user `if`. They do not share the markout that
follows it four lines later. So this is not a documented design decision that a
reviewer missed - the file's own load-bearing comment claims the invariant that is
violated.

## BR-4, why the severity is where it is

The df error only bites when G is small - t(14) vs t(8) matters, t(200) vs t(150)
does not. But the two cautions the defect silences are exactly the small-G
cautions. So the defect is loudest precisely where it matters most: a user with
9 estimation clusters and 6 junk units gets *no* few-cluster warning at all and a
reference distribution 6 df too generous, and the direction is anti-conservative
(t(14,.975)=2.145 vs t(8,.975)=2.306). Under `bootdisplay(wald)` that is the
printed interval and the printed p-value.

Offsetting, and worth stating: the draw universe error runs the other way. The
real clusters' total drawn mass is m ~ Binomial(G+P, G/(G+P)) instead of exactly G,
so conditional on m the replication is an m-out-of-n bootstrap and the variance is
inflated by E[G/m] = 1 + Var(m)/G^2 + O(.) - about +4% on variance at my 9-real /
6-phantom geometry. The two errors do not cancel in any principled way; they are
two independent wrongs in the same engine.

`bootstrap(bayes)` should be immune on the *weights* (the Dirichlet aggregation
property makes the renormalized real-cluster weights Dirichlet again), but its
`e(boot_df)` comes from the same inflated count, so the df half of the defect hits
both modes. I did not measure this; it is an analytic remark, not a result.

## BR-4 is DISPLAYED, not just stored

The run's own footer prints the reference distribution. Same estimation sample,
two arms:

```
ARM A (6 all-missing-outcome units left in file)  probe_br4.log:171
  "e(V) and e(df_r) carry the standard-deviation family referred to the t(14)
   distribution, so test and lincom run by hand reproduce the bootdisplay(wald)
   numbers ... and so do the tools that read e(V) without being asked: esttab,
   coefplot and honestdid export the standard-deviation family"
ARM B (same units deleted)                        probe_br4.log:243
  "... referred to the t(8) distribution ..."
```

Same two lines, t(14) vs t(8), in the singleton and missing-control routes too
(`probe_br4c.log:293` / `:225`, `:452` / `:384`). wooldid tells the user in as many
words which tools will consume the wrong number.

## NEGATIVE RESULT (mine): the documented `e(boot_rngstate)` replay works

`probe_rng.do` / `probe_rng.log`. sthlp (`rseed()`) promises that a run without
`rseed()` can be replayed exactly by `set rngstate` to `e(boot_rngstate)` before
repeating the call. Measured on a 30-unit panel, reps(60), no `rseed()`:
`e(boot_rngstate)` is a 5011-character string, `e(boot_seed)` correctly missing,
and after `set rngstate` the repeat run reproduces the draw matrix with
`mreldif(e(bootdist)_A, e(bootdist)_B) = 0` and a bit-identical `e(V)`
(`+1.2e9f02971741aX-003` both). The promise holds exactly. I attacked it because
it is the one reproducibility claim in this surface that nobody had tested.

## One more place the two counts collide (documentation)

`sthlp:557` (bootsave()) says the `w_` multiplicities sum "to the number of
clusters", and that under `bootalpha(size)` the `_rep = 0` row holds `1/N`, "N the
number of observations in the estimation sample". Both are computed off `smp`, not
off `e(sample)`: in `probe_br4.log` BLOCK 2 the 40 `w_` columns sum to 40 while
`e(N_clust)` is 30. So the saved file is self-consistent with the wrong count and
gives a downstream reader no way to notice.

## SK-2 (NEW, mine) -- permute() shares the same universe, and there it reaches
## the randomization design itself

Static, pinned; empirical confirmation in `probe_permuniv.do` / `.log`.
`wooldid_permrun` builds its cluster map at `ado:7479-7488` from
`((`ifr') == 1)` and `!missing(`clus')` -- the raw user `if` again, with the
`ado:682`/`ado:704` markout never applied. `wooldid_permctl`'s `insmp`
(`ado:5957`) is the same. So a unit with no estimable row is a cluster in the
SHUFFLE: it receives an onset in every replication, and an onset that belongs to
a data-bearing cluster can land on it.

That is worse than the bootstrap case in one specific way. The permutation
reference distribution is supposed to be the randomization distribution of the
OBSERVED design. Here each replication has a random number of data-bearing
treated clusters -- fewer, in expectation, than the observed run has -- so the
null draws are not draws of the observed design's statistic. The support count
`wdp_npd` that drives the "only one possible assignment" refusal (`ado:7547`) and
the `permexact` enumeration is likewise computed over G+P clusters, so `permexact`
enumerates assignments that are duplicates as far as the data are concerned.

Caveat I will not paper over: whether including a data-less unit in the shuffle is
*wrong* depends on whether that unit was part of the assignment mechanism
(attrition -> arguably right) or is a junk row (-> wrong). wooldid cannot know.
What is not defensible is that it is undocumented and that it is inconsistent
with `e(N_clust)` in the same `ereturn`.

## SK-3 (NEW, mine) -- sort-order reproducibility of the two legs

`probe_sort.do` / `.log`: the count path sorts its own map frame
(`ado:9240 qui sort `clusv'`) before drawing, so it should be invariant to the row
order of the user's data at a fixed `rseed()`. The physical and re-identification
legs `frame copy` the user's frame and `bsample` it unsorted (`ado:9246-9248`,
`ado:9344-9346`). Four arms, panel order vs a shuffled copy of the same rows.

### SK-2 measured (`probe_permuniv.log`)

12 real units (4 treated 2004, 4 treated 2006, 4 never) + 6 units whose outcome is
missing on every row (3 with onset 2004, 3 never). `permute(shuffle)`,
`permreps(40)`, `cluster(unit)`, `permsave()`:

```
                     e(N)  e(N_clust)   onset_ columns   att
ARM A (phantoms in)   120      12            18          +1.83bca602bfc4cX-001
ARM B (phantoms out)  120      12            12          +1.83bca602bfc4cX-001
onset_14 (a unit with no data) exists in ARM A and is NON-MISSING -- i.e. it was
assigned a treatment onset -- in 24 of the 40 replications.
```

18 shuffle slots for 12 estimation clusters. 11 of the 18 clusters carry an onset,
so in expectation 6 x 11/18 = 3.67 of the 11 treatment onsets land on a unit with
no data in every replication, leaving about 7.3 data-bearing treated clusters
against the observed run's 8. The null distribution is therefore the randomization
distribution of a *weaker* design than the one the observed statistic came from.
This is systematic, not Monte Carlo.

The reported interval moved: `lb_95ci` for att was **0.2272** in ARM A and
**0.4141** in ARM B (the null band's upper end 0.5301 vs 0.3432) on the same
estimate 0.75729865 and the same `p_perm` 0.049. I will not claim the 82%
difference in that lower bound is all systematic -- at `permreps(40)` a .975
quantile is almost the maximum draw and the two arms also run different RNG
streams. The deterministic part is the 18-vs-12 universe and the 24/40 onsets
assigned to a data-less unit.

### SK-3 measured, CLUSTERED arms (`probe_sort.log`) -- my hypothesis was WRONG here

30x8 panel, `rseed(777001)`, reps(80), panel order vs a shuffled copy of the same
rows:

```
count leg (att),        cluster(unit):  mreldif(e(bootdist)) = 1.641e-15
re-id leg (att att_i),  cluster(unit):  mreldif(att column)  = 1.745e-15
```

Both are order-invariant. `bsample, cluster()` sorts internally, so the unsorted
`frame copy` on the re-identification leg does not cost reproducibility after all.
A clean negative on the clustered half of my own hypothesis. (Cross-check on BR-3:
the two legs' clustered SEs here agree to ~6 ulps, `+1.c5c39cbbea5ecX-003` vs
`+1.c5c39cbbea5e6X-003`, independently reproducing BR-3's clustered leg.)

### SK-3 measured, UNCLUSTERED arms -- here it does fire

```
count leg (att), no cluster():  mreldif(e(bootdist)) = 0.48457418
  se_att 0.214508 (panel order) vs 0.209947 (shuffled) at the SAME rseed(777001)
```

Same rows, same seed, different row order, completely different draws. Cause is
structural and unavoidable in this form: with no `cluster()`, `wooldid_bootrun`
sets the resampling cluster to `_n` (`ado:8921 qui gen double `obsclus' = _n`), so
"cluster 7" means "whatever row is seventh today". `bsample, cluster()` sorts by
that id, so the same draw sequence lands on different observations.

I am NOT calling this a code defect: observation resampling is positional, and
Stata's own `bootstrap` prefix behaves the same way. What is a defect is the
silence. `rseed()` is documented as what makes the run reproducible, the *clustered*
runs really are order-invariant (1.6e-15 and 1.7e-15 above), so a user who learned
the habit on clustered specifications carries a false expectation into the
unclustered ones - and nothing in the help distinguishes the two. Remedy is one
sentence in `rseed()`, not code.

```
re-id leg (att att_i), no cluster():  mreldif(att column) = 0.77309209
  se_att 0.198795 (panel order) vs 0.189499 (shuffled) at rseed(777001)
```

Both unclustered legs move; both clustered legs do not.

## MORE NEGATIVE RESULTS (mine)

`probe_save.do` / `probe_save.log` -- the three records of one run agree.
`att att_i cluster(unit) bootstrap(standard) reps(50) rseed(90210) bootsave(sksv)`:
`e(bootdist)` is 50 x 2 named `main_att main_att_i`; the file has 51 rows,
`_rep` 0..50, the `w_` block, and `_dta[wooldid_engine]=bootstrap`,
`_dta[wooldid_seed]=90210`, `_dta[wooldid_rngstate]` present. The `_rep = 0` row
equals `e(b)` (att bit-identical `+1.05fe7226cd4a0X-002`, att_i within 2 ulps) and
the file's replication draws equal `e(bootdist)`'s to at most 1 ulp. Nothing wrong.

`probe_br4c.log` -- `e(boot_failed) = 0` and `e(boot_completed) = 30` in all four
routes. At this geometry and reps(30) the phantom clusters did not push any
replication into failure, so the "realized grid also moves" hypothesis I raised
did not show. Recorded because the check was made, not because it found anything.

## WHAT I DID NOT DO

- I did not independently verify BR-4's sub-claim that the count leg stops
  matching Stata's own `bootstrap` prefix replication for replication. The
  mechanism is clear from Stata's documented default (the prefix drops
  observations outside `e(sample)` before resampling unless `nodrop`), and it is
  the least consequential leg of BR-4, so I spent the budget on the df, caution
  and realism legs instead.
- BR-3's Monte-Carlo-not-bias conclusion: I reproduced the mechanism, the
  clustered agreement and the size of the gap, but the discriminating test
  (same-leg seed-to-seed spread at reps(300)) is the original prober's probe8. I
  read it; I did not re-run it. My single-seed numbers are consistent with it.
- SK-3's "Stata's own prefix behaves the same way" is by construction (unclustered
  `bsample` draws row positions), not something I measured.
