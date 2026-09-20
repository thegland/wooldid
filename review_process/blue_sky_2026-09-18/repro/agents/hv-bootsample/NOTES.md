# hv-bootsample — verification of "the bootstrap/permute resampling universe is not e(sample)"

## Code reading (pinned build, ado lines verified)

- `wooldid_bootrun` (ado:8908-8923): `smp = ((ifr)==1)`, then `replace smp = 0 if missing(clusv)`.
  NO markout of y / controls / weights / dose. `ifr` = `wdp_if` = the RAW user `if`
  (c_local at ado:615), NOT the internally screened `ifvar`.
- `ifvar`, the real estimation sample (ado:652-686): raw `if`==1 AND `!missing(i,t,y,w)` AND
  `w>0` AND markout of controls/fe/cluster/dose. So smp is a STRICT SUPERSET of e(sample)'s rows.
- `wooldid_permrun` (ado:7479-7488): byte-identical construction. Same mismatch.
- `wooldid_bootctl` (ado:6459-6501): `smp2 = ((ifr)==1) & !missing(clus)` — the SAME unscreened
  concept — drives THREE guards:
    * r(459) "bootstrap() requires at least 2 clusters **in the estimation sample**"  (ado:6466)
    * caution "only N clusters are available for the bootstrap" if N<10             (ado:6485)
    * caution "only K treated cluster(s) are available for resampling" if K<6        (ado:6498)
- `wdb_C` (ado:8956) = # distinct clusters in smp → `c_local wdb_r_nclus` (ado:9537) →
  `bootdf(clusters)` df = wdb_C - 1 (ado:6600, 6738, 9570) → posted `e(boot_df)` (ado:10233).
  Also the BCa jackknife loop `forval g = 1/wdb_C` (ado:9493) and bootsave's `_rep=0` row
  (1/wdb_r_nclus under equal, 1/wdb_r_nrows under size).
- `__wdbn` (ado:8946) = the cluster's row count IN THE MAP FRAME = Gamma shape AND divisor
  under `bootalpha(size)`. Rows out of e(sample) inflate it.

## Contradicted record
- ado:666 comment: "all six forwarded copies of the expression share it (... permrun's and
  bootrun's smp), so the estimation sample, the permutation design and the bootstrap frame
  AGREE." They agree on the exact-one rule only; they do NOT agree on markout.
- known_bydesign §4 (BS-37) is about the exact-one rule — NOT contradicted.
- ado:9367 shows the author knows the frame ≠ e(sample) ("the replication drops rows the frame
  still holds (missing y, an event-study window, singletons)"). So the G/df/guard consequences
  look like an oversight, not a ruling.

## Oracle design
Control = physically drop the dead rows (`drop if missing(y)`) and re-run identically.
wooldid's own markout means the OBSERVED estimates are identical either way, so ANY difference
in bootstrap/permute output is purely the sample mismatch. Do NOT expect bit-identity of the
replication draws (C changes → the Dirichlet/bsample stream shifts); compare df, guard text,
and SEs statistically (MC noise on a bootstrap SE ≈ SE/sqrt(2R)).

## Theory before measuring
- bayes + bootalpha(equal) + WHOLLY dead clusters: Dirichlet aggregation + aweight scale
  invariance ⇒ the normalized weights on the live clusters are EXACTLY Dir(1^G'). Predict the
  SE is RIGHT. The damage is df / guards / BCa, not the SE.  <-- expected clean negative
- bayes + bootalpha(size) + PARTIALLY dead clusters: shape and divisor are n_total, so the
  weight CV^2 shrinks by the live fraction f ⇒ SE ~ sqrt(f) times too small, ON TOP of the
  shrinkage the help already documents.
- standard (multinomial): wasted draws only ⇒ variance inflation 1+(1-p)/G', negligible.
  Real damage is df/guards + a higher rep failure rate.

## Probe ledger
(see p*.do / p*.log in this directory)

### Probe ledger (final)
| file | what it tested | result | finding |
|---|---|---|---|
| `p1_guards.do` / `.log` | dead clusters (all-y-missing states) vs oracle, bayes + standard | e(boot_df) 19 vs 15 and 11 vs 7; BOTH few-cluster cautions suppressed; identical e(N) and point estimate | HV-BOOT-01 |
| `p2_size_degen.do` / `.log` | A degenerate 1-live-cluster; B `bootalpha(size)` with 50% within-cluster missingness; C bootsave identities; D singleton unit on COMPLETE data | A: caught by a downstream guard (negative); B: SE ratio .684 vs predicted sqrt(f)=.699 — **32% understated**; C: `_rep=0` holds 1/336, e(N)=164; D: boot_df 10 vs 9 with identical e(N) | HV-BOOT-02, -03, -04 |
| `p3_perm_assert.do` / `.log` | A permsave design universe; D validate asserts. **B and C INVALID — see the header of the do-file; superseded by p5** | A: 20 assignment columns vs 16, 634,888,800 vs 43,680 possible assignments, treated-cluster caution suppressed; D: r(411) on the full file, rc 0 on the identical analysis sample | HV-PERM-01, HV-VAL-01 |
| `p4_df_bca.do` / `.log` | A 2-live-of-12 clusters; B `bootalpha(equal)` with wholly dead clusters at reps(600); C BCa jackknife; D ES window | A: boot_df 11 vs 1, both cautions suppressed; B: ratio 1.0165 vs 2.89% noise — **CLEAN**; C: "18 fits" for 13 clusters, effect confounded; D: window does NOT drop rows | HV-BOOT-05, two negatives |
| `p5_perm_refusals.do` / `.log` | A onset heterogeneity among dead rows; B `.a` on a dead row; C permonsetcohorts bijection; D permexact | A/B/C: hard r(459) on the full file, rc 0 once the single row is deleted, ordinary `wooldid` fine either way; D: exact p 0.1000 vs 0.4000, "40 of 42" pruned | HV-PERM-02, -03 |
| `p6_nestregime.do` / `.log` | T2/T3 nesting regime decided on `smp`: one cross-cluster dead row | att_i se .1071 vs .1188, p **0.040 vs 0.093**, CI flips across zero, identical e(N) and estimate | HV-BOOT-06 |
| `p7_nest_seeds.do` / `.log` | is HV-BOOT-06 systematic? four seeds at reps(400) + a real-straddle control | ratios .8135 / .8693 / .8884 / .9033, all < 1; control (stray row given a real y) se .10191 == the full-data .10183, far from the oracle .12517 — the regime flip IS the mechanism | HV-BOOT-06 |

### Verdict
The seeded claim is **TRUE**, mechanism as stated, and the consequence surface is wider than the
claim: three safety guards, bootdf(), the BCa jackknife, bootsave()'s scale, the permutation
design description, the `permexact` reference set, three hard permute refusals, and the
one-day-old T2/T3 nesting regime all read the unscreened sample. The related validate-assert
claim is also TRUE. `bootalpha(equal)`'s SE is the one thing measured CLEAN.
