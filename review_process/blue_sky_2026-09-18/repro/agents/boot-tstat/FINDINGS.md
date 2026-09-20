# FINDINGS -- boot-tstat prober (boottstat / bootdf() / bca / simultaneousbands)
Method: recompute-from-the-dump + shape invariants.
All probe do-files and .log files are in this directory.

---

## BT-01  simultaneousbands under boottstat produces a band NARROWER than the pointwise interval, which the help says is impossible
- **surface**: `simultaneousbands` + `boottstat` + `bootstrap(standard)`
- **severity**: high   **confidence**: certain   **realistic**: yes   **seeded**: yes (lead (a))
- **repro**: `repro_band_narrower.do`
- **claim**: At `level(95)` the simultaneous band for an event-study family can be
  strictly INSIDE the pointwise percentile-t interval at one of that family's horizons.
- **observed** (p04_bands.do, run S3; G=20 states, skewed log-normal outcome error,
  6 horizons, reps(399)):
  ```
  horizon pre3   pointwise percentile-t  [-9.64223368,  0.65471983]
                 simultaneous band       [-6.96720703,  4.39249989]
  ```
  The band's lower endpoint is 2.68 units INSIDE the pointwise one.
  `e(esband_engine) = bootstrap-studentized`, `e(esband_crit_att) = 9.5180674`.
- **expected**: `wooldid.sthlp:525` states, of exactly this engine: *"a genuine
  percentile-t, on the same scale as the pointwise intervals printed beside it,
  **so the band can never come out narrower than the interval it accompanies**"*.
  The same paragraph correctly warns that non-nesting DOES happen under
  `bootdisplay(percentile)/bc/bca` "where the printed pointwise interval is a rank
  interval and not symmetric".
- **mechanism**: the pointwise percentile-t interval is EQUAL-TAILED from the SIGNED
  t-quantiles, `[b - q_.975(t_r)*se, b - q_.025(t_r)*se]` (`wooldid_boottp`, ado:14389;
  wired at ado:9874). The band is SYMMETRIC, `b +/- q_.95(max_h |t_rh|)*se`
  (`wooldid_bootband`, ado:16866-16876). For a symmetric t_r distribution
  `q_.95(|t|) = q_.975(t)` and the band dominates; a SKEWED replication t-distribution
  breaks that equality, exactly as it does for the coefficient-draw families the help
  already flags. The percentile-t pointwise interval is itself a rank interval that is
  not symmetric, so the sentence exempting `boottstat` is unsound as written.
- **corroboration**: on a second design (p09, G=12, one high-effect cluster,
  non-studentized) the margin is razor thin in the same direction:
  `(q.975 - b)/sd = 2.1702` against `c = 2.1738` at horizon post0.

---

## BT-02  bootdf(clusters) counts clusters that contribute NOTHING to e(sample); the same printout says "N clusters = 17" and "t(19)"
- **surface**: `bootdf(clusters)` / the bootstrap cluster map frame
- **severity**: medium   **confidence**: certain   **realistic**: yes   **seeded**: yes
- **repro**: `repro_bootdf_deadcluster.do`
- **claim**: `G` for `bootdf(clusters)` is the number of clusters with a non-missing
  `cluster()` value inside the user `if`, NOT the number of clusters in the estimation
  sample. A cluster whose rows are all dropped by markout (missing outcome, missing
  control) is counted in the degrees of freedom, and is also drawn by `bsample` as an
  empty slot.
- **observed** (p03_df_nosub.do; 20 states, outcome suppressed in 3 of them -- an
  ordinary administrative shape):
  ```
  A (3 suppressed clusters left in the data): e(N_clust)=17  e(boot_df)=e(df_r)=19
  B (same estimation sample, 3 states dropped): e(N_clust)=17  e(boot_df)=e(df_r)=16
  identical point estimate in A and B.
  ```
  The printed table header of run A says `N clusters = 17` while the footer of the SAME
  run says `e(V) and e(df_r) carry ... the t(19) distribution`. invt(19,.975)=2.093
  against invt(16,.975)=2.120, so every Wald p-value and normal-theory interval in the
  run is referred to the wrong distribution, and a hand `test`/`lincom` inherits it.
- **expected**: `sthlp:554` -- *"bootdf(clusters) instead uses the t distribution with
  G-1 degrees of freedom, where G is the number of cluster() groups"*, and
  *"e(df_r) is aligned with the reference distribution actually reported"*. The run's own
  `e(N_clust)` is 17.
- **where**: `wooldid_bootrun` ado:8911-8955. `smp` is `(ifr==1) & !missing(clusv)` only;
  `wdb_C` is counted on it (ado:8954) and returned as `wdb_r_nclus` (ado:9537), becoming
  `nclus()` in `wooldid_bootwire` and hence `wdbdf = nclus - 1` (ado:9573). No markout
  against the estimation sample anywhere on that path.
- **RETRACTED sub-claim (measured and withdrawn)**: I first read a +4.4% bootstrap-SE gap
  between A and B at reps(99). At reps(499) over three seeds the ratio A/B is
  0.97312 / 0.99366 / 1.04714, mean 1.00464, sign not consistent -- that gap was
  Monte-Carlo noise and the standard error itself is NOT materially inflated by the dead
  clusters. What does survive: the dead clusters are still DRAWN, so 3, 2 and 2 of 499
  replications failed outright across the three seeds in A (B's count is 0 on the same
  design at reps(99)). The df leg is the finding; the resampling waste is a note.
- **not a defect**: under `boottstat` (the re-identification leg) the P7 internal-
  consistency assertion does NOT misfire on this data (p07 T1: rc 0, 98/99 complete).

---

## BT-03  under bootdisplay(wald) -- where the pointwise interval IS symmetric -- the simultaneous band still sits inside it, which the help's non-nesting caveat explicitly does not cover
- **surface**: `simultaneousbands` + `bootdisplay(wald)` + `bootdf()`
- **severity**: medium   **confidence**: certain   **realistic**: yes   **seeded**: no
- **repro**: `repro_band_ignores_bootdf.do` (sweep: p11_banddf.do)
- **claim**: `sthlp:525` scopes the non-nesting warning to *"where the printed pointwise
  interval is a rank interval and not symmetric - under bootstrap() with the default
  bootdisplay(percentile), or with bootdisplay(bc) or bootdisplay(bca)"*. Under
  `bootdisplay(wald)` the pointwise interval is `b +/- zc*sd`, symmetric, so on the help's
  own logic the band must nest it. It does not.
- **observed**, two independent mechanisms:
  1. *the band's critical value never sees `bootdf()`*. G=8, 2 horizons, reps(999):
     `e(boot_df)=7`, `invt(7,.975)=2.364624`, `e(esband_crit_att)=2.108943`; the band is
     inside the pointwise interval on BOTH sides at BOTH horizons. Same at G=12 (p09 K3):
     2.1738 against `invt(11,.975)=2.2010`, both sides, both horizons. `bootdf(clusters)`
     is the DEFAULT whenever `cluster()` is given.
  2. *different replication subsets*. p11 cell `W_G8_M3`: `c = 2.382340 > zc = 2.364624`,
     yet 4 of 6 horizons still violate, because the band's `sd_h` is computed on the
     LISTWISE-complete replications while the table's `se_h` is computed per column
     (`wooldid_bootfin` vs `wooldid_bootband`).
- **where**: `zc` at ado:9576 uses `wdbdf`; `wooldid_esband` (ado:10879) passes only
  `level()` to `wooldid_bootband`, whose critical value is `wooldid_bootqs(M, lev/100)`
  (ado:16866) with no degrees-of-freedom correction and no shared replication mask.
- **note**: `sthlp:554` does say `bootdf()` leaves "quantile-based" objects unaffected, so
  the band's insensitivity to `bootdf()` is arguably intended. The divergence is the
  under-scoped nesting caveat at `sthlp:525` plus the fact that a band narrower than the
  pointwise interval is not simultaneous at `level()` under any display.

---

## BT-04  bootsave()'s `_rep = 0` row for `jointW_*` holds an event-study COEFFICIENT, not the observed joint statistic
- **surface**: `bootsave()` + `boottstat` + `jointtests` (and the `advW*` columns)
- **severity**: medium   **confidence**: certain   **realistic**: yes   **seeded**: no (found here)
- **repro**: `repro_bootsave_jointW.do`
- **claim**: the observed-statistic row of the `jointW_*` columns is filled by POSITION
  from a matrix with a completely different column set, so it holds an observed
  pre-treatment event-study coefficient. Anybody who recomputes the joint p-value from
  the saved file, which is exactly what the file is sold for, gets a silently wrong
  number.
- **observed**:
  ```
  e(joint_horiz) reported by wooldid          = 0.682119205
  n replications with a usable W_r            = 150
  implied #{W_r >= W_obs}                     = 102  -> W_obs in (16.648, 17.128]
  bootsave _rep==0 value of jointW_att        = -0.33935852     <-- NEGATIVE
  min / max of the replication W_r            = 0.4055 / 8896.1
  -0.33935852 is exactly e(b_es_att)["pre3"]
  p recomputed from the file as documented    = 1.000000000
  ```
  A joint W statistic is a quadratic form and cannot be negative, so the wrong value is
  detectable in principle - but nothing in the file or the run says so.
- **expected**: `sthlp` `bootsave()` entry: *"Under boottstat the replication joint-test
  statistics are added as jointW_name and the adversarial statistic as advW, **with the
  observed statistics on the _rep = 0 row**"*, and the same entry sells the file as
  self-sufficient ("so the draws can be recentred from the file alone").
- **where**: `wooldid_bootsavemk` ado:14915-14923.
  `JW = e(boot_jointW)` has one column per FAMILY (`wdb_jfk`: att, att_it, ame, ame_it
  and `_s#` twins). `jobs` is built by `wooldid_bootjobs` (ado:14666) and has one column
  per (family, PRE-HORIZON) pair holding `e(b_es_<fam>)[1,j]`. The loop writes
  `JO[1, k]` under `JW`'s column index `k`. The engine itself never uses `JO` that way:
  `wooldid_bootjnt` gets the observed statistic through the NAME-keyed
  `wooldid_bootjwobs(fk, jobs, wobs)` (ado:15310).
- **scope**: the `jointW_*` leg is MEASURED. The `advW*` block six lines below
  (ado:14925-14931) has the identical positional shape (`aobs` is the same kind of
  observed-coefficient vector, built by the same `wooldid_bootjobs` with `allh=1`, while
  `awacc`'s columns are `<family>_a0`/`<family>_a1`), so the same defect is expected
  there -- but I did NOT run `adversarial`, so that half is READ, not measured.

---

## BT-05  with cluster() at or above `i`, every subgroups() row's analytic standard error is NUMERICALLY IDENTICAL to the pooled one
- **surface**: `subgroups()` analytic SE -> `e(boot_ase_obs)`, the `boottstat`
  studentization denominator and the scale of every printed percentile-t interval
- **severity**: high   **confidence**: probable   **realistic**: yes   **seeded**: no
- **repro**: `repro_sgse_clusteri.do` (and p12_sgse.do, p13_sgse2.do)
- **claim**: under a clustered VCE at or above the level of `i`, the standard error
  reported for each `subgroups()` level equals the pooled standard error to 15
  significant digits, on the att row and on every event-study row. The subgroup standard
  errors therefore carry no subgroup-specific information, while the estimates do.
- **observed**:
  ```
  repro_sgse_clusteri.do   24 states, i = state, cluster(state), subgroups(mod(state,2))
    main att b=1.567921674  se=0.126889554  (+1.03deab9e44e73X-003)
    sg0  att b=1.496542253  se=0.126889554  (+1.03deab9e44e75X-003)
    sg1  att b=1.639301096  se=0.126889554  (+1.03deab9e44e74X-003)
    ratios to pooled: 1.00000000000  1.00000000000   (5 of 6 table rows identical)
  CONTROL, same data, robust instead of cluster():
    main 0.226363662   sg0 0.267818831   sg1 0.244313253   ratios 1.183 / 1.079
  CONTROL, two separate runs restricted to each subgroup, cluster(state):
    reg2==0 alone b=1.465722071 se=0.221480605
    reg2==1 alone b=1.670121277 se=0.132589734
  p13_sgse2.do   i = district (60 cohorts), cluster(state) STRICTLY ABOVE i:
    main 1.153391171 se=+1.644d7ae0eb146X-004
    sg0  1.169219231 se=+1.644d7ae0eb141X-004
    sg1  1.137563111 se=+1.644d7ae0eb14aX-004    ratios 1.00000000000
  ```
- **expected**: two subgroups holding different cohorts, with different point estimates
  and (as the separate runs show) genuinely different sampling variability, cannot have
  the same standard error to the last two bits. The `robust` control uses the SAME weight
  vectors and gives three different numbers, so the collapse is a property of the
  clustered `e(V)`, not of the estimand's weights.
- **why it is on this surface**: `wooldid_bootwire`'s first walk (ado:9606-9625) fills
  `aseobs` from column 2 of `e(wdidmainresults)` and `e(wdidsgresults_*)` -- exactly these
  numbers -- and `aseobs` is `e(boot_ase_obs)`, the `boottstat` studentization denominator
  and the multiplier in `b - q*se_obs`. Percentile-t is partly protected (`se_obs` and
  `se_r` are the same functional, so the t-ratio is still comparable); the Wald and
  normal families, and the `t_pt` column itself, are not.
- **relation to what is already ruled**: known_bydesign section 1 rules that analytic SEs
  are UNDERSTATED (0.66-0.75 coverage) when `cluster()` is at or above `i`. That entry is
  about the LEVEL of the pooled standard error. It does not say the subgroup standard
  errors become numerically identical to it. I have not excluded a structural explanation
  (a clustered `e(V)` that is proportional to a constant matrix on the cell-dummy block
  would make every weight vector summing to one give the same variance), which is why the
  confidence is `probable` rather than `certain`.

---

---

## CLEAN (probed hard, nothing found)
See NOTES.md for detail. Summary:
- **percentile-t kernel**: recomputed n, t_obs, p, lb, ub for every coefficient from a
  `bootsave()` dump: BIT-EXACT (`RECOMP_MAXRD n=0 t=0 p=0 lb=0 ub=0`). Add-one
  `2*(min+1)/(n+1)` confirmed on studentized COEFFICIENT rows.
- **BCa**: `e(boot_accel)` reproduced exactly by an independent delete-one-cluster
  jackknife of my own; z0, both endpoints and the dual p reproduced exactly; jackknife
  count == G at G=20 and G=8; estimate inside the interval on every row.
- **Refusals**: `bootdisplay(bc|wald|bca)+boottstat` and `bootdisplay(bca)` without
  `cluster()` all r(198) with informative messages.
- **Analytic band family independence (D-9)**: `e(esband_crit_att)` and the whole band
  are bit-identical when `att_it` or `subgroups()` is added to the run.
- **Dump shape invariants**: `rowtotal(w_*) == G` on every replication; `#se_ vars ==
  #coefficient vars`; dump `_rep==0` == `e(b)` and `e(boot_ase_obs)` exactly.
- **WOOLDID_BOOTNOSUB**: does what its comment says (se_r/se_obs rms 0.927 -> 1.204).
- **The `kk` column cursor with subgroups()** (the hazard the code itself flags at
  ado:9596-9605): recomputed the studentized se/t/p/lb/ub for all 18 columns of a
  `subgroups()` run (6 main + 6 sg0 + 6 sg1) from the dump. `SGWALK_WORST se=0 t=0 p=0
  lb=0 ub=0` -- the subgroup block is aligned, no off-by-one. (p10_sgwalk.do)
- **Studentization legs** (p08_legs.do), mean(se_r)/se_obs over replications:
  unclustered fweight leg 0.945-0.966, unclustered+aweight phys leg 0.953-0.971,
  clustered reid leg rms 0.915-0.972 (the help predicts ~0.92 +/- 0.04).
  `boottstat + bootstrap(bayes)` is refused r(198). `bootdf(0)`, `bootdf(-5)`,
  `bootdf(2.5)` all refused; `bootdf(1e9)` is accepted and posts `e(df_r)=1e9`
  (absurd input, not filed).
- **P7 internal-consistency assertion**: does NOT misfire when clusters contribute no
  rows to `e(sample)` under the re-identification leg (p07 T1: rc 0, 98/99 complete).

---

---

## UNTESTED / NOT REACHED
- Seeded lead 5(c): the band-refusal path that is supposed to NAME the offending horizon
  was never triggered. My rare-horizon design (espostlength(5)) produced post0..post4,
  complete in every replication (e(esband_n_att) = e(boot_completed) = 393). Reading the
  code, the `fewer than two replications are complete over every horizon of the family'
  branch (ado:16836) is the one that does NOT name a horizon, contra the D-8 comment
  five lines above it; that is READ, not RUN.
- BCa's BS-17 blanking path (`1 - a*w <= 0') was not provoked: at G=8 the accelerations
  were 0.01-0.05 and every family was computed.
- The `advW*' half of BT-04 is read, not measured (see the scope note there).
- Multi-way cluster() with bootstrap(): the help says it is refused; not run.
- `bootdf(1e9)' is accepted and posts e(df_r)=1e9. Absurd input, not filed.
- Whether BT-05 also fires with cluster() STRICTLY FINER than i was not tested.
