# BR-4, LENS 0 (independent reproduction) — VERDICT: CONFIRMED, severity high

Built from the CLAIM ALONE. Own DGP, own geometry, own three exclusion routes.
Did not open either prior agent's do-file before my log was written.

## Falsification condition, stated before running
- REFUTED  if `e(boot_df) == e(N_clust)-1` in BOTH arms, or both cautions fire in both arms.
- CONFIRMED if `e(N_clust)` is EQUAL across arms and `e(boot_df)` is NOT.
- PARTIAL  if `e(N_clust)` itself moves.

## DGP (mine)
9 "counties", 2010-2017 (72 rows). 3 treated 2013, 2 treated 2015, 4 never treated.
=> 9 resampling clusters, 5 treated clusters: BOTH documented cautions (<10 total,
<6 treated) must fire on the clean file. `wooldid y id year onset, att cluster(id)
controls(x) bootstrap(standard) reps(60) rseed(777)`.
Five extra units added by three different routes; ARM A leaves them in the file,
ARM B `drop`s them. Nothing else differs.

## Result — all three routes, deterministic quantities (lens0_br4.log)

| route | arm | e(N) | e(N_clust) | e(boot_df) | e(df_r) | e(df_r_analytic) | att (%21x) | cautions |
|---|---|---|---|---|---|---|---|---|
| B1 singleton units (1 period each) | A in  | 72 | 9 | **13** | 13 | 8 | +1.8a7e8ddbfd100X-002 | **NEITHER** |
| B1 | B out | 72 | 9 | **8**  | 8  | 8 | +1.8a7e8ddbfd100X-002 | BOTH |
| B2 outcome missing on every row    | A in  | 72 | 9 | **13** | 13 | 8 | +1.9992e84e17a08X-001 | **NEITHER** |
| B2 | B out | 72 | 9 | **8**  | 8  | 8 | +1.9992e84e17a08X-001 | BOTH |
| B3 controls(x) missing every row   | A in  | 72 | 9 | **13** | 13 | 8 | +1.9992e84e17a08X-001 | **NEITHER** |
| B3 | B out | 72 | 9 | **8**  | 8  | 8 | +1.9992e84e17a08X-001 | BOTH |

CONFIRMED on the pre-stated condition, three times over, three independent routes
into the "excluded from estimation but present in the file" state.
att bit-identical across arms in every route. e(N) and e(N_clust) identical.

## The sharpest single line (NEW — neither prior agent stated it)
One and the same `ereturn` posts

    e(N_clust)         = 9     (estimation clusters)
    e(df_r_analytic)   = 8     (= G-1, computed by the ANALYTIC leg, correct)
    e(boot_df) = e(df_r) = 13  (= G-1 with G = 14, the file's unit count)

wooldid computes the correct resampling-cluster count in one leg of the same run and
the wrong one in the other. `e(df_r_analytic)` is the program's own witness against
`e(boot_df)`; you do not need a second run to see the contradiction.

## Displayed, not merely stored (lens0_br4.log:220 vs :307)
ARM A footer: "e(V) and e(df_r) carry the standard-deviation family referred to the
t(13) distribution ... and so do the tools that read e(V) without being asked: esttab,
coefplot and honestdid".  ARM B footer: the same sentence with t(8).
t(13,.975)=2.1604 vs t(8,.975)=2.3060 -> the Wald-family interval is 6.3% too narrow.

## Mechanism (lens0_mass.log, bootsave() on the ARM A singleton run, reps 120)
- bsv_a.dta carries **14** `w_` columns (w_1..w_9, w_101..w_105) for `e(N_clust)=9`.
  The five phantom singletons get their own resampling columns.
- Phantom mass drawn in **120 of 120** replications (mean 5.08, min 1, max 11).
- Total drawn mass over the nine REAL clusters: mean 8.92, **sd 2.00, min 3, max 13**.
  A correct 9-cluster bootstrap holds that total at exactly 9 in every replication.
- Distinct real clusters present per replication: mean 5.89, min 3, max 9.
=> the skeptic's correction (a) is upheld: the percentile family's DRAW UNIVERSE is
wrong too, not only the df. Magnitude is second order — real mass ~ Bin(14, 9/14),
Var 3.21, so E[G/m] ~ 1 + Var/G^2 = +4% on variance (+2% on SE) at this geometry,
swamped by MC noise at reps in the tens. I did NOT measure a systematic SE gap and do
not claim one: my arm-A SE was LARGER than arm B in B1 (0.260 vs 0.228) and SMALLER in
B2/B3 (0.170 vs 0.210). Direction is noise; the df is not.

## One consequence neither prior agent reported
`e(boot_failed)` = **1** (of 60) in ARM A for routes B2 and B3, and **0** in every ARM B.
The skeptic reported 0 in all four of their routes. Phantom clusters consume draw slots,
so some replications land on too few real/treated clusters to estimate. A second-order
harm, but it shows the phantoms are not inert.

## Files
- /private/tmp/.../boot-reproduce-lens0/lens0_br4.do  + .log   (three routes, two arms each)
- /private/tmp/.../boot-reproduce-lens0/lens0_mass.do + .log   (draw-mass mechanism)
- b1.dta / b2.dta / b3.dta / bsv_a.dta                          (my data, regenerable from the do-file)
Known cosmetic bug in my own block 4: `local wl : list clean "..."` errored r(198)
AFTER `desc, varlist` had already printed the 14 `w_` names; lens0_mass.do covers it.

## Framing
Does NOT contradict known_bydesign §4 (the exact-one `if` rule genuinely does apply
identically everywhere). What it contradicts is the load-bearing code comment at
ado:663-670 — "so the estimation sample, the permutation design and the bootstrap
frame agree" — and sthlp:554, which defines G as "the number of cluster() groups",
against a caution that says clusters "available for the bootstrap". Mechanism: ado:682
and ado:704 zero `ifvar` for missing y/i/t/w and for missing controls; ado:615 forwards
the RAW user `if` as `wdp_if`; bootrun (ado:8917) and bootctl (ado:6461) rebuild their
sample from that, markout never applied.
Adjacent, out of BR-4 scope, flagged not probed: permctl/permrun (ado:5957, 7479) build
their sample the same way, so permute() is likely to carry the same mismatch.

## Post-verdict checks (advisor-prompted)
- **Is the frame a RULING?** `ado:8909` labels it "sample concept (A7)". `grep -rn "A7"
  review_process/` returns **nothing** — A7 is not a documented ruling anywhere in the
  project record, and the only other `A7` in the ado (`ado:5766`) is an unrelated
  zero-variance probe. So this is not a knowingly-ruled frame whose consequences were
  left unreconciled; the ado:663-670 framing above stands as written.
- **Knife-edge, not proportion** (inferred from the mechanism at ado:6461, NOT run):
  `smp2` counts every unit in the file, so at G=9 estimation clusters ONE extra singleton
  (or one unit with a missing control) lifts the count to 10 and kills the <10 caution
  outright; one extra treated-but-unusable unit at 5 treated clusters kills the <6 one.
  The silencing is a threshold crossing, not a function of the phantom share.
- **sthlp:532 prefix-match leg: NOT reproduced here.** I did not run the Stata
  `bootstrap` prefix comparison. It is INFERRED: 14 draw columns against 9 estimation
  clusters means the draw sequence cannot align with a prefix that resamples e(sample),
  so the prober's MAXRELDIF result is necessary rather than merely observed.
- **e(boot_failed) softened.** B1A had the same 14-cluster universe and same rseed with
  0 failures; 1/60 vs 0/60 is consistent with phantoms consuming draw slots but is not
  evidence on its own. Do not cite it as a finding.
