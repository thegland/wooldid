# hv-cwm — verification notes (customweightmultiplier vs the cell-weighted estimands)

STATE: COMPLETE. Both halves of the seeded claim settled empirically.

## Files
- `repro_cwm_estimand.do` / `.log` — THE self-contained repro (arms A/B/C + per-cohort weight tables)
- `p1_cwm_estimand.do` / `.log` — 5 arms; raw-data oracle for the cell weights; gteffects identity check
- `p2_cwm_adjacent.do` / `.log` — growing cwm (F), 0/1 subpop cwm (G); H arms failed r(198) (missing contreatcontrols)
- `p3_cwm_collapse.do` / `.log` — fast collapse vs nofastcollapse under cwm (H/K); escomp-vs-gteffects (L)
- `panel*.dta` — generated panels (not project data)

## Mechanism (code, pinned build)
ado:1581 `local customweighttext * `customweightmultiplier''
ado:1619 `W_it  = total(esamp*w)  by(i,t)`      <- NO cwm
ado:1620 `W_i   = total(w*esamp)  by(posttreatment,i)` <- NO cwm
ado:1630/1637/1644/1650 numerators are `w `customweighttext'`.
=> cell total weight = sum_cell(w*cwm)/sum_cell(w). Equals cwm_cell when cwm is cell-constant.

## Measured (repro_cwm_estimand.log)
ARM A (no cwm):   att_it = 2.74890010 == simple mean of gteffects cells; w_att_it in {1}
ARM B (cwm=pop, 50:1 spread):
                  att_it reported          = 3.74381423
                  documented (simple mean) = 2.74890010     -> 36.2% relative error
                  pop-weighted mean        = 3.74381423     -> exact match, mechanism confirmed
                  w_att_it min/max         = 1 / 50
                  att reported             = 3.74124884     -> att_it now 0.07% from att
                  per-cohort w_att_i       = 1,1,50,50 (pre AND post blocks) vs 1,1,1,1 in arm A
ARM C ([aw=sw*pop], same multiplier folded into the SAMPLE weight): w_att_it min/max = 1/1 -> CANCELS

p1 arms: A w_att_it sd 1.2e-16; B min1/max50; C 1/1; D (within-cell cwm) 0.655-1.616 = w-weighted
mean of cwm per cell, matches the raw-data oracle at 2.1e-16; E ([aw=w*mfine], mfine varying
WITHIN cell) 1/1 -> cancellation is exact even for a cluster FINER than the cell.
gteffects identity (sthlp:1136) holds in every arm: rdW <= 4.8e-16.

## ARM D (repro, time-varying cwm = pop*1.15^year) - MEASURED
per-cohort share of its block's total w_att_i:  main .008400/.009057/.479530/.503014
                                        pretreatment .007422/.008119/.463144/.521315
i.e. the pretreatment att_i row weights cohort 1 11.6% less than the main att_i row does.
With no cwm both blocks are exactly 1 per cohort.

## Verdicts
HALF 1 CONFIRMED. HALF 2 REFUTED: the bootctl gate (ado:6430-6442) is right; bootstrap(bayes)
passes its Dirichlet draw as the sample weight (`repw = basew*drawv`, ado:9449-9451), which lands
in W_it as well and cancels exactly (arms C and E). cwm is numerator-only; the Dirichlet draw is
numerator-and-denominator. The recon's either/or is a false dichotomy — both statements hold.

## Engines x cwm (p4_cwm_engines.log) - CLEAN
N0 base (cwm, cluster(county)) -> N1 bootstrap(standard) reps(25): point estimates maxreldif = 0
N2 permute(shuffle) permreps(25): maxreldif = 1.19e-16   (sthlp:532 "point estimates are unchanged")
N3 bootstrap(bayes) + att_it + cwm: refused r(198), the cell-weighted gate fires as designed
N4 bootstrap(bayes) + att only + cwm: point estimate maxreldif = 0

## Clean (probed, nothing found)
- default cell collapse vs `nofastcollapse` under cwm, binary AND contreat(dose)+contreatcontrols(basic),
  including a cwm that is exactly 0 on the cell's FIRST row (the case ado:5636 names as having
  once broken): max reldif 1.3e-13 (contreat) / 3.2e-16 (binary). Comment at ado:5627-5644 verified.
- e(escomp) `wshare_att_it` reproduces the shares rebuilt from e(gteffects) `w_att_it` exactly
  (102/102/102/102 at post0-3, 101 at pre3 with only 3 cohorts).
- the runtime note "Note: all Estimates Shown Apply <var> Weights Post-Estimation" does fire
  (main table ado:11369 and per-subgroup ado:11529).
