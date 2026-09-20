# hv-contreat-sentinel — NOTES

## Task
Verify recon F5: contreat() ships a 0 point estimate with a -1e100 variance, and the two
exported copies (e(V) vs e(V_es_ame)) disagree.

## Build asserted
`which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado  (version 1.4 2026-08-19). GOOD.

## P1 (p1_repro.do / p1_repro.log) — CLAIM REPRODUCED, 5/5 configs
DGP: 100 units, 5 cohorts (4 treated onset 4/6/8/10, 1 never), T=1..12,
dose = dose0_i * 1{t>=ttre} (the canonical "exposure x post" continuous-DiD design),
y = ufe + 0.2t + 0.8*dose + N(0,1).  i = coh (20 units per cohort cell).

Spec: `att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3)
contreat(dose) contreatcontrols(...)`.

| cfg | contreatcontrols | vce           | e(V) sentinels | V_es_ame zero diag |
|-----|------------------|---------------|----------------|--------------------|
| A   | untreatedZero    | default (iid) | 10 of 40       | 3 of 6 (x2 fams)   |
| B   | basic            | default (iid) | 10 of 40       | 3 of 6             |
| C   | untreatedZero    | customvce(robust) | 10 of 40   | 3 of 6             |
| D   | untreatedZero    | cluster(unit) | 10 of 40       | 3 of 6             |
| E   | basic            | cluster(unit) + fe(unit) | 10 of 40 | 3 of 6          |

Sentinel rows (all b == EXACTLY 0, e(V) diag == -1.0000e+100):
pretreatment_ame, _ame_it, _ame_i, _ame_itime, ES_ame_pre4/pre3/pre2,
ES_ame_it_pre4/pre3/pre2.
The SAME horizons in e(V_es_ame)/e(V_es_ame_it) carry EXACTLY 0.

Both cautions DO fire at run time in all five configs:
 - "Caution: standard errors were not identified for 10 of 40 reported estimates..."
 - "Caution: e(V_es_ame) contains 3 of 6 horizons with zero or missing variance..."
So the asymmetry is announced. Recon lead 1 is REPRODUCED but is announced twice.

## Secondary observation spotted in P1 (self-found, unseeded)
In config A (iid, no cluster at all), ES_ame_post0 / post1 / post2 all carry the
IDENTICAL variance .1320351853 while the att horizons differ
(.1253/.0587/.0939/.0939/.0587/.0939). The help calls "one number repeated across
every ame quantity" the signature of the cell-containment pathology — but says it
arises under cluster() at/above i. Here it appears with NO cluster() at all, and the
AMESE sescan is gated on cluster()/customvce(cluster) (ado:4034), so nothing scans it.
=> follow up in P3.

## Next
P2 — realism: does the sentinel fire when dose VARIES IN EVERY PERIOD?
P3 — identical-ame-variance signature; all-zero e(V) under robust/lattice(all 3).
P4 — consequence surface: ereturn post V_es_ame -> zero-width CI; honestdid; est save.

## P2 (p2_realdose.do) — REALISM GATE
Dose varying in every period: 0 sentinels in 10/10 configs (basic/byTime/byTreated x
iid/robust/cluster(unit), plus untreatedZero). Trigger = pre-treatment dose identically 0.

## P4 (p4_consequence.do) — CONSEQUENCE SURFACE
esttab/lincom/test/estat vce/estimates save+use all render the sentinel honestly.
`ereturn post b_es_ame V_es_ame` -> rows show as `(omitted)`, NOT zero-width CI. Recon's
predicted harm refuted. Internal consumers (joint pre-trend, pretrendpower) correctly withhold.
`test ES_ame_pre4 ES_ame_post0` -> p=.0034 df=1 (one constraint silently dropped in a qui run).

## P5b + repro_ame_degenerate_preblock.do — THE TWO REAL FINDINGS
1. honestdid delta(rm) on e(b_es_ame)/e(V_es_ame) returns the SAME interval at M=0, .5, 1
   under iid/robust/cluster(unit). att family widens normally. delta(sd) dies r(3498).
2. adversarial ame W and q are BITWISE identical computed on the full 6x6 vs the post-only
   3x3: q=2 vs att q=5. Independent reimplementation reproduces both printed p's exactly.

See FINDINGS.md. STATUS: complete.
