# F3 - att_i coverage under finer-than-i clustering (prereg P4/P5)

truths (zero-noise dgpE readout): att_i = 0.8, att = 0.953846153846

nominal 0.95; central ~95% binomial region, normal approx rounded outward; MINN=30; MARGINAL corridor (0.9, 0.98)


**P4 - FINE**: cluster(unit), STRICTLY FINER than i=state. rows=150, rc!=0: 0. 
replications completed (min/median/max): 199.0/199.0/199.0; failed total 0; e(N_clust) values seen: [160]

| pred | cfg | family | estimand | covered/n | rate | band | Wilson | verdict |
|---|---|---|---|---|---|---|---|---|
| P4 | FINE | percentile (displayed) | att_i | 142/150 | 0.9467 | (137, 148) | [0.8983, 0.9727] | PASS |
| P4 | FINE | percentile (displayed) | att | 147/150 | 0.9800 | (137, 148) | [0.9429, 0.9932] | REPORTED (control: heterogeneity live) |
| P4 | FINE | normal/Wald | att_i | 145/150 | 0.9667 | (137, 148) | [0.9243, 0.9857] | REPORTED (secondary family) |
| P4 | FINE | normal/Wald | att | 147/150 | 0.9800 | (137, 148) | [0.9429, 0.9932] | REPORTED (control: heterogeneity live) |

**P5 - NEST**: cluster(state) == i, the regression leg. rows=150, rc!=0: 0. 
replications completed (min/median/max): 199.0/199.0/199.0; failed total 0; e(N_clust) values seen: [40]

| pred | cfg | family | estimand | covered/n | rate | band | Wilson | verdict |
|---|---|---|---|---|---|---|---|---|
| P5 | NEST | percentile (displayed) | att_i | 149/150 | 0.9933 | (137, 148) | [0.9632, 0.9988] | FAIL |
| P5 | NEST | percentile (displayed) | att | 149/150 | 0.9933 | (137, 148) | [0.9632, 0.9988] | REPORTED (control: heterogeneity live) |
| P5 | NEST | normal/Wald | att_i | 149/150 | 0.9933 | (137, 148) | [0.9632, 0.9988] | REPORTED (secondary family) |
| P5 | NEST | normal/Wald | att | 149/150 | 0.9933 | (137, 148) | [0.9632, 0.9988] | REPORTED (control: heterogeneity live) |

## FAILURES
- P5/NEST/percentile (displayed)/att_i: 149/150 = 0.9933, band (137, 148)
