## Run inventory

| cfg | arm | n rows | rc!=0 | mean secs | reps completed | reps failed | boot_N_V |
|---|---|---|---|---|---|---|---|
| AN | bayes | 23 | 0 | 114.8 | 199.0 | 0.0 | 199.0 |
| AN | plain | 25 | 0 | 1.1 | - | - | - |
| AN | std | 24 | 0 | 63.8 | 199.0 | 0.0 | 199.0 |
| AT | bts | 10 | 0 | 111.3 | 199.0 | 0.0 | 199.0 |
| AT | plain | 10 | 0 | 0.8 | - | - | - |
| BA | bayes | 14 | 0 | 136.2 | 199.0 | 0.0 | 199.0 |
| BA | plain | 15 | 0 | 0.7 | - | - | - |
| BA | stdbca | 15 | 0 | 83.2 | 199.0 | 0.0 | 199.0 |
| BI | attit | 28 | 0 | 76.3 | 198.9 | 0.1 | 198.9 |
| BI | plainit | 29 | 0 | 0.4 | - | - | - |

## Point estimate and the standard-error ratio

(ratio = mean reported SE / SD of the estimate across datasets; 1.00 = calibrated. The SD denominator is that arm's OWN estimates, because the event-study window changes the estimator.)

| cfg | arm | n | truth | mean b | bias | bias t | SD(b) | mean SE | ratio (SE) |
|---|---|---|---|---|---|---|---|---|---|
| AN | bayes | 23 | 0.00000 | 0.1193 | +0.1193 | +1.76 | 0.3244 | 0.3081 | 0.950 (0.148) |
| AN | plain | 25 | 0.00000 | 0.1189 | +0.1189 | +1.76 | 0.3370 | 0.3450 | 1.024 (0.157) |
| AN | std | 24 | 0.00000 | 0.1002 | +0.1002 | +1.48 | 0.3308 | 0.3453 | 1.044 (0.160) |
| AT | bts | 10 | 0.00000 | 0.2118 | +0.2118 | +1.66 | 0.4035 | 0.3324 *coefficient-draw SD | 0.824 (0.202) |
| AT | plain | 10 | 0.00000 | 0.2118 | +0.2118 | +1.66 | 0.4035 | 0.3368 | 0.835 (0.211) |
| BA | bayes | 14 | 0.79444 | 0.8066 | +0.0121 | +0.14 | 0.3248 | 0.3030 | 0.933 (0.199) |
| BA | plain | 15 | 0.79444 | 0.7821 | -0.0123 | -0.15 | 0.3270 | 0.2400 | 0.734 (0.146) |
| BA | stdbca | 15 | 0.79444 | 0.7821 | -0.0123 | -0.15 | 0.3270 | 0.3612 | 1.105 (0.228) |
| BI | attit | 28 | 0.79444 | 0.7695 | -0.0250 | -0.45 | 0.2947 | 0.3417 | 1.159 (0.172) |
| BI | plainit | 29 | 0.79444 | 0.7636 | -0.0308 | -0.57 | 0.2911 | 0.2045 | 0.702 (0.104) |

## Interval coverage of the true att, nominal 95 (percent, MCSE in brackets)

| cfg | arm | percentile | bias-corrected | BCa | Wald | displayed |
|---|---|---|---|---|---|---|
| AN | bayes | 87.0 (7.0) n=23 | 87.0 (7.0) n=23 | - | 87.0 (7.0) n=23 | 87.0 (7.0) n=23 |
| AN | plain | - | - | - | - | 96.0 (3.9) n=25 |
| AN | std | 87.5 (6.8) n=24 | 87.5 (6.8) n=24 | - | 91.7 (5.6) n=24 | 87.5 (6.8) n=24 |
| AT | bts | 70.0 (14.5) n=10 | 70.0 (14.5) n=10 | - | 70.0 (14.5) n=10 | 100.0 (0.0) n=10 |
| AT | plain | - | - | - | - | 90.0 (9.5) n=10 |
| BA | bayes | 92.9 (6.9) n=14 | 85.7 (9.4) n=14 | - | 92.9 (6.9) n=14 | 92.9 (6.9) n=14 |
| BA | plain | - | - | - | - | 86.7 (8.8) n=15 |
| BA | stdbca | 93.3 (6.4) n=15 | 93.3 (6.4) n=15 | 86.7 (8.8) n=15 | 100.0 (0.0) n=15 | 86.7 (8.8) n=15 |
| BI | attit | 96.4 (3.5) n=28 | 92.9 (4.9) n=28 | - | 100.0 (0.0) n=28 | 96.4 (3.5) n=28 |
| BI | plainit | - | - | - | - | **82.8 (7.0) n=29** |

## Interval coverage of the true att, nominal 90 (percent, MCSE in brackets)

(percentile / BC / BCa taken from the level-independent dual p-values, so this row exists only where the DGP is the null the p tests; Wald recomputed as b +/- z*se.)

| cfg | arm | percentile | bias-corrected | BCa | Wald |
|---|---|---|---|---|---|
| AN | bayes | 87.0 (7.0) n=23 | 87.0 (7.0) n=23 | - | 87.0 (7.0) n=23 |
| AN | plain | - | - | - | 88.0 (6.5) n=25 |
| AN | std | 87.5 (6.8) n=24 | 87.5 (6.8) n=24 | - | 87.5 (6.8) n=24 |
| AT | bts | 70.0 (14.5) n=10 | 70.0 (14.5) n=10 | - | - |
| AT | plain | - | - | - | 70.0 (14.5) n=10 |
| BA | bayes | n/a | n/a | n/a | 85.7 (9.4) n=14 |
| BA | plain | n/a | n/a | n/a | 86.7 (8.8) n=15 |
| BA | stdbca | n/a | n/a | n/a | 93.3 (6.4) n=15 |
| BI | attit | n/a | n/a | n/a | 92.9 (4.9) n=28 |
| BI | plainit | n/a | n/a | n/a | 82.8 (7.0) n=29 |

## Mean interval width by family (same datasets, so widths are directly comparable)

| cfg | arm | n | percentile | bias-corrected | BCa | Wald | displayed |
|---|---|---|---|---|---|---|---|
| AN | bayes | 23 | 1.2410 | 1.2056 | - | 1.2746 | 1.2410 |
| AN | plain | 25 | - | - | - | - | 1.4275 |
| AN | std | 24 | 1.3589 | 1.3355 | - | 1.4285 | 1.3589 |
| AT | bts | 10 | 1.3170 | 1.2846 | - | 1.3752 | 1.7543 |
| AT | plain | 10 | - | - | - | - | 1.3933 |
| BA | bayes | 14 | 1.1849 | 1.1688 | - | 1.3340 | 1.1849 |
| BA | plain | 15 | - | - | - | - | 1.0566 |
| BA | stdbca | 15 | 1.4299 | 1.3506 | 1.3440 | 1.5899 | 1.3440 |
| BI | attit | 28 | 1.3532 | 1.3176 | - | 1.5040 | 1.3532 |
| BI | plainit | 29 | - | - | - | - | 0.9001 |

## Paired standard-error ratios (per dataset, so free of sampling noise in b)

(each bootstrap arm's reported SE divided by the SAME dataset's analytic clustered SE, from the analytic arm run in the same config with the same event-study window.)

| cfg | arm | n | mean se_boot/se_analytic (SE) | min | max |
|---|---|---|---|---|---|
| AN | bayes | 23 | 0.942 (0.042) | 0.668 | 1.531 |
| AN | std | 24 | 1.039 (0.039) | 0.759 | 1.632 |
| AT | bts | 10 | 1.043 (0.091) *coefficient-draw SD from the hidden Wald table, not the displayed studentized se | 0.699 | 1.514 |
| BA | bayes | 14 | 1.276 (0.065) | 0.908 | 1.608 |
| BA | stdbca | 15 | 1.526 (0.102) | 0.932 | 2.274 |
| BI | attit | 28 | 1.723 (0.093) | 1.047 | 2.951 |

## p-value size under the null DGP (rejection percent, MCSE in brackets)

| cfg | arm | p_percentile @5 | @10 | p_bc @5 | @10 | p_bca @5 | @10 |
|---|---|---|---|---|---|---|---|
| AN | bayes | 13.0 (7.0) n=23 | 13.0 (7.0) n=23 | 13.0 (7.0) n=23 | 13.0 (7.0) n=23 | - | - |
| AN | std | 12.5 (6.8) n=24 | 12.5 (6.8) n=24 | 12.5 (6.8) n=24 | 12.5 (6.8) n=24 | - | - |
| AT | bts | 30.0 (14.5) n=10 | 30.0 (14.5) n=10 | 30.0 (14.5) n=10 | 30.0 (14.5) n=10 | - | - |

## Joint pre-trend test (parallel pre-trends hold in both DGPs, so this is size)

| cfg | arm | family | n | rej @5 | rej @10 |
|---|---|---|---|---|---|
| AN | bayes | fixedmetric | 23 | 4.3 (4.3) n=23 | 4.3 (4.3) n=23 |
| AN | std | fixedmetric | 24 | 0.0 (0.0) n=24 | 0.0 (0.0) n=24 |
| AT | bts | studentized | 10 | 0.0 (0.0) n=10 | 0.0 (0.0) n=10 |

## Cells more than 2 MC standard errors from nominal

- BI/plainit main coverage@95 = 82.8 (n=29, MCSE 7.0)
