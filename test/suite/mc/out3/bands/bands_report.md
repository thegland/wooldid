# Simultaneous sup-t band coverage (mc_bands, PS-45)

outdir `out3/bands`; 6 shard file(s); 450 rows after de-duplication on (arm, tau, shard, rep); 0 malformed line(s) dropped; 0 duplicate row(s) dropped.

The band's claim is SIMULTANEOUS: `band_all` = 1 iff EVERY reported horizon's truth is inside the band on that rep. Nominal 0.95 (mc_bands.do passes no level(), so wooldid's default level(95) is in force). Pointwise coverage is the CONTROL: each horizon near nominal, and the joint pointwise event `pw_all` BELOW `band_all`.

Judgment band = central ~95% binomial region, normal approx rounded outward (house convention, analyze_mc.py); the interval printed beside every rate is a Wilson 95% score interval. Cells with fewer than MINN=30 reps are reported, never judged.

| check | verdict | detail |
|---|---|---|
| rc == 0 on every row | PASS | rc!=0 on 0/450; threshold 0.00 (driver ACCEPTANCE: rc = 0 on every row) |

## Run inventory

| arm | engine(s) | rows (rc=0) | rc!=0 | band_absent | mean secs | mean crit | ncomp (modal) | nfail>0 |
|---|---|---|---|---|---|---|---|---|
| an | analytic-simulated | 150 | 0 | 0 | 2.28 | 2.6438 | - | 0 |
| bs | bootstrap | 150 | 0 | 0 | 209.72 | 2.6584 | 199 | 0 |
| bt | bootstrap-studentized | 150 | 0 | 0 | 433.94 | 5.0493 | 199 | 0 |

## 0. Integrity (per-row facts; judged at any n - these are not statistical)

| check | verdict | detail |
|---|---|---|
| an band present on every rep | PASS | band_absent>0 on 0/150 (a rep with no e(esband_att) means the kernel refused the whole family on one degenerate horizon, so that rep contributes NO band event; driver ACCEPTANCE: band_absent = 0) |
| an nh == 8 | PASS | nh != 8 on 0/150 (observed [8.0]); the window is esprelength(4) espostlength(4) esfixedbaseperiod |
| an pw_n == nh | PASS | pw_n != nh on 0/150 (the pointwise control would then be measured on a different horizon set from the band) |
| an band_hits == sum(bd*) | PASS | mismatch on 0/150 |
| an pw_hits == sum(pw*) | PASS | mismatch on 0/150 |
| an band_all == 1[band_hits == nh] | PASS | mismatch on 0/150 (band_all IS the simultaneous event) |
| an row count == M (150) | PASS | 150 rc-0 rows of 150 planned |
| bs band present on every rep | PASS | band_absent>0 on 0/150 (a rep with no e(esband_att) means the kernel refused the whole family on one degenerate horizon, so that rep contributes NO band event; driver ACCEPTANCE: band_absent = 0) |
| bs nh == 8 | PASS | nh != 8 on 0/150 (observed [8.0]); the window is esprelength(4) espostlength(4) esfixedbaseperiod |
| bs pw_n == nh | PASS | pw_n != nh on 0/150 (the pointwise control would then be measured on a different horizon set from the band) |
| bs band_hits == sum(bd*) | PASS | mismatch on 0/150 |
| bs pw_hits == sum(pw*) | PASS | mismatch on 0/150 |
| bs band_all == 1[band_hits == nh] | PASS | mismatch on 0/150 (band_all IS the simultaneous event) |
| bs ncomp == B (199) | PASS | ncomp != 199 on 0/150 (a lost replication means the critical value was simulated from fewer draws than requested) |
| bs nfail == 0 | PASS | nfail>0 on 0/150 |
| bs row count == M (150) | PASS | 150 rc-0 rows of 150 planned |
| bt band present on every rep | PASS | band_absent>0 on 0/150 (a rep with no e(esband_att) means the kernel refused the whole family on one degenerate horizon, so that rep contributes NO band event; driver ACCEPTANCE: band_absent = 0) |
| bt nh == 8 | PASS | nh != 8 on 0/150 (observed [8.0]); the window is esprelength(4) espostlength(4) esfixedbaseperiod |
| bt pw_n == nh | PASS | pw_n != nh on 0/150 (the pointwise control would then be measured on a different horizon set from the band) |
| bt band_hits == sum(bd*) | PASS | mismatch on 0/150 |
| bt pw_hits == sum(pw*) | PASS | mismatch on 0/150 |
| bt band_all == 1[band_hits == nh] | PASS | mismatch on 0/150 (band_all IS the simultaneous event) |
| bt ncomp == B (199) | PASS | ncomp != 199 on 0/150 (a lost replication means the critical value was simulated from fewer draws than requested) |
| bt nfail == 0 | PASS | nfail>0 on 0/150 |
| bt row count == M (150) | PASS | 150 rc-0 rows of 150 planned |

## 1. SIMULTANEOUS coverage - the headline (nominal 0.95)

| arm | engine | band_all rate (n, Wilson 95%) | count band |
|---|---|---|---|
| an | analytic-simulated | 139/150 = 0.9267 [0.8735,0.9586] | [137,148] |
| bs | bootstrap | 139/150 = 0.9267 [0.8735,0.9586] | [137,148] |
| bt | bootstrap-studentized | 146/150 = 0.9733 [0.9334,0.9896] | [137,148] |

| check | verdict | detail |
|---|---|---|
| an (analytic-simulated) simultaneous coverage | PASS | band_all 139/150 = 0.9267 [0.8735,0.9586] vs nominal 0.95, count band [137,148], MARGINAL corridor [0.90,0.98] |
| bs (bootstrap) simultaneous coverage | PASS | band_all 139/150 = 0.9267 [0.8735,0.9586] vs nominal 0.95, count band [137,148], MARGINAL corridor [0.90,0.98] |
| bt (bootstrap-studentized) simultaneous coverage | PASS | band_all 146/150 = 0.9733 [0.9334,0.9896] vs nominal 0.95, count band [137,148], MARGINAL corridor [0.90,0.98] |

## 2. POINTWISE coverage per horizon - the control (nominal 0.95 each)

| arm | rel -5 | rel -4 | rel -3 | rel -2 | rel 0 | rel 1 | rel 2 | rel 3 | n |
|---|---|---|---|---|---|---|---|---|---|
| an | 141/150 | 139/150 | 141/150 | 142/150 | 138/150 | 144/150 | 146/150 | 139/150 | 150 |
| bs | 141/150 | 141/150 | 139/150 | 141/150 | 137/150 | 142/150 | 145/150 | 140/150 | 150 |
| bt | 145/150 | 140/150 | 141/150 | 141/150 | 141/150 | 143/150 | 144/150 | 141/150 | 150 |

(* = outside the central ~95% binomial region for that horizon.)

| check | verdict | detail |
|---|---|---|
| an per-horizon pointwise coverage (8 horizons) | PASS | all horizons inside their central ~95% binomial region |
| bs per-horizon pointwise coverage (8 horizons) | PASS | all horizons inside their central ~95% binomial region |
| bt per-horizon pointwise coverage (8 horizons) | PASS | all horizons inside their central ~95% binomial region |

## 3. The multiplicity gap - joint POINTWISE coverage must sit BELOW the band's

(With 8 horizons, intervals that only claim pointwise coverage cover ALL of them well under 95% unless the horizons are nearly perfectly correlated. If pw_all is not below band_all the band is buying nothing.)

| arm | pw_all rate (n, Wilson) | band_all rate (n, Wilson) | gap |
|---|---|---|---|
| an | 105/150 = 0.7000 [0.6224,0.7676] | 139/150 = 0.9267 [0.8735,0.9586] | +0.2267 |
| bs | 106/150 = 0.7067 [0.6294,0.7736] | 139/150 = 0.9267 [0.8735,0.9586] | +0.2200 |
| bt | 107/150 = 0.7133 [0.6364,0.7797] | 146/150 = 0.9733 [0.9334,0.9896] | +0.2600 |

| check | verdict | detail |
|---|---|---|
| an pw_all share below band_all share | PASS | pw_all 105/150 = 0.7000 [0.6224,0.7676], band_all 139/150 = 0.9267 [0.8735,0.9586], gap +0.2267 (directional check; no band - the sign is the claim) |
| bs pw_all share below band_all share | PASS | pw_all 106/150 = 0.7067 [0.6294,0.7736], band_all 139/150 = 0.9267 [0.8735,0.9586], gap +0.2200 (directional check; no band - the sign is the claim) |
| bt pw_all share below band_all share | PASS | pw_all 107/150 = 0.7133 [0.6364,0.7797], band_all 146/150 = 0.9733 [0.9334,0.9896], gap +0.2600 (directional check; no band - the sign is the claim) |

## 4. Secondary (information only, never fails the component)

| arm | mean crit | sd crit | mean band_hits/nh | mean pw_hits/pw_n | mean secs |
|---|---|---|---|---|---|
| an | 2.6438 | 0.0161 | 0.9900 | 0.9417 | 2.28 |
| bs | 2.6584 | 0.0861 | 0.9875 | 0.9383 | 209.72 |
| bt | 5.0493 | 0.6479 | 0.9967 | 0.9467 | 433.94 |

The critical value is the band's whole content: `crit` well above the pointwise 1.96 (or the engine's t) is what widens the band, and a `crit` that does NOT move with the engine would mean the simultaneity is nominal only.


## Verdict summary

- checks FAILED: 0
- checks not judged for want of reps (MINN=30): 0
