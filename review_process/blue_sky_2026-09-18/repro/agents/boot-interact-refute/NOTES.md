# boot-interact-refute — running notes (SKEPTIC)

Build asserted: r00_setup.log line 54 -> .../scratchpad/bluesky/build/wooldid.ado, "version 1.4 2026-08-19".

## My own data (all probes)
30 states x 2 counties = 60 counties, 16 periods, 10 treated states with onsets 7..12
(one every third state), y = 2 + .05*t + unit effect + treatment effect + N(0,1).
**No pre-trend in the DGP.** i = county, cluster(state), G = 30.
(r00 used a 4-county/18-period variant; same conclusions, slower.)

## r01_mech.do / r01_mech.log — BOOTX-1 mechanism: CONFIRMED exactly
- p_printed == p_hand = Ftail(keff, e(df_r), W/keff) with W = b'S^-1 b on the
  BOOTSTRAP pre block, at EVERY arm (20/50/100/300 reps and bootdf clusters/normal/3).
  Identity holds to every printed digit. e(df_r) = 29 = G-1 at every reps().
- BUT on well-behaved data at realistic reps the numbers are undramatic:
  analytic p .000676 (keff 5/6, rank caution) | boot20 .5848 | boot50 .7775 |
  boot100 .7393 | boot300 .6206; Hotelling .7389/.8196/.7646/.6264; wooldid's own
  rank joint test .7143/.7843/.7129/.5800. Wald and rank AGREE at realistic reps.
- bootdf moves it, as claimed: .77752 (clusters) / .78247 (normal, chi2) / .764631 (3)
  at identical W = 3.20688. Magnitude tiny in the body of the distribution.
- **ANALYTIC keff = 5 of 6 with the rank caution; EVERY bootstrap arm keff = 6 of 6
  with no caution — at reps(300) on an ordinary 6-pre-period spec.** This is the
  realistic face of BOOTX-2 and it reproduces on my own plain panel.

## r03_oracle.py / r03_oracle.out — size + coverage oracle (no Stata)
b ~ N(0,I), S ~ Wishart(I,n-1)/(n-1) independent; wooldid's p = Ftail(k, dfr, W/k).
size at nominal 5% (k, n, dfr=G-1):
 (6,20,29) .175 | (6,50,29) .0628 | (6,100,29) .0395 | (6,200,29) .0312 |
 (6,999,29) .0248 | (11,50,29) .1087 | (11,60,39) .1012 | (8,50,29) .0782 |
 (6,50,199) .0963 | (16,200,29) .0270.  Hotelling is .049-.051 EVERYWHERE.
=> the distortion is NOT monotone as the prober claims: anti-conservative below
   n ~ 70, CONSERVATIVE above (half nominal size at reps(999) with 30 clusters).
   More clusters at fixed reps makes the anti-conservative side WORSE (.0963 at
   dfr=199 vs .0628 at dfr=29, same n=50, k=6).
Roth slope CI coverage at nominal 95%: (6,50) .928 | (11,50) .887 | (6,20) .843 |
   (6,200) .954 | (11,60,dfr39) .901.
mdt_wald understatement factor sqrt((n-k-2)/(n-1)): .926 at (6,50), .869 at (11,50).
   NOT the 3.8x the prober reports — that figure is a reps(7) artifact.

## r05/r07/r06/r08/r09/r10b/r11/r12b — see FINDINGS.md for the full write-up
Harness trap worth recording: a `mata: ... end` block nested inside a Stata `foreach { }` loop
or inside a `prog def` aborts the do-file (--Break-- r(1), or "command end is unrecognized"
r(199)). Keep mata blocks at top level. Second trap: a Stata launched from a bash call that
RETURNS is killed with --Break-- when the harness reaps the process group; keep the calling
shell alive until Stata exits.

## Final state (all probes complete)
FINDINGS.md holds the five verdicts (BOOTX-1 PARTIAL/medium, BOOTX-2 PARTIAL/medium with the
headline mechanism REFUTED, BOOTX-3 PARTIAL/low, BOOTX-4 CONFIRMED/low, BOOTX-5 PARTIAL/low)
and three new findings (BX-N1 medium/certain, BX-N2 medium/probable, BX-N3 medium/probable).
Decisive probes: r07_nulldir (refutes BOOTX-2's mechanism), r13b_disc (part A = BX-N3
discriminator, part B = BOOTX-4 reproduction), r11_boottstat (BX-N1), r03_oracle + r03_limits
(BOOTX-1 size).
