# fineclus-cellwt — gate archaeology on the refusals T2/T3/T1 removed

## Gate inventory (read from the pinned build, verified line by line)

| old site | old rc | what it refused | status on pinned build |
|---|---|---|---|
| ado:5990 permctl | 459 | permute() with cluster not nesting i | **STILL LIVE** at ado:6078, guarded by `if "`clus'" != "`ivar'" & "`onc'" == ""` — lifted only under `permonsetcohorts` |
| ado:6353 bootctl | 198 | bootstrap(standard) + cell-weighted estimand + NO cluster() | GONE. Replaced by `wdb_reid` set anyway + `wdb_nest=0` keep-i + plain `bsample` |
| ado:6364 bootctl | 459 | bootstrap(standard) + cell-weighted + cluster not nesting i | GONE. Replaced by the two-regime rule at ado:8987-8998 / 9348 |
| ado:6365 bootctl | 459 | bootstrap(standard) + boottstat + cluster() not nesting i | GONE, same replacement |

ado:6451-6456 says so explicitly: "Nothing is screened in its place, deliberately."

## What the replacement actually is (ado:8987-8998, 9341-9350)
```
local wdb_nest 0
if "`reid'" != "" & "`clus'" != "" { ... assert cmin==cmax by(i) ... wdb_nest = (_rc==0) }
...
if "`clus'" == "" bsample
else bsample, cluster(`clusv') idcluster(__wdbns)
if `wdb_nest' qui gegen double __wdbnewi = group(__wdbns `i')
else          qui gen   double __wdbnewi = `i'
```
Consequences I am testing:
1. The predicate is ALL-OR-NOTHING over the whole sample. ONE unit straddling two
   clusters flips the entire run from the multiplicity bootstrap to the presence
   bootstrap. Old behaviour: r459. New behaviour: silent.
2. Under keep-i the equal-cell weighting counts a twice-drawn cluster ONCE, so the
   estimand's cluster multiplicity variance is lost; only 0-vs->=1 PRESENCE survives.
   Asymptotic prediction: the composition part of the variance is scaled by
   (1-p)/p with p = 1-1/e, i.e. 0.582 in variance, 0.763 in SE.
   *** This is exactly the collapse the design note warns about in the OTHER regime. ***
3. `wdb_nest` is posted NOWHERE. `e(perm_irebuilt)` records T1's rebuild; T2/T3's
   regime has no e() record, no display line, nothing.

## Confirmed by reading git (no Stata needed)
- RF-2 (commit 813c414, 2026-09-13) MOVED the reps()-resolution caution from an
  unconditional site (old :6011) to INSIDE `if "`clus'" != "" {` (now :6459-6501).
  So `bootstrap(standard) reps(<40)` with NO cluster() prints no resolution caution.
  T2 made that configuration newly reachable for cell-weighted estimands.

## Runs
- probe1_straddle.do — ONE-CELL regime flip. RUNNING.
- probe2_noclus.do   — unclustered cell-weighted: rc, hand bootstrap, MC SD. QUEUED.
- probe3_overcov.do  — mechanism of the 149/150 over-coverage. QUEUED.

## probe1_straddle.do RESULT (2026-09-17)
10 industries x 3 units, T=6, i=unit, cluster(ind), att+att_i, bootstrap(standard)
reps(99).  Arm B = arm A with ONE CELL changed (unit 1's t=6 row in industry 2).
Arm C = arm A with unit 1 dropped.  nest: A=1, B=0, C=1.  All rc 0.

mean se(att_i) over 2 rseeds:
  dseed 101  A .2268  B .2098  C .2289     B/A = .925
  dseed 102  A .2963  B .2609  C .2935     B/A = .880
  dseed 103  A .2329  B .2381  -           B/A = 1.022
mean se(att)   A .2809  B .2815            B/A = 1.002   <- the control behaves
b(att_i) is BIT-IDENTICAL across A and B (the cluster variable does not enter the
point estimate), so the ONLY thing the one-cell change moves is the SE.

Reading: on a noisy design with only 3 units per cluster the regime flip costs ~6% in
the att_i SE on average and is not consistently signed; the observation-weighted att
SE is untouched, which is the predicted signature (att carries multiplicity through
duplicate ROWS, att_i does not).  probe1c isolates the component.

## FINAL STATE 2026-09-17
All probes complete.  See FINDINGS.md.  Findings: F-A (mechanism, deterministic),
F-B (one-cell regime flip, -23% on se(att_i)), F-F (phantom straddler through a row
the estimator discards, same -23% with e(N) and se(att) bit-identical), F-E (mechanism
of the 149/150 over-coverage, and the record's stated cause has the wrong sign),
F-C (RF-2 regression: reps caution unreachable without cluster()), F-D (no e() record
of the regime), F-G (permute gate scope).  CLEAN: C-1 unclustered cell-weighted SE is
right-sized (pre-registered band), C-2 point estimates never move, C-3 se(att)
untouched, C-4 all three retired gates genuinely run.

Files: probe1_straddle.do/.log/p1.csv, probe1b_mechanism.do/.log,
probe1c_isolate.do/.log/p1c.csv, probe1d_phantom.do/.log/p1d.csv,
probe2a_caution.do/.log, probe2b_noclus_se.do/.log/p2b.csv,
probe3_overcov.do/.log/p3.csv (+p3_partial_reps199.csv),
probe4_permgate.do/.log, nm.do/.log (column-name reconnaissance).
