# SKEPTIC NOTES — refuting PI-1 (phantom cluster in the permute() design)

## THE LOAD-BEARING FACT (found before any Stata run)

`wooldid.sthlp` lines 430-431 of the pinned build document PI-1's mechanism VERBATIM,
inside the `permute()` option entry:

> "Requires cluster() with a single variable - the randomization unit (**the randomization
> design consists of the clusters having at least one row where the *if* expression
> evaluates to exactly 1 with a nonmissing cluster value; rows *wooldid* later removes by
> its deterministic markout, e.g. missing outcomes or zero weights, still count toward the
> design and its possible-assignment count**) - and every cluster must sit inside exactly
> one treatment-onset group ..."

So:
 - the design population predicate (if==1 & !missing(cluster)) is DOCUMENTED,
 - "rows removed by markout still count" is DOCUMENTED,
 - **missing outcomes are named as the example**,
 - "and its possible-assignment count" documents the e(perm_npossible) difference,
 - the three exclusion routes' asymmetry is exactly what this sentence defines.

PI-1's EXPECTED rests on "sthlp:61 states only that the `if` rule applies identically ...
nothing says a missing outcome is treated differently from a missing cluster id."
That is FALSE. sthlp:430-431 says exactly that. PI-1 cited sthlp:61 and missed 430.

## PLAN
E1  permexact (RNG-free) A/B/C/D/E variants — does the reference SET differ, not just the
    RNG stream? (the prober's 30-rep p-differences are confounded with the shuffle stream)
E4  NEW: does the permctl NESTING check (run on `insmp`, which includes markout-dead rows)
    REFUSE r(459) on a panel that estimates cleanly without permute()?
E5  NEW: measure the bootstrap extension PI-1 asserted from code alone. Theory:
    bootstrap(bayes) is EXACTLY invariant in distribution (normalized Dirichlet sub-vector
    over the real clusters is Dirichlet over their own alphas; aweights are scale-free),
    so the shared predicate has NO consequence there.

## STATUS
- [x] read brief, known_bydesign, ado 7481-7490 / 5959-5964 / 8912-8923, sthlp 430-431
- [x] E1  - [x] E4  - [x] E5  - [x] E6  - [x] E6B

## E1 RESULT (e1_exact.do / e1_exact.log) — RNG-FREE, permexact

8 states x 3 counties x 8 years; state 1 = 2004, state 2 = 2006, rest never;
state 8 = the phantom.  `att cluster(state) permute(shuffle) permexact`.
Every arm: e(N)=168, e(N_clust)=7, e(perm_failed)=0.

| arm                                   | npossible | exact p   | 95% shift CI        | b (hex)              |
|---------------------------------------|-----------|-----------|---------------------|----------------------|
| A  phantom, y missing on all its rows |    56     | .07142857 | [-.5312, .8277]     | +1.a7c7591ced275X-001|
| B  phantom rows dropped               |    42     | .04761905 | [-.5130, .5876]     | +1.a7c7591ced275X-001|
| C  phantom cluster value set missing  |    42     | .04761905 | [-.5130, .5876]     | same                 |
| D  `if state != 8`                    |    42     | .04761905 | [-.5130, .5876]     | same                 |
| E  phantom killed by missing CONTROL  |    56     | .07142857 | [-.4588, .9437]     | +1.e331bfc86f07cX-001|
| E' same, phantom rows dropped         |    42     | .04761905 | [-.4571, .5168]     | same b               |
| J  phantom is the SOLE 2008 adopter   |   336     | .07142857 | [-.5312, .8277]     | same as A            |
| J' same, phantom rows dropped         |    42     | .04761905 | [-.5130, .5876]     | same                 |

So PI-1's NUMBERS reproduce, and harder than PI-1 showed them: with exact
enumeration there is no RNG confound at all (PI-1's 30-rep contrast was
confounded — 10 clusters shuffle differently than 9).  The p CROSSES .05
(.0714 vs .0476).  B == C == D bit-for-bit: the three exclusion routes really
do split 1-vs-2, exactly as PI-1 says.  Arm E shows the same holds when the
killing variable is a CONTROL, not the outcome — the help's "e.g." covers it.

BUT: every one of these is the behaviour sthlp:430-431 SPECIFIES, naming
missing outcomes as its example and naming "and its possible-assignment count"
as the thing that moves.  So the claim's "silently / with no caution /
nothing says a missing outcome is treated differently from a missing cluster
id" is FALSE, and the EXPECTED it is built on collapses.
=> VERDICT PARTIAL, severity low (documented, deliberate design; the numbers
   are right, the characterisation as an undocumented defect is not).

## E4 RESULT (e4_nesting.log / e4b_nesting.log)  — NEW, mine
State 3 adopts 2006; its county 9 is a merge miss (ttre==. AND y==. on all 8 rows).
Onset constant within every i, so ado:772's within-i screen does not fire.
 V1d: keep vs drop -> e(N)=184 BOTH, att +1.8da894f951c6cX-001 BOTH, identical=1
 V1b: keep + permute(shuffle) permexact -> r(459) "clusters ... mix onset values
      or mix treated and never-treated rows" naming 24 in-sample observations
      (the whole cluster; only 8 of them are the dead ones)
 V1c: drop + permute -> runs, npossible=168
 V2b: dead county carrying a stale onset 2010 instead -> r(459) as well
=> a runnable specification becomes a hard error because of rows no fit can use.

## E6 RESULT (e6_bootscreens.log) — NEW, mine
bootstrap()'s cluster-count screens (ado:6460-6500) count over
smp2 = ((if)==1) & !missing(cluster) -- the phantom-inclusive population.
 A1 9 real + 1 phantom : N=216 Nclust=9, "only 9 clusters are available for the
    bootstrap" caution NOT printed
 A2 phantom dropped    : N=216 Nclust=9, caution IS printed
 B1 5 real treated + 1 phantom treated : N=264 Nclust=11, "only N treated
    cluster(s) are available for resampling" caution NOT printed
 B2 phantom dropped    : N=264 Nclust=11, "only 5 treated cluster(s)" printed
 C1 1 real + 1 phantom : the two-cluster screen (whose message says "in the
    estimation sample") is PASSED, and the cautions print the inflated count
    ("only 2 clusters"), but a downstream guard catches it:
    "Error: cluster(state) leaves 1 cluster(s) in the estimation sample" r(459).
    => no live consequence for (c); report (a)+(b) only.

## E5 RESULT (e5_boot.log) — NEGATIVE, refutes PI-1's bootstrap extension
9 clusters (4 treated) + 1 phantom vs phantom dropped, att cluster(state) reps(400),
3 seeds each. mean bootstrap SE:
  bayes    kept .072886  dropped .072812  ratio 1.0010
  standard kept .089842  dropped .090492  ratio 0.9928
Both inside MC noise. Structural reason: aweights are scale-free and the
normalized weight sub-vector over the real clusters has the same distribution
with or without the phantom (Dirichlet sub-vector / multinomial reweighting).
=> bootstrap() shares the PREDICATE but not the consequence.

## E6B RESULT (e6b_controls.log)
D1 phantom treated killed by a missing CONTROL: N=264 Nclust=11, few-treated
   caution NOT printed.  D2 phantom dropped: same N/Nclust, "only 5 treated
   cluster(s) are available for resampling" printed.

## STATUS: COMPLETE.  FINDINGS.md written.

## E7 RESULT (e7_physleg.log) — the PHYSICAL-resample legs, also invariant
att_i under bootstrap(standard) (physical resample + re-identification),
3 seeds: mean SE kept .082352 vs dropped .082811, ratio 0.9945.
boottstat att reps(200), 2 seeds: p .050/.040 kept vs .040/.050 dropped;
intervals A [.245,1.182] [.310,1.217] vs B [.378,1.123] [.278,1.211] -
within-arm seed spread >= between-arm spread, so nothing structural.
