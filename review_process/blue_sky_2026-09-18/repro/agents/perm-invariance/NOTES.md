# perm-invariance — metamorphic probing of permute() at fixed rseed

## PREDICTIONS STATED BEFORE RUNNING (2026-09-17)

Structural reading of the pinned build (ado 7467-7900, mata 13628-13790, 16510):
- `wooldid_permrun` builds a one-row-per-cluster MAP FRAME from `frame put clus ttre strata if smp`,
  `bysort clus: keep if _n==1`, then `sort permstrata clus`. Cluster keys are unique after the dedup,
  so that sort has NO ties -> the generation order is a deterministic function of the cluster VALUES,
  not of the dataset's row order.
- `wooldid_permx` reads that frame in generation order and, in the sampled branch, calls
  mata `jumble(ons[idx])` once per stratum per rep. RNG consumption depends only on the per-stratum
  vector LENGTH. `set seed rseed` happens first.
- In the EXACT branch (`permexact`) NO random number is consumed at all: `wooldid_permx_enum`
  enumerates every distinct multiset permutation deterministically.

### TIER 1 — the assignment matrix (the `onset_*` block of permsave())
Pure function of (rseed, per-stratum onset multiset, ascending (stratum, cluster-value) order).
MUST be bit-identical under: dataset row shuffle; option-list reorder; the same command run twice in
one session; a fresh session; relabeling of `i`; any transformation of `y`; an affine remap of
(t, ttre); duplicating every row; `if` vs physically dropping; weight rescaling.
MAY legitimately differ under: a relabeling of CLUSTER ids that changes their ascending order
(different clusters then sit at different generation rows, so a fixed seed draws a different map).

### TIER 2 — the null distribution e(permdist_*) / the accumulator
Identical wherever tier 1 is identical AND the estimator is equivariant.
Under y -> a + b*y (b>0) every permuted coefficient scales by b and the intercept shift drops out.

### TIER 3 — reported p-values and intervals
Identical wherever tier 2 is identical up to a common positive scale.
TWO EXTRA invariances that must hold even where tier 1 moves:
 (a) under `permexact` the enumerated column SET is every distinct assignment, so the null
     distribution is a multiset invariant: p MUST be bit-identical under ANY cluster relabeling
     and any row order.
 (b) the default equal-tailed p is 2*min(p_lo,p_hi): it must be invariant under y -> -y
     (which swaps p_lo and p_hi).

A failure in tier 1 but not tier 3 = latent reproducibility bug. A failure in tier 3 = worse.

## STATUS
- [x] read brief + known_bydesign
- [x] read permute engine source
- [ ] probes

## RESULTS (running log)

### probeA — order and label transformations, sampled draws (60 reps, 40 units / 10 states)
All nine variants rc=0.  Tier-1 (permsave assignment block, values only):
 rerun, rowsh, revrow, optord, irelU, irelD, crelU : IDENTICAL
 crelD (order-REVERSING cluster relabel)           : DIFFERS  <- PREDICTED, legitimate
Tier-2/3 (every e() scalar and every perm_* matrix, rtol 1e-10): ALL EQUAL for
 rowsh, revrow, optord (bitwise), irelU, irelD, crelU.  crelD differs (MC noise).
=> PREDICTION CONFIRMED EXACTLY.  No order- or label-dependence.

### independent p-value oracle (python, from the dumped e(permdist_b)/e(perm_comparator_b))
Sampled run: reported p reproduces 2*min((#ge+1)/(n+1),(#le+1)/(n+1)) for all 14
coefficients to 1e-12.  permexact run: reported p reproduces the PLAIN share
2*min(#ge,#le)/n (n=90), i.e. the observed assignment is inside its own reference set
and is NOT double counted.  Shift interval verified as [b_obs - qhi, b_obs - qlo].

### probeC — permexact (90 assignments)
e(perm_reps)=90 = 6!/(2!2!2!) as predicted.
C1rowsh, C2crelU, C4irelD : ALL EQUAL (rtol 1e-10) including every p.
C3crelD (order-reversing cluster relabel): ONLY e(permdist)/e(permdist_b) differ --
the ROW ORDER of the stored enumeration moves -- every p, interval, band and kernel
entry is bit-identical.  Exactly the predicted multiset invariance.

### RNG hygiene
 - documented replay `set rngstate e(perm_rngstate)` then repeat the call: BITWISE
   identical (H3a vs H3b).  The documented claim at sthlp:542 holds.
 - user `set seed` before a no-rseed run is honoured and reproducible (H4a==H4b).
 - a sampled permute() run MOVES c(rngstate) (documented).
 - after a run with rseed(), a subsequent runiform() is identical regardless of the
   user's own pre-run seed: rseed() fully determines the post-run stream (standard
   `set seed`-inside-a-command behaviour, as Stata's own prefixes).

### probeB — data transformations (60 reps)
 B1dup  (every row duplicated): assignment IDENTICAL; every coefficient, p-value,
        interval and band IDENTICAL.  Only N / escomp_* / ll / rss / rmse / r2 move.
 B2if vs B3drop (`if state<=9` vs physically dropping those rows): ALL EQUAL (1e-10),
        assignment identical to each other.  BS-37's claim holds for the permute design.
 B4tsh  ((t,ttre) + 50): ALL EQUAL (1e-10); the assignment block differs only by the
        +50 carried in the onset VALUES.
 B5yadd (y + 1000): ALL EQUAL (1e-10).
 B6ymul (y x 1000): every p-value and t-ratio BIT-IDENTICAL; sd/lb/ub/band scale by
        exactly 1000 (rel 1e-16).
 B7yneg (y -> -y): every p-value BIT-IDENTICAL; the null band mirrors correctly
        (bandlo(-y) = -bandhi(y)).  The equal-tailed convention is sign-symmetric.
 B8w vs B9w1k (aweights x 1000): every reported number equal within 1e-10; the
        replication accumulator e(permdist_b) moves by 4.6e-10 relative at one cell,
        which is reghdfe's tol(1e-9) convergence floor, not a defect.  e(escomp_refW)
        and e(escomp_trtW) report the RAW weight sum, so they scale by 1000.

### probeD — fidelity and further invariances (40 reps)
 FIDELITY: for reps 1, 7, 23 and 40, the assignment read out of ps_D0.dta's onset_*
 columns, written back onto the panel and re-estimated with the engine's own rep option
 string (att att_i fastnose cluster(state) esprelength(2) espostlength(3)
 esrelativeto(-1)) reproduces the stored replication coefficients exactly or to 1 ulp.
 => the permsave() assignment block is correctly labelled by cluster.
 permreps PREFIX: the first 40 assignment rows of a permreps(80) run are bit-identical
 to a permreps(40) run at the same seed -- RNG consumption per rep does not depend on R.
 D2attonly: dropping att_i from the estimand list leaves the att p bit-identical.
 D3str1:    permstrata(a constant variable) is bit-identical to no permstrata.
 D4noes:    dropping the event-study window leaves the att p bit-identical.

### probeH / probeE2-J — THE ONE FINDING (PI-1)
A cluster all of whose rows have a missing OUTCOME stays in the permutation design.
 H1ymiss vs H2ydrop: identical estimation sample (N=360, G=9) and bit-identical point
 estimates, but 10 vs 9 clusters in the design, e(perm_npossible) 16800 vs 5040, and
 different p-values and intervals.  H3cmiss (missing CLUSTER value) == H2ydrop, and
 J2 (`if state != 7`) == J1 (dropped): so missing-cluster and `if` remove the cluster,
 missing-outcome does not.  H4partial == H4bdrop: partial missingness is neutral.
 J0: with the phantom the SOLE carrier of onset 2008, a data-carrying cluster received
 the 2008 onset in 29 of 30 replications although the observed estimation sample has no
 2008 cohort at all.  perm_failed=0, perm_blank=0 -- the mechanism is a different design,
 not broken replications.
 SCOPE: wooldid_bootrun uses the same predicate (ado:8912-8923); the repro shows
 bootstrap e(V)[1,1] = 0.01910091 (rows kept) vs 0.01501891 (rows dropped).

### probeE — permonsetcohorts (shipped 2026-09-16)
 E0/E1irelD/E2crelU/E3rowsh: reversing the i codes, relabelling the clusters and
 shuffling the rows all leave every number equal to 1e-10.  e(perm_irebuilt)=1.
 E4onc vs E5noonc: with cluster()==i the option is a documented no-op -- the two dumps
 differ ONLY by the presence of the scalar e(perm_irebuilt)=0.  CLAIM CONFIRMED.

### probeE2 — interleaved strata under permexact (the NFPERM-02 geometry)
 8 clusters, strata {1,3,5,7} and {2,4,6,8}, so ascending cluster order crosses stratum
 order.  e(perm_reps)=36 = 6x6 as predicted.  Reversing the STRATUM labels, and reversing
 the CLUSTER labels, each leave every p, interval, band and kernel entry bit-identical;
 only the row order of e(permdist)/e(permdist_b) moves.  The shipped NFPERM-02 fix holds.

### probeI — permute(probabilities) fidelity
 A table with eight deterministic rows pointing each cluster at an onset it does NOT have,
 written (a) ascending, (b) descending, (c) with the cluster ids relabelled 100-state:
 in every case every cluster drew exactly the onset its own row specifies, across all 24
 replications.  The permprobs()/onset pairing is value-keyed, not positional-by-listing.

### fresh-session determinism
 dump_base.txt (probeA session) vs dump_B0.txt (probeB session): identical at rtol 0.

### NOT a finding, recorded to avoid re-derivation
 - HYG_EXACT_STATE_UNCHANGED in probeC was CONFOUNDED (the previous run left the stream at
   the same rseed state).  It is NOT evidence that permexact leaves c(rngstate) untouched.
 - H5future: never-treated coded as an out-of-window onset (2020 on a panel ending 2009)
   gives results bit-identical to plain missing, design included.  Clean.

### probeF — conventions under outcome transformations
 permtstat with y, 1000y, 1e-6 y: every e(perm_p) bit-identical, t differs 2e-14.
 permabsp / permside / permstepdown / default under y -> -y: reported p bit-identical,
 p_lo and p_hi swap (as they must), step-down adjusted p's identical.
 permside(left) on y == permside(right) on -y.  Sanity: permside(left) and permabsp DO
 move the p relative to the default, so the passes are not vacuous.
 simultaneousbands under permute() is refused with an informative r(198).

### probeG — subgroups() label swap
 sg 1<->2: G0's wdidsgresults_1 == G1's wdidsgresults_2 bit for bit, and conversely;
 the pairwise contrast flips sign, mirrors its interval, keeps its p.  Order-preserving
 relabel to {5,9} renames the keys only.  No surface crossed the subgroups.

## STATUS: COMPLETE.  One finding (PI-1), everything else clean.
