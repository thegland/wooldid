# perm-exact probe notes  (agent: permexact vs independent brute-force enumeration)

Build asserted: `which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado, "version 1.4 2026-08-19".

## Harness facts learned (cost me time; record for a successor)
1. wooldid's Mata functions are NOT visible in the session after an ordinary `wooldid` call
   (`wooldid_permx() not found`, r(3499)). To drive the enumerator directly I extracted the
   three functions VERBATIM:  `sed -n '13624,13785p;13790,13840p' build/wooldid.ado >
   permx_src.mata.do`  (md5 e5cdfc392613bb2c1b11a61a077cf7c2) and `do` that.
2. `wooldid_permx`'s `pname`/`lname` must be NON-EMPTY names (Mata `&` has no short-circuit,
   so `st_matrix("")` errors r(3300)). The real caller always passes tempnames.

## p02_enum.do  -- M (log-factorial count) vs A (recursive enumerator).  SEEDED LEAD: CLEAN.
20 designs. For EVERY plain-missing design: M == cols(A) == distinct(A), every column
preserves the per-stratum multiset, and the observed assignment appears EXACTLY ONCE.
Ties, singleton strata, an all-never-treated stratum, uneven strata, no-never-treated,
negative onsets, fractional onsets, 3 strata, 12-cluster heavy cancellation, 7560 columns:
all reconcile. Boundary: M=18900 and M=40320 both refused with the documented r(198) text
and ceiling 10000 -- NO fallback to sampling. Matches the help (caution >1,000, refuse >10,000).
=> The two independent algorithms do NOT disagree on any reachable design.

EXCEPTION (matches the code's own comment at ado:13695-13707): extended missings.
  ons = 4 7 .a . . .  -> M=1   but cols(A)=240 (30 distinct)
  ons = 4 .a .b . . . -> M=0   but cols(A)=324 (6 distinct)
Both are refused at parse (ado:5963-5969, r459), so this is unreachable from the command
line -- TO BE VERIFIED empirically, incl. under safety(off).

OPEN: permstrata() with MISSING values silently makes "missing" its own stratum
  (c11: str = 1 1 . . 2 2 -> M=2). Check whether the parser guards that.

## Still to run
p01_conv.do (6 conventions, dumps for the python oracle), p03_oracle.do (independent
brute-force refit of all 30 assignments), p05 guards/boundary, permonsetcohorts.

---

## p01_conv.do + p01_oracle.py -- every documented convention, checked against an oracle
6-cluster design, M = 30, 30/30 completed, `e(perm_addone)` and `e(perm_incomplete)` empty.
All six runs (default / permabsp / permrawp / permside(left) / permside(right) / permtstat)
reproduce EXACTLY (60 checks, 0 disagreements) from the draws alone:
  default   p = min(1, 2 min(m0,m1)/n)          = 2/30 = .0667
  permabsp  p = ma/n                            = 1/30 = .0333
  permrawp  p = ma/n                            = 1/30 = .0333   (identical to permabsp
            under an enumeration -- absp==1 takes ma/n when exact|memb, absp==2 always)
  left      p = m0/n = 1.0 ;  right  p = m1/n = .0333
  permtstat same formulas on b_r/se_r
plus e(perm_p_lo/p_hi), the population SD convention under `exact`, e(perm_nullband),
e(perm_ci) = b_obs - [qhi, qlo] under wooldid_permq's rule, and e(perm_p_se) left MISSING
under the enumeration (as sthlp:512 promises).

## p03_oracle.do -- independent brute force, BITWISE
Enumerated all 30 assignments myself (cluster a gets onset 4, cluster b != a gets 7) and
refit on the engine's own replication path. Sorted draws are BITWISE identical to
e(permdist_b) for the coefficient path, and to e(permdist_b)/e(permdist_se) for permtstat.
e(perm_comparator_b) is bitwise the identity assignment's own fastnose fit.

## p04_rich.do + p04b_brute.do -- 54 permuted elements, incl. INCOMPLETE ones
8 clusters, M = 56, esprelength(2) espostlength(3), 4 estimands, subgroups(), jointtests.
54 elements; 18 complete (n = 56) and 36 estimable in only n = 44 assignments.
270 checks, 0 disagreements; every one of 54 x 56 draws bitwise identical to my refits.
The incomplete branch is right: exk = 0, memb = 1, exact share over n = 44, floor 2/44,
e(perm_n) = 44.

## p05_strata.do -- REGRESSION test of NFPERM-02 (crossing strata). CLEAN.
Generation order (strata, cluster) = 5 6 7 8 1 2 3 4 while frlink re-sorts to 1..8.
M = 16, reconciles bitwise; every saved assignment respects its OWN stratum; the observed
assignment is found (memb = 1) and p = 2/16. The __wdprow fix holds.

## p06_onc.do -- permexact + permonsetcohorts (T1, shipped one day before this review)
i = cohort, cluster = unit (strictly finer), M = 30, e(perm_irebuilt) = 1.
Reproduced all 30 replications by hand (rebuild i as "the units sharing the drawn onset",
observed cohort codes) -- BITWISE identical, p exact, no add-one fallback.

## p07/p07b_edges.do -- guards and boundary
E1 extended missing onset -> r(459);  E2 same under safety(off) -> r(459) (NOT bypassable)
E3 onset varying within cluster (i = unit, cluster = state) -> r(459), and under
   safety(off) too. My "silent onset collapse" hypothesis is GUARDED (ado:5971-5982).
E4 permstrata() with missing values -> r(459);  E5 permstrata() not constant within
   cluster -> r(459) (ado:6138-6149)
E6 M = 18900 -> r(198) with the ceiling text; NO fallback to sampling, e(perm_reps) unposted
E7 permexact+permreps() -> r(198);  E8 permexact+permute(probabilities) -> r(198)
E9 degenerate M = 1 -> r(459)
E10/E11 customvce(cluster ...) and E12 unconditionalse -> REFUSED r(198) with permute().
   This closes the only route I could see to a bitwise mismatch between the rep-0
   comparator and the identity replication (the brief names those two paths as the
   non-bit-reproducible ones).
E13 poisson, E15 aweights, E16 permstepdown -> run, e(perm_addone) empty
E14 a cluster whose y is entirely missing still counts toward M (NPOSS = 30 on 6 clusters
   with one wholly dropped) -- exactly as sthlp:431 states.

## p08_more.do
F1 cluster(state t) (two-way) is accepted by the ORDINARY path (rc 0);
F2/F3 permute()+permexact / permreps -> informative r(198) "permute() requires cluster()
   with a single variable"; F3b bootstrap() likewise r(198). Informative, not a defect.
F4 `if unit<=6` restricts the DESIGN: NPOSS 56 -> 30. Correct.
F7 permtci verified against the oracle: q = wooldid_permq(sort(|t_r|), level/100) = the
   29th of 30 order statistics; reported interval matches to display precision.

## p09_exactness.py -- EXACT SIZE, no simulation
For each of the 56 assignments taken in turn as the observed one, applied wooldid's own
exact-share rule and counted rejections. For all 54 elements, including the n = 44 ones
whose reference set is selected by the assignment, P(p <= alpha) <= alpha at
alpha = .05/.10/.20/.3333 under BOTH the default and the absolute-value conventions.
Nothing is anti-conservative.

## p11_joint.do -- the joint pre-trend kernel, independently recomputed
Recomputed the fixed-metric W_r series in numpy from the (already bitwise-verified) draws:
matches e(perm_jointW) to 1e-9, W_obs is exactly the 12th largest of 56, and
e(joint_horiz) = 12/56 = .2142857142857143 exactly. The identity draw carries W = W_obs
exactly, so the exact-share/membership logic is right here too.

---

## LATENT, NOT REACHABLE -- recorded for a successor
`wooldid_permtail` takes the plain share whenever `exact | memb` (ado:15992), where
`exact` is COMPLETENESS; but `wdp_addone` (ado:16147) is set on `!memb` alone, and its
caution (ado:8297) then says the element's "p-values use the add-one form (1 + count)/
(n + 1) rather than an exact share". If an element could ever have exk = 1 AND memb = 0,
the caution would be false AND the p would be a plain share over a reference set that
does not contain the observed statistic -- which can be 0 when the observed statistic is
extremal. I could not reach it: memb was TRUE for every one of 54 + 30 + 16 + 1 elements
across six designs, and the two estimation paths the brief names as non-bit-reproducible
(`unconditionalse`, `customvce(cluster ...)`) are REFUSED with permute() (E10-E12), which
removes the obvious trigger. Anyone who finds a permute()-compatible estimation path whose
rep-0 comparator and identity replication can differ by one ulp should re-open this.

## p12_permtail_unit.do -- Mata unit test, filed as FINDING PX-2
The "latent" item above is now demonstrated at the kernel level (exact=1, memb=0 gives
p = 0 where the caution promises 2/(n+1) = .0645). Reachability from the command line is
still unproven; see FINDINGS.md PX-2.

## UNTESTED (for a successor)
- permute(probabilities): its own support count (wooldid_permxp) was never oracle-checked;
  permexact is refused there by design, so the M-vs-support question is separate.
- permexact + contreat() / lattice() -- the continuous-treatment cfx cell kernels
  (wooldid_permcellwire) got ZERO coverage.
- permstepdown's per-subgroup matrices (_s#) and the att_it family (_i) were never formed
  in my designs.
- The adversarial p under permexact.
- Designs with 1000 < M <= 10000 were checked only at the enumerator level (p02 c13,
  M = 7560 reconciles), never run end to end, so the "this may be slow" caution above 1000
  and the behaviour of a 7560-replication enumeration are untested.
- permexact + poisson was run (rc 0, no add-one) but never reconciled against brute force.
