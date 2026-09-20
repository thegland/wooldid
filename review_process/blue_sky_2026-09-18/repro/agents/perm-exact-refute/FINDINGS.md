# SKEPTIC report -- perm-exact-refute
Build: .../bluesky/build/wooldid.ado  md5 444892481be8d50ec962a102d6abef51
Every claim below is reproduced from its own terms on my own data; I did not run the
original prober's do-files as evidence.

## VERDICTS ON THE SEEDED FINDINGS

### PERMEX-PX-1 -- CONFIRMED (low, documentation-only). Correction to the framing.
My repro: r1b_px1.do (self-contained; 5 clusters, onsets 4,4,6,.,. -> M = 30).
 - permexact leg: e(perm_stepdown) p_raw / p_singlestep / p_stepdown are
   .06666666666666669 and .03333333333333330 -- EXACT multiples of 1/30 (p*30 = 2.000 and
   1.000 to 10 decimals; p*31 = 2.0667 and 1.0333). Plain share.
 - permreps(199) leg: the same three rows are .045 and .015 -- EXACT multiples of 1/200
   (p*200 = 9.000 and 3.000; p*199 = 8.955 and 2.985). Add-one.
 - The footer printed verbatim, in the same run as the permexact leg:
     "... where p_raw is the two-sided |t| share (add-one) in the step-down's own metric,
      not the table's convention ..."
   and sthlp:514 says "the two-sided |t| share, add-one".
So the convention really does switch under permexact and the sentence really is
unqualified. The claim is true.
CORRECTION to the prober's framing (why I do not raise the severity):
 - The number wooldid prints is the statistically CORRECT one. Under a complete
   enumeration the observed assignment IS one of the M draws, so nex already counts it;
   (1+nex)/(n+1) = 2/31 would count it twice. The code comment at ado:16290 says exactly
   this. The prober's EXPECTED says so too, but its OBSERVED leads with 0.0645161 as "the
   value the footer and help promise", which reads as though the printed number were wrong.
 - The practical floor gap is 1/30 = .0333 (actual) versus the generic add-one floor
   1/31 = .0323 a reader would infer, i.e. 3%, not "a factor of two". The 2/(M+1) = .0645
   figure is what the code would print if the branch flipped, not what a reader infers
   from the word "add-one".
 - The footer states the governing rule FOUR lines above the offending sentence: "The
   assignments are enumerated rather than sampled, so the p-values are exact shares with
   no Monte-Carlo error" (verified in r1b_px1.log).
Severity stays LOW; documentation-only; the class is right.

### PERMEX-PX-2 -- PARTIAL. Code contradiction CONFIRMED; runtime reachability REFUTED
###                 on every route I could construct, including the one the prober missed.
My repro: r2_px2.do (self-contained, three designs).
The textual half stands and needs no run:
 - ado:15992 `if (exact | memb)` takes the plain share; ado:16127 `exk = (M < . & sum(s :< .) == M)`
   is COMPLETENESS; ado:16130 flags wdp_addone on `!memb` ALONE; ado:8297 then asserts the
   add-one form is in use. ado:15974-15979's own comment says "`memb' is MEMBERSHIP --
   which is what the SHARE convention turns on", and the code ORs `exact` in anyway.
The reachability half is where I differ. Three further designs, none reached it:
 - T1 permonsetcohorts + permexact with cluster strictly finer than i (i rebuilt in every
   replication, e(perm_irebuilt)=1, M=60): the rep-0 comparator statistic occurs bitwise
   exactly ONCE among the 60 draws for every element -> memb = 1.
 - T2 permtstat + permexact (b and se travel separate paths to s = b/se): same, hits = 1.
 - T3 unbalanced panel, M=20, two horizons estimable in only 16 assignments: exk = 0 with
   memb = 1 is REACHED (the "not estimable in every enumerated assignment" caution fires),
   which is the reachable half of the OR and shows membership detection works end to end.
 - The wdp_addone caution printed ZERO times across all runs in this directory.
AND I closed, BY MEASUREMENT (r6_guard.do), the one structural route that would have
reached it: an onset varying WITHIN a cluster would make the identity replication (which
writes one cluster-level onset to every row of the cluster) fit different data from the
rep-0 comparator (which fits ttrevar as given), so the observed statistic would not be a
member of its own reference set while the enumeration stayed complete. Measured on
cluster(state) with two units of one state at onsets 4 and 6:
   G1 default safety  -> r(459) "permute() requires clusters to be subsets of
                                 treatment-onset groups"
   G2 safety(off)     -> r(459), the SAME error; the guard is not liftable
   G3 control (same state made onset-homogeneous) -> runs, e(perm_exact)=1, M=6
The guard is ado:5975-5982, on the design sample, and covers the treated/never-treated-mixed
case as well.
Correction: PX-2's OBSERVED is a unit test of kernels called with inputs wooldid does not
produce. That is legitimate evidence about the kernel and no evidence about the command.
Severity LOW (latent inconsistency between three code sites; no user-visible consequence
demonstrated by the prober or by me on 3 + 6 = 9 designs).

## NEW FINDINGS (mine)

### SKEP-FLOOR-1 -- the p-value MAP states the wrong floor for permexact's DEFAULT
###                 convention: it says 1/M, the measured floor is 2/M.
Repro: r5_floor.do (self-contained; 6 clusters, onsets 5,5,7,.,.,. -> M = 60).
sthlp:607, in the table the help introduces with "This table is the map" under the column
header "convention, and the smallest p it can print":
   "permute() + permexact:  Plain shares over the M enumerated assignments ... Floor 1/M."
Measured on M = 60, a design whose observed statistic is the extreme of its own enumeration:
   default (no convention option): e(perm_p) main_att = .03333333333333330 = 2/60
                                   e(perm_p_hi)      = .01666666666666670 = 1/60
   permabsp:                       e(perm_p)         = .01666666666666670 = 1/60
   permside(right):                e(perm_p)         = .01666666666666670 = 1/60
   permrawp:                       e(perm_p)         = .01666666666666670 = 1/60
The default is 2 min(P_lo, P_hi) on plain shares, so 2/M is a hard floor there; 1/M is not
attainable. The help's own prose has it right at sthlp:512 -- "under permexact with M
assignments it is 2/M under the default and 1/M under the absolute-value conventions" --
so the MAP contradicts the prose, and the map is the surface a user consults. Corroborated
independently at M = 30 (r1b_px1.log, att p = .0667 = 2/30) and M = 48 (r4_oracle.py,
p = .0416667 = 2/48).
Why it matters: permexact is what wooldid's OWN runtime caution recommends when fewer than
six clusters are ever treated (ado:6070, "permexact enumerates the assignments and reports
the exact floor"), i.e. exactly where M is small and the floor binds. At M = 30 the map
promises a reachable .033 where the truth is .067 -- the two straddle .05, so a user
sizing a randomization-inference design is told it can reach significance when it cannot.
Decisive corroboration from wooldid's own runtime text: the incomplete-enumeration
caution at ado:7876 says "the smallest attainable p is 2/n under the default convention
(1/n under the absolute-value and one-sided conventions) rather than 2/`wdp_R'" -- i.e.
the CODE states the default complete-enumeration floor as 2/M, convention-aware, in the
message it prints to the user. Only the help's map row says 1/M.
Secondary, same row-family, weaker: the map's "permute(), joint and equality tests" row
says add-one with floor 1/(n+1) unconditionally; under permexact ado:16244 uses nex/n
(measured: the joint pre-trend p in r5_floor.log run A is 0.9833 = 59/60, a multiple of 1/M, floor 1/M not 1/61).
Realistic: yes. Severity: MEDIUM, and the reason it outranks PX-1 although both are
documentation-class: PX-1's unqualified word "add-one" changes no number a user computes
(the printed p is correct and the inferred floor is 1/31 against an actual 1/30, a 3%
gap), whereas this row IS a number, is wrong by a factor of two in the OPTIMISTIC
direction, sits in the table the help itself designates as the map under a column header
reading "the smallest p it can print", and at realistic M straddles .05. Confidence:
certain (measured at M = 30, 48 and 60, on three separate designs).
PRIOR-ART CHECK: not a rediscovery. review_process/permute_coverage_2026-09-02.md:830 and
research_2026-09-08/engine_parity.md:40,78 all state the exact floor as 2/M, i.e. the
project's own record agrees with the measurement and against the map row. The 2026-09-10
help sweep (review_process/research_2026-09-10/sweep_audit_raw.json:148) enumerates the
map's thirteen rows and quotes the floors of the others, but records the permexact row by
line number only -- it was audited for existence, not for correctness.

## NEGATIVE RESULTS (attacked hard, clean)

### N1 -- the permexact enumeration is exactly right, including under crossing strata.
r3_enum.do + r3_audit.py. wooldid's own comment (ado:13692-13706) says the possible-
assignment COUNT (log-factorials) and the ENUMERATION (a separate recursion) are
"two INDEPENDENT algorithms" and "nothing reconciles them", and names the two failure modes
(all-control replications, or a non-exhaustive enumeration still flagged exact). I audited
them against each other and against a Python brute force on an INTERLEAVED-strata design
(stratum = mod(unit,2) over 8 clusters, so ascending cluster order crosses stratum order --
the NFPERM-02 hazard): M from the formula = 48 = wooldid's reps = brute-force count; all 48
distinct; the observed assignment present exactly once; the enumerated SET equals the brute
force exactly; and no replication moves an onset across strata. Nothing wrong.

### N2 -- every reported permexact quantity reproduces from the draws to 1e-14.
r4_oracle.do + r4_oracle.py. Recomputed, in Python, from the permsave() draws alone:
e(perm_n), e(perm_p), e(perm_p_lo), e(perm_p_hi), e(perm_sd) (population SD under the
exact branch), e(perm_ci) (wooldid_permq's quantile rule, placed about the displayed
estimate) and ALL FOUR rows of e(perm_stepdown) (fixed-metric studentization, single-step
max and the monotone step-down). Every value matches. If the observed statistic were
ranked on a different path from the draws, or the step-down's ordering or monotonicity
were wrong, this would have caught it.
