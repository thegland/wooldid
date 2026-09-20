# SKEPTIC / perm-exact-refute -- running notes

Task: refute or confirm PERMEX-PX-1 and PERMEX-PX-2; hunt new defects.
Build: /private/tmp/.../bluesky/build/wooldid.ado  md5 444892481be8d50ec962a102d6abef51

## Code facts established by reading (no run needed)
- ado:15992 `if (exact | memb)` selects the PLAIN share in wooldid_permtail.
  `exact` = COMPLETENESS (n == M); `memb` = membership of the observed row.
- ado:15974-15979 comment says: "`exact' is COMPLETENESS and still governs the
  population-SD convention below. `memb' is MEMBERSHIP ... which is what the SHARE
  convention turns on."  -> the code ORs in `exact`, contradicting its own comment.
- ado:16127 `exk = (M < . & sum(s :< .) == M)`; ado:16130 flags wdp_addone on
  `M < . & !memb & sobs < . & sum(s :< .) >= 2`  (i.e. on !memb ALONE).
  So exk=1 & memb=0 => caution printed AND plain share used. Contradiction is textual.
- ado:16085 comment: "an element gets plain shares only when its n equals M" -- that is
  the `exk` half; it omits the memb half.  Both comments are half-right.
- ado:15940 permmemb returns 0 when M >= . (sampled draws) -> under sampled draws memb
  is ALWAYS 0, so the add-one branch is taken.  PX-1's scope (permexact only) is right.
- ado:16290-16320 wooldid_permstep uses the same `exk | memb` rule for all three p rows.
- ado:200 `wdp_OR = cond(permexact != "" & wdp_r_obsrow != "", wdp_r_obsrow, "0")`
  and ado:13642 defaults wdp_obsrow to the STRING "0" -> non-empty, so orow=0 is passed
  when the identity column is not located.  Identity search at ado:13714-13720.
- ado:13692-13706 comment: the assignment COUNT (log-factorial) and the ENUMERATION
  (wooldid_permx_enum) are two independent algorithms and "nothing reconciles them".

## Probes
- r1_px1.do  : permexact vs permreps, step-down p rows; does the footer say "add-one"?

## Results so far
- r1b_px1.do / .log  PX-1 REPRODUCED on my own data: permexact -> step-down p rows are
  exact multiples of 1/30 (plain share); permreps(199) -> exact multiples of 1/200
  (add-one). Footer printed verbatim: "p_raw is the two-sided |t| share (add-one) in the
  step-down's own metric".  So the label is unqualified where the code is exactness-aware.
  BUT the number printed (1/30) is the statistically CORRECT exact p; 2/31 would
  double-count the observed assignment, which is already in the enumeration.
- r2_px2.do / .log  PX-2 reachability NOT reached on 3 further designs:
  T1 permonsetcohorts+permexact (i rebuilt per draw, e(perm_irebuilt)=1, M=60):
     bitwise hits of the comparator among the 60 draws = 1 for every element.
  T2 permtstat+permexact (separate b and se paths): hits = 1 for every element.
  T3 incomplete enumeration (unbalanced panel, M=20, n=16 on two horizons):
     exk=0 reached, memb=1, the "not estimable in every enumerated assignment" caution
     fired; the wdp_addone caution did NOT.
  The wdp_addone caution ("did not carry ... add-one form") printed 0 times in any run.
- STRUCTURAL ROUTE CLOSED: onset varying within a cluster (the obvious way to make the
  identity replication differ from the rep-0 comparator) is a HARD guard at ado:5978
  ("permute() requires clusters to be subsets of treatment-onset groups"), r(459).
- r3_enum.do + r3_audit.py  NEW, CLEAN: independent audit of the permexact enumeration
  on an INTERLEAVED-strata design (stratum = mod(unit,2), 8 clusters, M = 12*4 = 48).
  wooldid reps = 48 = independent log-factorial count; all 48 distinct; observed present
  exactly once; the enumerated SET equals a Python brute force exactly; no onset leaks
  across strata (NFPERM-02's row-keying holds end to end).
- r4_oracle.do + r4_oracle.py  NEW, CLEAN: recomputed e(perm_p), e(perm_p_lo/hi),
  e(perm_sd), e(perm_ci), e(perm_n) AND all four rows of e(perm_stepdown) from the saved
  draws in Python. Every one matches to <=1e-14 relative.
- NEW FINDING (mine): sthlp:607, the row of the p-value MAP for permute()+permexact,
  states "{bf:Floor 1/M}". The DEFAULT (equal-tailed 2 min) permexact floor is 2/M --
  measured 2/48 in r4 and 2/30 in r1b -- and the help's own floors paragraph (sthlp:512)
  says "2/M under the default and 1/M under the absolute-value conventions".
- r5_floor.do / .log  measured the permexact floor by convention at M=60:
  default 2/60, permabsp 1/60, permside(right) 1/60, permrawp 1/60; joint pre-trend
  p = 0.9833 = 59/60 (plain share, floor 1/M, where the map's joint row says 1/(n+1)).
- r6_guard.do / .log  MEASURED the mixed-onset-cluster guard: r(459) under default safety
  AND under safety(off); an onset-homogeneous control runs (M=6).
- Prior-art check for the floor finding: the project's own records (permute_coverage
  2026-09-02 l.830; engine_parity 2026-09-08 l.40,78) say 2/M. The 2026-09-10 help sweep
  logged the map's permexact row by line number without quoting its floor. Not ruled.
STATUS: complete. Verdicts: PX-1 CONFIRMED (low), PX-2 PARTIAL (low). New: SKEP-FLOOR-1.
