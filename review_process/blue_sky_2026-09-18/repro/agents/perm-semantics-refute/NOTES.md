# SKEPTIC notes -- refuting PERMSEM-01 / PERMSEM-02

## Code reading (pinned build, ado line numbers)
- `wooldid_bootdist` (ado:14644): `ok = rowsum(A :< .) :> 0` -> e(permdist_b)/e(permdist_se)
  are COMPACTED to rows with any non-missing coefficient. Called for permute at ado:8191/8195.
- `wooldid_permqf` (ado:16193) and `wooldid_permadv` (ado:16386): `Wr = J(rows(A),1,.)`
  -> e(perm_jointW)/e(perm_advW) are at FULL accumulator length with missing rows.
- `wooldid_bootsavemk` (ado:14871): `st_addobs(R+1)`, `_rep = 0::R`, then
  `st_store(., idx, B[1,k] \ A[.,k])` -> permsave writes ALL R accumulator rows,
  failed ones all-missing. permsave calls it with 5 args (ado:8859), so NO jointW/advW
  columns in the permutation file (the bootstrap call at ado:6546 does pass them).
  ==> permsave's _rep IS the accumulator row index, i.e. exactly the row basis of
      e(perm_jointW)/e(perm_advW).  The prober's "neither is indexed by permsave's _rep"
      is WRONG for the jointW/advW pair.
- `wooldid_permtail` (ado:15945): `out[6] = (exact ? sqrt(variance(v)*(n-1)/n) : sqrt(variance(v)))`
  where `exact` = exk = COMPLETENESS and `memb` = MEMBERSHIP governs the share (BS-23).
  out[6] -> e(perm_sd) (ado:8143) and e(perm_sd_`est') (ado:222), the latter PRINTED in the
  footer at ado:11734.
- `wooldid_permstep` (ado:16272) has the SAME completeness-tied population rescale
  (`if (exk) sd = sd :* sqrt((rows(CC)-1)/rows(CC))`) on the printed t_obs row of
  e(perm_stepdown) -- a common scalar per column so ranks/p are unaffected.

## Probe A -- pA_permsd.do / pA_permsd.log  (PERMSEM-02, own data)
9 counties, 2010-2019, 3 treated 2015, subgroups(big) on 3 units, permute(shuffle) permexact.
M=84, completed=84, failed=0.
RESULT: every COMPLETE column (n=84) posts the POPULATION SD (divisor n); every INCOMPLETE
column (n=83, the SG0 family) posts the SAMPLE SD (divisor n-1). ratio samp/pop = 1.006006
vs 1.006079 -- the incomplete rows sit 0.6% above the convention the complete rows use.
EVERY column's p is the EXACT SHARE (P_IS_EXACT_SHARE), complete and incomplete alike;
e(perm_addone) is empty; e(perm_incomplete) names the SG0 elements.
==> PERMSEM-02 CONFIRMED on independent data and independent geometry. Severity low.
Nuance: "no marker" is overstated -- e(perm_incomplete) + e(perm_n) let a user tell which
elements took which branch. What is missing is the DOCUMENTED RULE, not the marker.

## Probe B/B2/B3 -- pB_rowalign, pB2_rowalign, pB3_rowalign (PERMSEM-01, own data)
- pB  (16 counties, 4 late reporters, shuffle, 100 reps): failed=0. Everything aligned
  100/100/100, permsave 101 rows. So with no failures there is NO misalignment.
- pB2 (18 counties, four distinct onsets, late reporters AND early leavers): failed=0 too.
  wooldid tolerates a cohort with no pre-period / no post-period without failing the rep.
- pB3 (10 districts, 2 adopt 2017, 3 stop reporting after 2014; a shuffle that puts BOTH
  onsets on early leavers leaves no treated observation -> that rep fails; P = 3/45):
  requested=100 completed=93 FAILED=7.
    e(permdist_b) 93 x 6, e(permdist) 93 x 1        -> claim (i) CONFIRMED
    e(perm_jointW) 100 x 1 with exactly 7 missing   -> claim (ii) CONFIRMED
    PB_ALIGN agree=100 of 100: the missing rows of e(perm_jointW) are EXACTLY the
      _rep rows of the permsave file whose coefficient columns are all missing
      -> claim (iii) "neither is indexed by permsave's _rep" is REFUTED.
    e(permdist_b) == the compacted permsave coefficient block, maxabsdiff = 0
      -> permdist_b is recoverable from the file too; the file carries the whole key.
    Naive (1+#{W>=Wobs})/(N+1) on the full-length e(perm_jointW): p = 0.5148 vs the
      posted/correct 0.4787; overcount = 7 = e(perm_failed) EXACTLY -> claim (iv) CONFIRMED.
  TRAP for anyone repeating this: with permsave() the file ALSO carries onset_<id> columns,
  which are populated on FAILED replications too. "Row is all missing" must be evaluated on
  the coefficient columns only.

## Probe D -- pD_draworacle.do/.log  (NEW, negative result: the draws are faithful)
Took three successful replications out of the permsave file, relabelled the onset variable
by hand from the onset_<id> columns, refitted wooldid with no permute(), and compared:
  rep 1  reldiff 1.04e-15 ; rep 17  reldiff 0 ; rep 57  reldiff 0
The failed replication 24 (onsets on districts 6 and 10, both early leavers) errors r(459)
on a hand refit -- the failure is real, not an accounting artefact.
_rep = 0 row == e(perm_comparator_b) bit for bit.

## Probe C -- pC_conventions.do/.log  (NEW, four clean negatives)
 C1 multi-onset enumeration: onsets {2014,2016,2016} on 7 clusters (4 never treated)
    -> e(perm_npossible) = 105 = 7!/(1!2!4!), completed = 105. The log-factorial counter
    and the recursive enumerator (two independent algorithms, ado:13773 flags the hazard)
    agree on a multiset design.
 C2 e(perm_p) lies inside e(perm_p_ci) for every element on a 200-rep sampled run (0 of 4 out).
 C3 the permabsp / permrawp / permside() mutual-exclusion guard fires r(198) on all three pairs.
 C4 e(perm_p) == min(1, 2*min(p_lo, p_hi)) for every element (0 of 4 bad).
 C5 under permexact + membership, e(perm_p_se) and e(perm_p_ci) are withheld (0 of 4 posted).

## Probe E -- pE_onsetcohorts.do/.log (NEW ground, negative result)
permonsetcohorts (T1, 2026-09-16, assessment suite never run): 20 plants, 4 onset cohorts
of 5, i = cohort, cluster(plant) strictly finer. e(perm_irebuilt)=1; with cluster()==i it is
0 (documented no-op). Hand-rebuilt three drawn assignments by the documented map and
refitted: reldiff 0 / 4.9e-16 / 0. A refit of the same drawn onsets WITHOUT rebuilding i
errors r(459) each time -> the rebuild is load-bearing, not cosmetic.

## Probe F -- pF_permsave_failedrow.do/.log (NEW FINDING PERMREF-01)
sthlp:493 says a failed replication's permsave row is "all missing". It is not: the
onset_<id> columns carry the drawn assignment on failed rows too.
  PF_FAILED_BY_HELP_READING = 0 of 7 ; PF_FAILED_BY_COEF_COLUMNS = 7 of 7.
_rep 24: all six coefficient columns missing, onset_6 = onset_10 = 2017.

## FINAL STATE
FINDINGS.md written. Verdicts: PERMSEM-01 PARTIAL (low) -- the two row bases are real and
the over-count arithmetic is exact, but perm_jointW/perm_advW ARE indexed by permsave's
_rep and permdist_b is recoverable from the same file, so the correct framing is a
documentation gap about missing rows in e(perm_jointW). PERMSEM-02 CONFIRMED (low) with
three scoping corrections. One new finding (PERMREF-01) and six negative results.
