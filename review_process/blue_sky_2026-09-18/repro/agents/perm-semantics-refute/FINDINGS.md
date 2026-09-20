# SKEPTIC report -- refutation pass on PERMSEM-01 / PERMSEM-02

Build: scratchpad/bluesky/build/wooldid.ado (sha c1d8979), loaded with `adopath ++`.
All data generated here; the project test suite was not touched.
Every probe do-file and its .log are in this directory.

---

## VERDICT PERMSEM-01 -- PARTIAL (severity stays low)

My repro: pB3_rowalign.do / .log (controls: pB_rowalign.do, pB2_rowalign.do -- zero
failures, zero misalignment).

Failure mechanism, mine and realistic: 10 districts, 2 adopt in 2017, 3 stop reporting
after 2014 (attrition). A shuffled assignment landing both onsets on early leavers leaves
the replication with no treated observation, so it fails. P = C(3,2)/C(10,2) = 3/45.
No permprobs. Measured: requested 100, completed 93, FAILED 7.

CONFIRMED
  e(permdist_b)   93 x 6    compacted to completed replications
  e(permdist)     93 x 1    same compacted basis
  e(perm_jointW) 100 x 1    full accumulator length, 7 missing rows
  e(perm_advW)   100 x 1    full accumulator length
  permsave file  101 rows   _rep 0..100, full accumulator
The "sharper case" arithmetic is exactly right: (1+#{W_r>=W_obs})/(n+1) applied to the
posted full-length e(perm_jointW) WITHOUT dropping missing rows gives p = 0.51485 against
the posted/correct 0.47872 -- an over-count of 7 = e(perm_failed) EXACTLY (Mata's :>=
ranks missing above every finite value).  See the PB_RECON lines.

REFUTED -- "neither is indexed by permsave's _rep".
False for the jointW/advW pair, which is the half that matters because that pair is the
help's advertised reconstruction route. wooldid_bootsavemk (ado:14871) does st_addobs(R+1),
_rep = 0::R, st_store(., idx, B[1,k] \ A[.,k]) -- the file is the whole accumulator -- and
wooldid_permqf/permadv write Wr = J(rows(A),1,.) on that same accumulator. Measured row by
row, not by counts:
    PB_ALIGN agree=100 of 100 ; JWmiss=7 FILEmiss=7
    PB_VERDICT3 perm_jointW IS row-indexed by permsave _rep
The missing rows of e(perm_jointW) are EXACTLY the _rep rows of the permsave file whose
coefficient columns are all missing. And e(permdist_b) is recoverable from the same file:
PB_EQ permdist_b == compacted permsave coefficient block, maxabsdiff = 0. So the run does
supply one row basis for everything -- the permsave file -- and it reconciles the two e()
bases.

CORRECTION. This is not a row-alignment defect against permsave; it is a DOCUMENTATION GAP
ABOUT MISSING ROWS. sthlp:1348-1349 advertises e(perm_jointW)/e(perm_jointWobs) with "the p
is reconstructable from the pair" and never says the first is at full REQUESTED length with
missing rows while e(permdist_b) beside it is compacted. A user doing the reconstruction the
sentence invites gets a silently different p whenever any replication failed. select() on
the mismatched pair errors loudly, so the only quiet shapes are a positional loop and a
by-row merge -- the prober had that right.

ENGINE-AGNOSTIC (source read only): the same gap sits on the bootstrap side. sthlp:1340
describes e(boot_jointW) as "the replication joint statistics, one named column per test"
with no missing-rows statement either, and wooldid_bootqf (ado:16345) writes
Wr = J(rows(A),1,.) exactly as wooldid_permqf does. The prober's contrast -- "the help
states this for the bootstrap side" -- is about bootsave()/e(bootdist_se), not about
e(boot_jointW). Not measured here.

REALISM: true but narrow -- needs failed replications AND a hand reconstruction. The
prober's own route (permprobs with a 14% chance of zero treated clusters) is contrived;
mine (administrative attrition) is not. Severity low stands.

---

## VERDICT PERMSEM-02 -- CONFIRMED (severity stays low)

My repro: pA_permsd.do / .log. Own data, own geometry: 9 counties 2010-2019, 3 treated in
2015, subgroups(big) on 3 counties, permute(shuffle) permexact. M = C(9,3) = 84,
completed 84, failed 0.

  main_att              n=84=M   posted +1.7784ba9c41b5cX-001   POPULATION 1/n
  pretreatment_att      n=84     posted +1.5d2808019274dX-001   POPULATION
  ES_att_post0/post1    n=84                                    POPULATION
  main_SG0_att          n=83     posted +1.a5d4f9434c5dfX-001   SAMPLE 1/(n-1)
  pretreatment_SG0_att  n=83                                    SAMPLE
  ES_att_SG0_pre2/pre1/post0/post1  n=83                        SAMPLE
(the four columns tagged NEITHER in the log sit 1 ulp off my own reconstruction, which
recomputes from the compacted e(permdist_b) WITHOUT sorting, while wooldid_permtail does
sort(v,1) before variance() -- summation order, nothing more. They are unambiguously on the
same side: the pop/sample gap is 0.6%, five orders of magnitude above a ulp.)

On the same columns every p is the EXACT SHARE (PA_P ... P_IS_EXACT_SHARE for all twelve),
with e(perm_addone) empty and e(perm_incomplete) naming exactly the SG0 family. So the same
reference set is a population for its p and a sample for its SD, inside one table.
Source: wooldid_permtail ado:15968,
  out[6] = (exact ? sqrt(variance(v)*(n-1)/n) : sqrt(variance(v)))
where `exact' is exk (COMPLETENESS) while `memb' (MEMBERSHIP) governs the share.

THREE CORRECTIONS
1. "no marker" is overstated. e(perm_incomplete) names the elements that lost assignments
   and e(perm_n) gives each element's n, so a user CAN tell which branch an element took.
   What is missing is the DOCUMENTED RULE, not the marker (sthlp:1345 says only "SD of the
   permutation distribution").
2. Scope it right. "sampled draws -> sample SD, complete enumeration -> population SD" is a
   defensible rule and should not be attacked. The defect is precisely that an INCOMPLETE
   ENUMERATION falls to the SAMPLED branch for the SD while BS-23 puts it in the EXACT
   branch for the p.
3. "feeds no p-value or interval, so nothing downstream moves" is right about inference but
   wrong about visibility. The same kernel row 6 is also posted as the legacy scalar
   e(perm_sd_<est>) (ado:222) and that scalar is PRINTED in the run footer (ado:11734,
   "att: estimate ... SD ..."). In my run main_att is complete so the printed number is the
   population SD. Whether a POOLED estimand can be incomplete while other elements are not
   was NOT established here: under permexact a whole-row failure costs every element the
   same rows, so exk = 0 everywhere and the table is uniformly sample-SD with no switch.
   Per-element loss on main_att would need it unestimable on a replication that otherwise
   estimated; I did not reach that case, and e(perm_incomplete) would name it if it arose.

SECOND SITE (source read only, not measured here): wooldid_permstep ado:16294 carries the
same completeness-tied rescale -- if (exk) sd = sd :* sqrt((rows(CC)-1)/rows(CC)) -- on the
printed t_obs row of e(perm_stepdown), while its p_raw/p_singlestep/p_stepdown rows use
(exk | memb). exk is one flag per stepdown table there, so the switch is across tables in a
run rather than within one, and ranks are unaffected (a common scalar per column).

---

## NEW FINDING (found while refuting, not seeded)

### PERMREF-01 [low/certain] permsave()'s failed rows are not "all missing", and the help says they are

Repro: pF_permsave_failedrow.do / .log

sthlp:493, describing permsave(): "(_rep = 1 to permreps(), or 1 to M under permexact; a
failed replication keeps its row, ALL MISSING) plus a _rep = 0 row carrying the observed
estimates". The same paragraph adds "one onset_<clustervalue> column per cluster holding the
onset that cluster carried in each replication", and those columns carry the drawn
assignment on a FAILED replication's row too. Measured on a file with e(perm_failed) = 7
of 100:
    PF_ALLCOLS = 16   PF_COEFCOLS = 6
    PF_FAILED_BY_HELP_READING  = 0    (every column missing)      <- finds none of the 7
    PF_FAILED_BY_COEF_COLUMNS  = 7    (coefficient columns only)  <- finds all 7
    PF_FAILED_REPS = 24 53 71 80 81 82 90
and _rep 24 lists as main_att . pretreatment_att . ES_att_* . beside onset_6 = 2017,
onset_10 = 2017 -- the assignment that failed, intact. Hand-refitting that assignment errors
r(459) (pD_draworacle.do, PD_FAILREFIT_RC=459).

Consequence: a user identifying failed replications the way the sentence invites finds none
of them on a file that has several -- and it is exactly the workflow PERMSEM-01 is about. It
cost me one wrong probe (pB2_rowalign.do reports PB_SAVE_FAILEDROWS = 0 because it swept the
onset_* columns into the test). Keeping the assignment on a failed row is arguably the
better behaviour; the sentence is what needs fixing, plus a note that the failed-row test is
on the coefficient columns. bootsave() at sthlp:557 says the same thing about a file that
also carries w_<cluster> resample columns; not tested here.

---

## NEGATIVE RESULTS (probed hard, clean)

N1  The permutation draws ARE the relabelled refits. pD_draworacle.do. Took three completed
    replications out of the permsave file, rebuilt the onset variable by hand from the
    onset_<id> columns, refitted wooldid with no permute(), compared main_att:
    reldiff 1.0e-15, 0, 0. The _rep = 0 row equals e(perm_comparator_b) bit for bit. The
    failed replication's assignment genuinely errors r(459) on a hand refit, so the failure
    accounting is real and not a bookkeeping artefact.

N2  permonsetcohorts (T1, shipped 2026-09-16, assessment suite deferred and never run)
    rebuilds i faithfully. pE_onsetcohorts.do. 20 plants in 4 onset cohorts of 5, i =
    cohort, cluster(plant) strictly finer. Rebuilt each drawn assignment's i by the
    documented map (drawn onset -> observed i code of the cohort with that observed onset;
    never-treated -> the never-treated cohort's code, ado:7578-7590) and refitted by hand:
    reldiff 0, 4.9e-16, 0 on three replications. e(perm_irebuilt) = 1 with cluster() finer,
    0 with cluster() == i (the documented no-op). A refit of the same drawn onsets WITHOUT
    rebuilding i errors r(459) every time, so the rebuild is doing real work.

N3  Multiset enumeration count. pC_conventions.do C1. The possible-assignment count and the
    enumerator are two independent algorithms and ado:13773 flags that nothing reconciles
    them. On onsets {2014, 2016, 2016} over 7 clusters (4 never treated),
    e(perm_npossible) = 105 = 7!/(1!2!4!) and the run completed 105 replications.

N4  p-value bookkeeping. pC_conventions.do C2/C4/C5.
      e(perm_p) lies inside e(perm_p_ci) for every element on a 200-rep sampled run (0/4 out)
      e(perm_p) == min(1, 2*min(p_lo, p_hi)) for every element (0/4 bad)
      under permexact + membership, e(perm_p_se) and e(perm_p_ci) are withheld (0/4 posted)

N5  The permabsp / permrawp / permside() mutual-exclusion guard fires. pC_conventions.do C3:
    r(198) on all three pairs, matching sthlp:1355.

N6  Replication failure is hard to provoke on realistic unbalanced panels. pB_rowalign.do
    (16 counties, 4 reporting only from 2016) and pB2_rowalign.do (18 counties, four
    distinct onsets, late reporters AND early leavers, singleton cohorts) both completed
    EVERY replication. wooldid tolerates a cohort with no pre-period and a cohort with no
    post-period without failing the replication. It took removing the treated observations
    altogether to make one fail. That bounds how often PERMSEM-01 can bite in practice.

N7  e(permdist) and e(permdist_b) keep rows on DIFFERENT tests and happened to agree.
    The legacy matrix (ado:7882) keeps a row on !missing(res[r,1]) -- the FIRST requested
    estimand -- while e(permdist_b) keeps on "any non-missing coefficient"
    (wooldid_bootdist, ado:14644). Both came out 93 x . in pB3, so they agree here. They
    would diverge only on a replication that estimated but lost its first pooled estimand;
    not reached. Flagged for the record, not filed.
