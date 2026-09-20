# FUZZ agent — NOTES (keep current)

WORKDIR: .../probes/fuzz
Build: adopath ++ .../bluesky/build  (verified `which wooldid` = pinned build)

## Harness
- `data/` 19 realistic datasets (ds_*.dta) + 6 zero-noise exactness datasets (ds_ex_*.dta).
  Built by build_data.do / build_exact.do from the recon generators. manifest.csv has counts.
- `assert.do` — wdfz_snap + wdfz_assert: automated invariants on every rc-0 run
  (SE>=0, p in [0,1], CI order/containment, t=b/se, e(b)/e(V) conformable+symmetric+PSD,
   e(N)>0, level name on CI column, data signature/sort/varcount unchanged).
- `runner.do` (+ mkbatch.sh, runbatch.sh) — generic spec-file driver, pipe-delimited spec:
  id|ds|y|i|ifclause|wclause|levtxt|eng|opts. Writes res_<tag>.csv and a full .log with
  @@BEGIN/@@END sentinels.
- `exact_runner.do` (+ truth.do, mkexact.sh, runexact.sh) — ORACLE arm. On zero-noise data
  the design truth of all four estimands is closed-form (truth.do = program form of the
  recon gen_truth.do). Spec: id|ds|y|i|if|wc|wvar|opts. Asserts |est-truth|/max(1,|truth|).
- `gen_fuzz.py` — random compat-aware option sampler (phase A/B specs).
- `gen_exact.py` — identification-preserving sampler (core) + identification-changing (observe).
- `classify.py` — parses a batch log, dedupes failures by (rc, first Error line),
  splits ATTRIB (wooldid "Error:" message) vs RAW (bare Stata error).

## Status
- A1 (50 runs): rc 0=33, 198=13, 459=4. ALL failures ATTRIB with informative wooldid errors.
- B1 (50 runs): rc 0=40, 198=5, 322=1, 459=4. ALL failures ATTRIB.
- X0 smoke: oracle reproduces design truth to 1e-15 on 10 baselines. ORACLE IS SOUND.
- B2..B5 running; X1..X3 (core exactness, 150 runs) + O1..O2 (observe, 60) running.

## Known false positives in my own assertions (do not report)
- `BV_DIM(bN/V0x0)` whenever e(properties)=="b": fastnose / permute-without-permanalyticse
  legitimately post no covariance; Stata hands back a 0x0 for e(V). Confirmed by
  probe_fastnose.do: test/lincom refuse with r(321), which is correct.
- `P_ZERO` on bootstrap PERCENTILE coefficient rows and under permrawp: ruled design (§2).

## Findings so far
- (see FINDINGS.md)

## Update (mid-run)
- 200 phase-A/B fuzz runs done (A1,B1-B3): rc 0=147, 198=29, 322=7, 459=17.
  EVERY non-zero rc carried an informative wooldid "Error:" line EXCEPT two:
  * B0062 r(459) `inconsistent estimation sample levels 0 and 1 of factor __000004`
    -- a raw Stata `margins' error carrying a wooldid TEMPVAR name.  Minimising:
    dropping `semielasticity' from the option list makes it go away.
  * B0061 r(198) `pretreatment_SG480000000_att_itime invalid name` -- name overflow past
    Stata's 32 chars from a 9-digit subgroup value.  NOT a finding: wooldid prints an
    explicit Caution predicting exactly this crash and naming the remedy.
- ORACLE arm X1+X2 (71 runs so far): every deviation from the design truth is an i=cohort
  run.  All 39+ i=unit runs match to <=2e-15.  i=cohort deviations are the documented
  unit-fixed-effect aggregation problem (gen_core header: "not a bug"); ex_noai datasets
  (gaisd=0) are queued as the control.
- probe_vneg: sthlp:1616-1617's claim about the -1e100 placeholder HOLDS exactly --
  test drops the constraint, lincom/nlcom propagate missing.  CLEAN.
- probe_fastnose: fastnose DOES post e(b_es_att) and e(numpre_att) and no e(V);
  sthlp:1325 listing fastnose under "covariance not returned, coefficients still are"
  is ACCURATE.  (§7 item 11 resolves CLEAN on this build.)
- probe_replay: `estimates replay <name>` fails r(198) with "command line input should
  feature all of y i t ttre" -- FINDING (repro/REPRO_fuzz_estreplay.do).

## Final tally (see FINDINGS.md for the write-ups)
Arms and totals:
 * A/B random option fuzz: 300 runs over 19 datasets. rc 0 = 217, 198 = 35, 322 = 11, 459 = 24 (approx).
   Exactly TWO failures out of ~83 were RAW (no wooldid "Error:"): the semielasticity family
   (FZ-2) and the 9-digit subgroup name overflow (documented by an explicit Caution, not filed).
 * X oracle arm: 160 zero-noise runs. EVERY i=unit run (113) reproduced all four design
   estimand truths to <= 2e-15. All 17 deviations were i=cohort, i.e. the unit-fixed-effect
   aggregation problem the generator's own header calls "not a bug".
 * O observe arm (30): timetrends, interactivecontrols(), coarsecohortcontrols(),
   esrelativeto(-2) and esfixedbaseperiod are all EXACT at i=unit -- they do not move the
   estimand on a TWFE DGP. impvalues moves it 5-15% (the documented cluster-equal-weight
   substitution). anticipation() moves it because the post set changes.
 * P metamorphic arm (100+ pairs so far): point estimates invariant to 1e-16 under row
   shuffle, reverse sort, id relabel, time shift, y shift, y negate, y rescale, weight
   rescale, if-vs-keep, a collinear control and row duplication. Three pairs moved DERIVED
   columns (se/t/p/CI) by 1e-9..2.6e-8 under sortrev/idperm/yscale; regtol(12)/regtol(13)
   does not shrink it and the point estimate stays exact, so it is the HDFE projection path,
   not an estimand bug. Reported as an observation, not a finding.
 * R refusal arm (80): every deliberately invalid combination was refused with an informative
   wooldid "Error:" line, except dots(-5) under permute() (FZ-5), cfxplottypes(zzz) /
   cfxlattice()/latticeignoreweights/makecfxplotsbysg/cfxplotnose/ccnoabsorb/contreatpoly()
   silently accepted when inert, and contreatelasticitytype() (FZ-1).
 * probe_gtwt (84 configs): every reported estimand equals sum(w_e*cell)/sum(w_e) over the
   e(gteffects) cells to 1e-15, all aggregation weights nonnegative, every estimand inside the
   convex hull of its cells. CLEAN.
 * probe_civ (94 configs, 1504 interval rows): (ub-lb)/(2*se) == invttail(e(df_r), alpha/2)
   exactly, at level 80/90/95/99, under iid/robust/cluster(region|unit|cohort). CLEAN.
 * probe_hull (69 configs): no estimand outside its cell hull. CLEAN.
 * probe_ident: adding estimands, subgroups(), diagnostic flags and level() never move the
   other rows; subgroups(sgall) reproduces the main table exactly. The ONE mover is the
   requested event-study window (see the note in FINDINGS).
 * probe_vneg: sthlp:1616-7's -1e100 placeholder claim holds exactly (test drops, lincom
   propagates missing). CLEAN.
 * probe_fastnose: fastnose posts e(b_es_att)/e(numpre_att) and no e(V); sthlp:1325 accurate.
