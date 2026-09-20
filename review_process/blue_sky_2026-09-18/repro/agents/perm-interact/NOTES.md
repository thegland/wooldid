# perm-interact — permute() x everything else (interaction matrix, hunting SILENT NO-OPS)

Agent label for finding ids: `PX-`.

## State (live)
- Build asserted path: pinned build at .../bluesky/build (p01 log line `which wooldid`).
- Data: `base.dta` built by `mkdata.do` from recon gen_core (40 units, 8 periods,
  onsets 2003/2005, dose, unit weights, count outcome, zbig strong covariate).
- Harness: `pchk.do` defines `pchk <tag> [, opts() wt() ifq() yv() reps() iv() cl()]`.
  Every call is `wooldid <y> unit t ttre [wt] [if] , att cluster(unit) permute(shuffle)
  permreps(30) rseed(777) <opts>` and prints, per tag:
    PCHK|tag|RC=..|SEC=..
    PCHK|tag|NCMP=..|MAXRELDIFF=<hex>|AT=<colname>     <- comparator check (see below)
    PCHK|tag|PDIST=RxC|NMISS=..|S1=<hex>|S2=<hex>      <- null-distribution fingerprint
    PCHK|tag|P_ATT=<hex>|B_ATT=<hex>|SD=<hex>

## THE PRIMARY DETECTOR (novel, not seeded)
`wooldid_permrun` runs a mandatory **rep 0 on the TRUE assignment** under the rep option
string `wdp_ro`, and posts that vector as `e(perm_comparator_b)` (ado:8169). The code
comment at ado:7660ff states that under the fast-mode default "the comparator then equals
it [the observed e(b)] bitwise". **Nothing asserts it.**
So: `max |e(perm_comparator_b) - e(b)| / |e(b)|` over name-matched columns is a direct,
per-option test of `repopts == the user's options`. Any nonzero value under the
fast-mode default = the replications are fitting a DIFFERENT MODEL than the table shows.
This catches dropped options AND mutated ones, with no RNG dependence.
Secondary detector: toggle option X on/off at identical rseed and compare the
`PDIST` fingerprint; identical fingerprint + different `e(b)` = X never reached the reps.

## Code-reading results so far (to be confirmed by running)
- ado:125-161 (permute rep string) vs ado:321-357 (bootstrap rep string). Diff:
  * bootstrap carries `contreatwithin`; permute does not — but **permctl REFUSES
    contreatwithin outright** (ado:5832), so this is a correct refusal, not a gap.
  * bootstrap carries `gteffects`; permute does not — deliberate: cells are keyed by
    cohort onset which the reps reassign. Display line ado:11726 + sthlp:1150 document it.
  * permute carries `customvce()` (ado:128) while bootstrap REFUSES it (ado:6414).
    BUT: validate refuses `customvce()` together with `cluster()` (ado:4667) and permctl
    REQUIRES a single `cluster()` (ado:5966). => permute+customvce looks UNREACHABLE and
    ado:128 looks like DEAD CODE. Must confirm empirically (p03: cvce_nocl / cvce_cl).
- permctl refusals found by reading: impvalues, unconditionalse, contreatwithin,
  contreatcontrols(untreatedZero), oldsyntax, pretrendpower(), suppresssgprimaryeffects,
  multi-way cluster, permtstat/permtci+fastnose, permreps<20, degenerate design.
  simultaneousbands is refused in the worker at ado:4906.
- e(V) under permute: permwire posts NO covariance of its own. Fast mode => no e(V) at all
  (`wdid_vcvstate` = absent, e(vce)/e(vcetype) cleared). Under permtstat/permanalyticse
  e(V) and e(V_es_*) stay ANALYTIC, and the part(perm) footer says so verbatim
  (ado:11676-11678: "e(V) and e(V_es_*) hold the analytic covariance, so test or lincom
  run by hand on them are analytic, not permutation, quantities"). Seeded lead (2) looks
  CLEAN on documentation; still to verify the line prints and that honestdid is not
  silently fed an analytic VCV without that line being visible.
- ES-window forwarding: `esprelength` is emitted RAW (A8) because ado:5153 bumps it under
  esfixedbaseperiod; `espostlength`/`esrelativeto` are emitted post-validate but the two
  mutations that exist (ado:4851 and ado:4856) are idempotent under the re-parse. Expect
  clean; check with the comparator under esf + subgroups + esprelength.

## Harness note
Project Stata cap of 8 was saturated for a long stretch; p01 queued behind it.

## RESULT p01 (2026-09-17) — 14 configs, ALL CLEAN
`grep "^PCHK|" p01_base.log`. Determinism verified (det1 == det2 bitwise on every
surface). Comparator `max|e(perm_comparator_b) - e(b)|/|e(b)|` = **exactly 0** on all 14
configs (`controls()`, `fe()`, `timetrends`, `subgroups()`, ES window, `esfixedbaseperiod`,
`esrelativeto()`, `anticipation()`, `oostreatedcontrols()`, `robust`, `nofastcollapse`,
`customweightmultiplier()`, `gteffects`). No option reaches the replications in mutated form.
Null-distribution fingerprints MOVED (option reached the reps) for: controls, timetrends,
subgroups, esprelength/espostlength, esfixedbaseperiod, customweightmultiplier, anticipation.
Fingerprints did NOT move for fe(region), oostreatedcontrols(2), robust, nofastcollapse,
gteffects, esrelativeto-without-esf — and in every one of those the OBSERVED e(b) did not
move either, so they are genuine model no-ops on this data, not forwarding gaps
(esrelativeto prints "Caution: esrelativeto() has no effect without esfixedbaseperiod").

## RESULT p02/p04/p12 — comparator clean everywhere (65+ configs)
p02 (weights aw/pw, `if`, level(90), each estimand and all four together,
interactivecontrols, safety(off), regtol(6), emptycellsoverride, jointtests,
adversarial, permstepdown): MAXRELDIFF exactly 0 on every config that ran.
p04 (21 configs with an ES window, incl. subgroups x esfixedbaseperiod, anticipation,
timetrends, oostreatedcontrols, controls, coarsecohortcontrols, fe on a TIME-VARYING
group, interactivecontrols, aw, `if`, customweightmultiplier, nofastcollapse, all four
estimands): MAXRELDIFF exactly 0 on all, EXCEPT the deliberate positive control
`es_ase` (permanalyticse) at 7.3e-11 on a near-zero pre-period coefficient — which is the
documented fastnose-vs-full F9 contract and proves the detector has power.
p12 (cluster() STRICTLY COARSER than i: 20 paired clusters, unit-level i — the geometry
permute() exists for): same sweep, MAXRELDIFF exactly 0.
Refusals verified informative: fw (users may not pass fweights), multi-way cluster,
impvalues, unconditionalse, contreatwithin, contreatcontrols(untreatedZero),
suppresssgprimaryeffects, pretrendpower(), oldsyntax, simultaneousbands,
permonsetcohorts+subgroups, cluster-not-nesting-i without permonsetcohorts,
permabsp+permrawp, fastnose+permtstat, fastnose+permtci, bootstrap()+fastnose,
customvce() (unreachable under permute either way).

## FINDING PX-01 (self-found, confirmed) — see FINDINGS.md
fastnose + permanalyticse under permute(): permanalyticse is a SILENT NO-OP; the run posts
e(perm_analyticse)=1 and prints a header and footer asserting analytic SEs and
e(V)/e(V_es_*) that are not there. ado:5396 already computes the right predicate (`fnperm`)
and uses it to refuse fastnose+jointtests/adversarial in exactly this case.
Repro: PX01_repro.do (self-contained, generates its own panel).

## FINDING PX-03 (self-found, confirmed) — see FINDINGS.md
semielasticity on a non-strictly-positive outcome fails inside margins and is never
diagnosed. Four modes; the worst (rc 0, printed table, NO e(b)) fires on the fastnose path,
which is exactly what permute() selects by default. +10 shift of the same outcome fixes all.
Repro: PX03_repro.do (all four modes plus both controls).

## FINDING PX-02 (seeded lead, contradicted) — customvce()+permute() is unreachable
Repro: PX02_repro.do.

## Probes still in flight at the time of writing
p14_pconv.do (p-value conventions + p/CI coherence), p06_oracle.do (re-run after fixing
e() clobbering; the gteffects arm already matched the hand enumeration exactly),
p15_item11.do (known_bydesign s7 item 11: does fastnose post e(b_es_att)/e(numpre_att)),
PX02_repro.do. chain7.out / chain8.out track them.

## RESULT p15 — known_bydesign s7 item 11 SETTLED (no divergence)
fastnose posts e(b_es_att) (matrix) and e(numpre_att) (scalar 2) WITHOUT e(V_es_att),
under the documented caution; impvalues posts none of the three and is already off the
sthlp:1325 list on the pinned build. So the list is correct as it stands and removing
fastnose from it would have made the help wrong. See FINDINGS.md.

## RESULT p08 — safety() does NOT shape the permutation null on the design tried
200 reps, 45% late entry / 25% early exit, e(perm_failed)=0 under safety(on), safety(off),
emptycellsoverride and with an ES window; p and SD bit-identical. Negative result.

## RESULT p09/p18 — replay, state leakage, honestdid: CLEAN
Replay after a permute run reproduces e(b)/e(perm_p)/e(perm_ci)/e(perm_kernel) bitwise
(p09's apparent mismatch was MY bug: a local's decimal round-trip; redone in scalars in
p18_replay.do). honestdid refuses r(198) in fast mode and runs on the analytic VCV after
permanalyticse, which the footer discloses. No perm state leaks into the next plain run.

## RESULT p06 — exact-enumeration oracle
M=60 enumerated by hand; the gteffects arm matched e(permdist_b) with MAXRELDIFF exactly 0
(60 vs 60 draws). The controls()/timetrends arms initially read a clobbered e(); patched to
snapshot e(permdist_b) into PDSNAP before the oracle loop and requeued as chain8.

## RESULT p14 — p-value conventions MEASURED CLEAN
15 runs x every column (108 column-checks): MISMATCH=0, COHERENCE_VIOL=0, worst abs diff
0.000e+00. All five conventions, both exact and sampled branches, on b and t scales, and
an 18-column subgroup run. e(perm_p_lo)/e(perm_p_hi) match the hand shares exactly.

## RESULT p06 — oracle COMPLETE
plain / controls(zbig) / timetrends / gteffects all 60 vs 60 with MAXRELDIFF exactly 0;
negative control (controls run vs no-controls oracle) separates by 0.946.

## STATE: report delivered. FINDINGS.md is final.
