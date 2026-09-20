# ereturn-contract prober — NOTES

Target: the e() contract (names, shapes, PSD, 32-char limit, staleness).
Method: CONTRACT TESTING against `recon/ereturn_map.md` + automated assertions over a broad
configuration sweep.

## Files
- `mkdata.doi`  — D1 unit panel (120x12, 4 cohorts + never-treated), D2 state-year admin panel
  (50 states x 16 yrs, unbalanced + missing y), D3 repeated cross-section with a genuine
  continuous dose that varies in EVERY period (the F5 discriminator), D4 = D1 + subgroup
  variables designed as collision bait (values {1,10,11}, negative, fractional, FIPS-like).
- `ckc.doi`    — `ckc <tag>`: ~20 automated assertions on the live e(). Writes
  `VIO|tag|CODE|detail` / `INF|tag|CODE|detail` to $CKCOUT.
  Assertions: conformability, symmetry, PSD (min eigenvalue), rank, duplicate e(b) colnames
  (the silent-padding signature), colsof==#colnames for EVERY e() matrix, name length <= 32,
  se>=0, p in [0,1], lb<=est<=ub, t == est/se, CI halfwidth == se*t(df,level),
  relyr-vs-name-token, and FOUR cross-export agreements:
  e(b) vs wdidmainresults col1; diag(e(V)) vs se^2; e(b_es_<fam>) vs the ES_ cells of e(b);
  diag(e(V_es_<fam>)) vs the matching e(V) diagonal; wdidsgresults_<v> vs the _SG<v>_ rows.
- `estate.doi` — one-line post-call e() state dump + user-object survival check.
- `p01_smoke.do` / `p02_stale.do` / `p03_sweepA.do` / `p04_sweepB.do`
- `launch.sh`  — cap-respecting launcher (waits for <8 stata-mp, kills at a deadline).

## Status log
- 02:1x  harness written. Heavy contention: 8-9 other agents' Stata jobs, my launcher waits.
- READ SO FAR (code, not yet run):
  * `wooldid_validate`'s subgroup screens (ado 12217-12240) ERROR 411 on negative and on
    non-integer subgroup values reaching the treated estimation sample. So the "negative /
    fractional subgroup value breaks `^(.+)_s([0-9]+)$`" collision route is GUARDED. The
    regexes that would break are at ado:10906 (esband), 11057 (vsplice), 14637 (mata
    bootveskey) — all assume `[0-9]+` subgroup values.
  * `matrix drop _all` / `clear mata` appear only under the opt-in `cleanmatrices` (ado 591,
    4534) — the failed-run epilogue drops only its own named tempname estimates, so the
    "a failed wooldid destroys the user's own matrices" hypothesis is NOT supported by code
    reading (still tested empirically in p02).
  * bootsave's Mata writer truncates names to 32 chars (`substr(...,1,32)`, ado 14889/14901/
    14913/14921/14929); the w_/jointW_/advW loops guard with `_st_varindex(nm) < .` (silent
    SKIP of a colliding column), the `se_` loop does NOT (st_addvar would error). File
    surface, not e(); parked.

## Static results (no Stata needed) — completed while waiting for a slot
1. **documented-vs-code name diff**: extracted every `ereturn scalar|local|matrix` name from the
   pinned ado (322 sites, 274 distinct names -> `code_names.txt`) and every `e(...)` token from
   the sthlp (`doc_names.txt`). Every LITERAL e() name in the code appears in the help. Of the
   constructed names, the only VISIBLE ones whose literal prefix is absent are
   `roth_slope_ll_<type>` and `roth_slope_ul_<type>` — and **both ARE documented**, at
   `sthlp:1381`, written as `e(roth_slope_ll/ul_{it:type})`. **This CONTRADICTS recon F11**,
   which filed `e(roth_slope_ul_<type>)` as "appears nowhere in the help
   (`grep -c roth_slope_ul wooldid.sthlp` -> 0)": the grep misses it because the help writes the
   pair with a slash. Do not file F11's first bullet.
2. **recon F10's "9991-9994 have no strlen guard" is also WRONG.** `ado:9990` opens
   `if strlen("boot_pctl_sgresults_`sg'") <= 32 {` and the block CLOSES at 9996, so
   `boot_wald_`, `boot_pctl_`, `boot_bc_` and `boot_bca_sgresults_<sg>` are all inside one
   guard. The guard measures the *pctl* name (20 chars + len(v)); `boot_bc_` (18) and
   `boot_bca_` (19) are shorter, so they are over-skipped rather than crashed — conservative,
   and exactly what `sthlp:1293-1299` promises.
3. No e() name is posted with two different kinds (scalar vs matrix vs local), and none is
   posted both hidden and visible. No intra-package name-kind collision exists statically.
4. `matrix drop _all` / `clear mata` occur ONLY under the opt-in `cleanmatrices` (ado 591/4534),
   and `sthlp:1114` describes exactly that ("clear out all matrices stored in STATA's main
   memory and in mata"). No doc-vs-actual divergence there.
5. Subgroup values: `wooldid_validate` errors r(411) on NEGATIVE (ado 12218) and on NON-INTEGER
   (ado 12234) subgroup values reaching the treated estimation sample, which is what the
   `^(.+)_s([0-9]+)$` regexes at ado:10906 / 11057 / 14637 depend on. Both screens have a
   NARROWER domain than the regexes' consumers, though: the negative screen is
   `ifstatement & ttre < .` and the integer screen / glevelsof are `ifstatement & t >= ttre`.
   p05 probes the residue (a negative or fractional value living only on never-treated units).
6. `ado:4045-4047` confirms recon F3: `matrix colnames `eb' = `rownamespool'` has NO
   length guard, while the analogous `4110` block DOES (`colsof == sizeof`). p06_ragged hunts
   the geometries where one loop appends and another does not.

## 02:30-02:40 — first launch failed; two harness bugs found and fixed
1. **My own generator bug (not wooldid's).** `gen y = ... + 0.5*post*(t-ttre+1)` makes y MISSING
   on every never-treated row, because Stata's `0 * .` is `.`. wooldid then correctly reports
   "no never-treated control group", r(459). Fixed by computing the dynamic term under
   `if post==1`. Verified in `p00_diag.log`: the same panel shape without the bug runs rc 0
   (VA/VB/VB2/VC all 0).
2. **`foreach ... { cmd }` on one line is illegal in Stata** ("program error: code follows on the
   same line as open brace"). `nameslog` had three of them. Fixed in p03/p04/p06.
3. **Hardening after recon's own warning:** `RUN` now takes the run's rc and SKIPS every
   assertion when the run failed, because a failed wooldid leaves the PREVIOUS run's e() in
   place (F1) and a naive assertion pass would silently grade the previous config twice.

## More recon corrections found by reading the pinned build
- **recon §4b item 3 ("e(perm_p_<est>) is read out of a FIXED column `wdpc` at ado:211") is a
  misreading.** `ado:210` is `local wdpc = colnumb(`wdpk', "main_`wdpe'")` — a NAME lookup into
  `e(perm_kernel)`, guarded by `if `wdpc' < .`. A subset estimand set cannot mis-map by
  construction. p04's `PERMX` assertion checks it empirically anyway (scalar vs table vs
  `e(perm_p)` vector, for every estimand in `e(wdid_estimands)`).

## Free win: my assertions applied to recon's 81-config `edump.txt` (no Stata slot needed)
`edump_audit.txt`. Over the 73 configs whose e() was NOT the stale-reghdfe surface:
- **0 duplicate e(b) colnames** (the silent-padding signature, recon F3) in any config.
- **0 e(b) colnames over 32 characters.** The longest is exactly 32:
  `pretreatment_SG1000001_att_itime` on C05_sg7 — i.e. 7-digit subgroup values sit ON the
  boundary and 8-digit crash, exactly as the ruled ledger says.
- **e(b)/e(V) conformable in every config**; `VASYM = 0` in every config (perfect symmetry).
- The only negative/sentinel diagonals are the three `contreat()` configs C01/C02/C03
  (10 sentinel rows each, mineig = -1e100) — the ruled sentinel.
- Every `ES_<fam>[_SG<v>]_<horizon>` family occupies a CONTIGUOUS block of e(b) columns in
  every config (no interleaving that would break a `keep()`/positional handoff).
So the *static* half of the name/shape contract is clean across recon's whole surface; my own
sweep adds the value-level cross-export identities that edump could not see.

## Measured results ledger (see FINDINGS.md for the writeups)
- sweepA (A01-A28, 22 configs with a posted e()): **ZERO violations** of conformability,
  symmetry (VASYM=0 everywhere), PSD, rank, duplicate-colname, per-matrix shape-vs-namelist,
  se>=0, p in [0,1], CI containment, t==est/se, CI-halfwidth==se*t(df_r,level), relyr-vs-name,
  numpre-vs-pre-column-count, escomp-vs-b_es alignment, pwsgs-vs-e(b) differences, and all four
  cross-export identities. The only VIO was my own false positive on `esrelativeto(-2)` (fixed).
- p05_edge: 0 violations. `predict` after a SUCCESSFUL wooldid is r(111); `margins` r(111);
  `estat vce`, `estat summarize`, `estimates store`, `estimates table` all rc 0.
- p06_ragged: 27 violations, all explained -- 24 BMISS on the `amevoid` structural zeros (the
  RULED SI-A4/R6 design: e(b) carries 0 where the table carries missing), 1 VRANK on the same
  run (the -1e100 sentinel), and 2 VRANK on ANALYTIC runs which became EC-4.
- CONTREAT SENTINEL, F5's open question answered: on D3 (a repeated cross-section whose dose
  varies in EVERY period, i.e. a genuine intensity measure) `contreat(dose)
  contreatcontrols(untreatedZero)` gives SENT=0, ZERO=0, 0 missing cells in e(V), cholesky
  clean, and BOTH e(b_es_ame) and e(V_es_ame) posted. The -1e100 sentinel is specific to the
  degenerate `dose = dose0*post` DGP recon used; it does NOT fire on a realistic dose.

## Final state (04:30)
Findings EC-1 … EC-9 are written up in FINDINGS.md with self-contained repro do-files
`repro_EC1_bootrank.do`, `repro_EC2_level.do`, `repro_EC4_singular.do`,
`repro_EC5_stale.do`, `repro_EC7_sgnames.do`, `repro_EC8_restore.do`.
Probes that ran: p00 p02 p03 p05 p06 p07 p07b p09 p10 p11 p12 p13. Probes written but NOT run
(slot contention + deadlines): p04_sweepB (bootstrap/permute/contreat sweep + the PERMX
per-estimand p-value assertion), p03b_sweepA2 (A29-A35: oldsyntax, interactivecontrols,
att_itime alone, summarystats, nofastcollapse), p08_vce (the VCE bookkeeping grid, e(sample)
semantics vs the four sample concepts, the cfx 22-character name surface, svy abuse).

## Advisor-requested controls (p14_controls.do / .log), all three resolved
1. EC-2 mechanism NAILED: a Stata `real' macro-expands with %18.0g. `local x = 99.9' ->
   `99.90000000000001'; 90.1 -> `90.09999999999999'; 99.99 -> `99.98999999999999'; while
   68.3, 99.95, 92.5, 97.5, 99.5 all round-trip clean. The three that fail wooldid are exactly
   the three whose expansion shows more than two decimals, which `set level' refuses.
2. EC-4 DOWNGRADED from high to medium: the singularity is confined to cluster() with a
   unit-level `i'. On the design the help endorses (i = cohort + fe(id)) the same data and
   window give rank 8 of 8 and cholmiss 0 under cluster(), robust and iid alike. The help
   already flags that geometry in prose, so what survives is the marker defect: es_ndegen is
   absent in both cases and the help says its absence marks a clean covariance.
3. EC-6 SHARPENED: [pw=] SEs are identical to [aw=] robust SEs to nine decimals
   (0.210058487 both), so the estimator is correct and only e(vce)/e(vcetype)/e(wtype)/e(wexp)
   are wrong/absent. Cause named: the single `ereturn post' at ado:4051 carries no weight clause.
4. honestdid was NEVER successfully invoked (my syntax error both times: vcv() resolved to
   e(V)). No claim is made about it. `pretrends power' did run rc 0 on the singular block.

## p15_perm (permute), clean
With clusters correctly nested inside onset groups, permute(shuffle) permreps(40) runs rc 0 on
five different estimand SUBSETS (att; att+att_i; att_itime alone; att_i+att_itime; all four)
and on an event-study spec. In every one, e(perm_p_<est>) equals e(perm_p)[main_<est>] exactly.
recon sec 4b item 3's worry about a fixed column `wdpc' mis-mapping a subset is CLEAN -- and the
static reading confirms why (ado:210 is a colnumb() name lookup, not a fixed index).
e(permdist) vs e(permdist_b): 4 vs 4 columns when K=4, but 4 vs 16 on the event-study spec --
the documented deprecated/true-twin split (ado:241-249), not a defect.
ckc found ZERO violations on every permute configuration.
