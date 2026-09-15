# Old-wooldid bug attribution report - full assessment run, 2026-09-14

**Claim this report makes:** every old-vs-new difference observed in the 2026-09-14
full assessment run traces to a documented, verified bugfix in old wooldid, or to a documented
intended change ruled by the author. Each part below gives, per family: what old wooldid did
wrong and whom it affected; the mechanism located in the v1.3 source; the fix commit with its
verification evidence; the configs and ledger rows that demonstrate it in THIS run with observed
magnitudes; and a runnable two-build reproduction whose decisive lines are quoted.

Builds: **old** = published wooldid v1.3 (`git main`, `301b706`); **new** = the sha named on the
`prov[new]` line below, NOT whatever `git HEAD` happens to be when you read this. The two differ
as soon as any fix lands after the run, and a reader who checks out HEAD gets the wrong tree.

Run verdicts, read from the comparison this report documents:

    prov[new] = new|1a52a8e
    VERDICTS EXPLAINED=411 MATCH=26 RC_EXPLAINED=122
    COMPARED 559 of 559 configured
    SUMMARY unexplained=0 explained=4631 additions=17078 stale_ledger=0 unscoped=0 prov_fail=0

Every part was drafted by an independent agent working from the ledger, the dumps, the v1.3
source and `wooldid_solved_issues.md`, and each ran its own two-build reproduction. Where a
ledger row's wording and the record disagreed, or where a reproduction could not be made to
show the defect, the part says so under its own open-questions heading rather than smoothing it
over; those are collected in the final part.

---


# Part 1 - the pre-2026-08-23 bug families, re-verified at HEAD


> ## Headline of record and evidence provenance (added 2026-09-13; identical in all six sections)
>
> **Superseding note added 2026-09-14.** Two things below have since moved and are corrected
> here rather than rewritten in place, because the block is accurate about the run it describes.
> (i) `4280f3b` is no longer `HEAD`; the report's own header names the run of record by its
> `prov[new]` sha, and **that** is the build the report's claims attach to. (ii) The closing run
> of this assessment is the **2026-09-14** one over a **559-config** universe against the same
> v1.3 old side, reported in Part 7.5: **559 of 559 compared, zero DIFFERS, unexplained 0,
> stale_ledger 0**. The 2026-09-11 run described below is where this section's EVIDENCE was
> gathered and its magnitudes measured; the 2026-09-14 run is what verifies that nothing still
> differs. Neither supersedes the other - they answer different questions.
>
> **Run of record.** The 2026-09-11 dual-build run: **old = `301b706`** (published v1.3),
> **new = `14c30d3`**, whose `wooldid.ado` is byte-identical to this section's pin `08ec0b4`
> and to current `HEAD` apart from a nine-line comment block and a shortened `dis as text`
> permutation caution (10 insertions, 1 deletion; no estimation code).
>
> **The report judges those dumps against the ledger and comparer at `4280f3b` = `HEAD`.** That
> suite defines a **538-config universe**; the 2026-09-11 dumps cover **507** of them:
> **507 compared — EXPLAINED = 374, MATCH = 25, RC_EXPLAINED = 108 — unexplained 0,
> stale_ledger 0, unscoped 0, prov_fail 0, COVERAGE_OK.** The other 31 configs are new at
> `4280f3b`, have no 2026-09-11 dumps, and are outside every claim in this report.
>
> An earlier headline — **507/507 compared, EXPLAINED = 399, RC_EXPLAINED = 108**, at ledger
> `26ce421` — is quoted by three of the six drafts (`legacy` and `intended` give the full
> 399/108 line; `inference` gives the 507/507). It is the **same dumps and the same 507 configs**
> judged by the pre-`4280f3b` comparer. The two reconcile exactly: `399 = 374 + 25`, the 25 being
> the configs on which MATCH became reachable once the comparer began comparing across batch 3's
> macro→scalar kind change (`4280f3b`, verbatim: *"It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)"*; the same paragraph
> puts the cells affected at *"361 values on 103 configs"*). **Both numbers are measured and
> neither is withdrawn.** See `LEDGER_DELTA.md` §0.
>
> **Evidence provenance.** Every `test/suite/out/run_old|run_new/<config>.csv` cited below is a
> **regenerable** artifact, not a stored one: `test/suite/out/` is gitignored (`test/.gitignore`
> line 5) and the 2026-09-11 dumps were overwritten in place by the 2026-09-13 run. Regenerate
> them with `./run_full.sh oldnew` from `test/suite`, which stages old from
> `git show main:wooldidinstall/wooldid.ado` and new from `git show HEAD:…`. `compare_verbose_0911.txt`
> — cited by four sections, twice as `scratchpad/compare_verbose_0911.txt`, which is **not**
> repo-relative and never was — is likewise a recipe: rebuild it with `compare_oldnew.py … --verbose`,
> passing `--ledger`/`--configs` from whichever commit's judgement you want. The command, the build
> table and the one harness caveat (`4280f3b` changed `mkrun.py`'s `dumpres`, so a regenerated dump
> may differ in its `files|` line or gain a `matnames|…|UNREADABLE` line — **value cells are
> unaffected**) are in **`review_process/assessment_2026-09-13/repro/README.md`**.
>
> **Reproductions.** The do-files this section cites "in this directory" are now parameterised —
> repo root and scratch workdir as positional arguments — and live at
> **`review_process/assessment_2026-09-13/repro/<section>/`**. Their exact invocations are in each
> file's header and in that README. Bulk artifacts (dumps, logs, copied `.ado` builds) are not
> copied: they are regenerable, and the README says how.
>
> **Per-edit log:** `EDITS.md`. **Ledger delta and closed open questions:** `LEDGER_DELTA.md`.


**Scope.** The old-wooldid bug families whose proofs already exist in
`wooldid_old_bug_report.md` sections 1.1, 1.3, 1.4, 1.5, 1.6, 1.7 and 1.17, plus the
re-check of section 2 (the eight "not visible as old-vs-new mismatches" items). A1
(B1/B2), A5-cfx (B7), A6-ext (B9), A7-drop / A7-kit (B10), A12-ttb (tcen), A15-mint
(mint staleness), G10B-esf (esf support-first). TM2 and DDG are handled by other agents.

**Builds.** old = published v1.3, git `main` = `301b706` (3,450 lines); new = git HEAD =
`08ec0b4` (16,475 lines). Both extracted with `git show <sha>:wooldidinstall/wooldid.ado`
into `old/` and `new/` of this directory and run in separate Stata processes.

**Run of record.** The 2026-09-11 dual-build run (dumps in `test/suite/out/run_old` and
`run_new`), re-judged verbosely in `compare_verbose_0911.txt` (**at ledger `26ce421`**;
the HEAD-ledger judgement is in the headline block above and in `LEDGER_DELTA.md` §0): **507/507 compared,
EXPLAINED = 399, RC_EXPLAINED = 108 (66 designed rc changes + 42 both-error parity
pins), 5,121 explained difference groups, 16,024 informational new-only additions, 0
unexplained, 0 stale ledger rows, 0 unscoped, 0 provenance failures.**

**Verification pass 2 (2026-09-13).** A second agent instance re-verified this report
from disk without re-running Stata: every HEAD line number below was re-located in
`new/wooldid.ado`, every `git log -S/-G` claim re-run, every cited cert suite confirmed
present at `08ec0b4`, the four reproduction logs confirmed complete (sentinels, rc 0,
`which wooldid` resolving to `old/` and `new/`), and the reproduction dumps re-tied to
`test/suite/out/run_{old,new}` with `test/compare.py --rtol 0` (basic 177/177 old,
278/278 new; trig_tt_big 96/96 old, 106/106 new numeric cells identical — the only
"difference groups" are `matnames` rows the reproduction's dumpres did not emit). Four
corrections landed: the L.4 cell labels (main_att AND main_att_it are bit-identical;
the movers are the att_i/att_itime cells), the HEAD-side mechanism of the L.5 `e(V)`
drop, the L.8 item-6 wording (the new-side T-M2 value is posted as a *scalar*, the
macro is removed), and the `con_p75` probe's second hit. Pass 2's `> NOTES.md`
overwrote pass 1's running notes; the resume record is in NOTES.md.

**What "still holds" means here — three checks per family, all performed:**

1. *The fix is still in the code at HEAD* — located by reading `new/wooldid.ado`
   (line numbers below are HEAD `08ec0b4`; v1.3 numbers were re-located in
   `old/wooldid.ado`, never trusted from the older report).
2. *Nothing in the 2026-08-24..09-11 commits changed the class* —
   `git log --since=2026-08-24 -S<keyword> -- wooldidinstall/wooldid.ado` is **empty**
   for every fix-line keyword, and each keyword is proven sensitive by its positive
   control (the same `-S` without `--since` returns exactly the fix commit):
   `e(sample) & \`ttre' == .` → `015ae88` (B1); `con_p75` → `0188b8f` (B2; its only other hit is `fc17051`, the 2023 original upload that introduced the local);
   `cfxseopt` → `3df74db` (B7); `\`t' >= . | \`y' >= .` → `2adf447` (B9);
   `total(\`esamp' * \`w')` → `40e1913` (B10; the two later hits `40414c7`, `d87adfd`
   are pre-08-24 line-count touches); `tcen` → `f7019c9`; `nopre` → `0c03010`;
   `VAL-G10B: the calendar-support screen runs FIRST` → `7ba1e16` (`controlcheck2`
   itself is NOT a sensitive keyword for G10B — the reorder moved the block without
   changing its occurrence count — so the fix commit's own comment text is the probe).
   The only post-08-24 commits that touched these blocks are
   two byte-identical relocations into subprograms — `308472c` (2026-09-08:
   `wooldid_reltime`, `wooldid_ccscreen`, `wooldid_sgdom`, `wooldid_hist`) and `6f604b9`
   (2026-09-08: `wooldid_sumstat`) — each certified in its commit message by a
   twenty-leg branch fingerprint and a zero-difference battery — and, for the cfx
   block only, `82934b3` (D-1, 2026-09-11), which changed the `cfxmatendpoint` line's
   permute() carve-out and nothing else (section L.2).
3. *The ledger rows fire on the same configs with the same magnitudes* — read from
   `compare_verbose_0911.txt`; every headline magnitude of the 2026-08-23 report
   recurs unchanged (2.92e-01 basic, 1.68e-01 adm_base, 3.33e-01 cont_ss, 5.07e-02 /
   1.36e-01 trig_extmiss, 8.22e-02 / 2.50e-01 trig_allest_drop, 2.26e-01 trig_tt_big,
   N 296 → 290 trig_mintsent_off, att 0.6149 → 1.4106 trig_esfsup_ref). Given check 2
   that is the expected outcome, and it is recorded here as a verified check rather
   than assumed.

**Runnable reproductions (this run).** Two families were re-run on both builds, as the
task asked: `repro_A1A2.do` (B1/B2 on config `basic`) and `repro_A12TTB.do` (tcen on
config `trig_tt_big`). `run_repro.sh` runs each do-file on the old build and then on the
new build in **separate Stata processes, strictly one at a time** (it waits while eight
or more `stata-mp` processes are running and kills a leg after 300 s); the do-files
mirror `mkrun.py` exactly (`set seed 98765`, `set sortseed 98765`, `ereturn clear`, the
config's own dataset and option string) and load the build by `adopath ++`. Every
wooldid call is `cap noi` with an `_RC=` line, every value is printed in `%21x`, and the
sentinels `REPRO_A1A2_DONE` / `REPRO_A12TTB_DONE` close each leg. `repro_cross.py`
decodes the logs and ties them to this run's dumps. All four legs completed rc 0 in
6-8 s each on 2026-09-12 22:32-22:33; the old-leg and new-leg dumps agree with
`run_old/*.csv` and `run_new/*.csv` **byte-for-byte on every common scalar/matrix
cell** (basic: 0/177 old, 0/278 new differ at rtol 0; trig_tt_big: 0/96 old, 0/106 new),
so the reproductions are the run, not a look-alike of it.

---

### L.1 Summary statistics: "control group" rows were pre-treatment-treated rows (B1) and sub-sample IQRs were the full-sample IQR (B2) — ledger `A1-B1a..e`, `A2-B2a..e`

- **Symptom / impact:** every `y_control`/`c_control` row of `e(summarystats)`
  duplicated the `*_pretreat` row, and every pre/post/control IQR cell plus the derived
  `e(*_iqr)` scalars carried the full-sample IQR. Users citing these descriptive tables
  published wrong control-group and dispersion statistics. Estimation was never touched.
- **Mechanism (v1.3, re-located):** the control-group `if` at old line 1379,
  `if e(sample) & \`ttre' <. & \`t' < \`ttre'`, is character-identical to the
  pre-treatment block's condition at old line 1353 (never-treated units are `ttre == .`);
  the three IQR lines 1359 / 1372 / 1385 read the in-sample locals
  (`\`\`sumvar'_p75' - \`\`sumvar'_p25'`) instead of the `pre_`/`post_`/`con_` ones, and
  line 1359's value is what `e(ypre_iqr)` returns.
- **Fix + proof:** `015ae88` (B1) and `0188b8f` (B2), 2026-07-19 — verified by 3-assessor
  + 3-verifier waves, only the affected cells changed, estimation bit-identical
  (`wooldid_solved_issues.md`, B1/B2). **At HEAD** the block lives in the subprogram
  `wooldid_sumstat` (moved verbatim by `6f604b9`): the control-group condition is
  `if e(sample) & \`ttre' == .` (new line 12024, with the comment "units never treated
  in-window"), and the IQRs read `pre_p75 - pre_p25` / `post_` / `con_` at new lines
  12004, 12017, 12030.
- **This run (2026-09-11):** `A1-B1a`/`A2-B2a` on `basic` — `summarystats` 10 cells,
  worst reldiff **2.92e-01** ([4,1] control mean 2.1024 → 2.2825; [2,4] pre IQR 1.5684 →
  1.3347; [3,4] post IQR 1.5684 → 1.3413), `ypre_iqr` 1.5684 → 1.3347. `A1-B1b`/`A2-B2b`
  on `pw_ss` (10 cells, 2.92e-01; `ypre_iqr` 1.4997 → 1.2985). `A1-B1c`/`A2-B2c` on
  `cont_ss` (20 cells incl. the contreat rows, 3.33e-01; `ypre_iqr` 11.081 → 7.669).
  `A1-B1d`/`A2-B2d` on `adm_base` (10 cells, 1.68e-01; 1.5179 → 1.2635). `A1-B1e`/`A2-B2e`
  on `ext_bin_ext_support` (5.88e-01; 0.7767 → 0.3787), `ext_cont_ext_support` (11 cells,
  5.46e-01; 1.8241 → 0.9055) and `ext_contract10_summary_control` (5.88e-01). Also
  bundled into `cov_kitchen` under `A7-kit-m` (summarystats worst 1.48e-01; `ypre_iqr`
  3.7133 → 3.7775). Identical to the 2026-08-23 magnitudes on every config the earlier
  report quoted. Co-firing rows on these configs (`CORE10-cmdline` — **retired at `4280f3b`**,
  `LEDGER_DELTA.md` §A — `PS62-esblock`,
  `AMESE-R1-*`, `X-est-*`) belong to other families.
- **Runnable reproduction — `repro_A1A2.do`, config `basic` on `test/data/dta_binary.dta`.**
  The hand truth is computed on the SAME realized `e(sample)` with `gstats sum` (the
  user `ttre` = 2005/2008 never exceeds max `t` = 2011, so the user variable equals the
  internal onset). Decisive lines (`repro_cross.out`):

  ```
  [old] control_mean=2.102447086016719 pretreat_mean=2.102447086016719 TRUTH_control_mean=2.282459182105959
  [old]   control==pretreat: YES (reldif 0.0e+00);  control==TRUTH: NO (reldif 7.9e-02)
  [old] insample_iqr=1.568425595760346 pre_iqr=1.568425595760346 post_iqr=1.568425595760346 control_iqr=1.568425595760346
  [old]   TRUTH pre_iqr=1.334739089012146 post_iqr=1.341322660446167 control_iqr=1.372958660125732
  [old]   pre_iqr==insample_iqr: YES;  pre_iqr==TRUTH: NO (reldif 1.5e-01);  post_iqr==TRUTH: NO (1.4e-01);  control_iqr==TRUTH: NO (1.2e-01)
  [new] control_mean=2.282459182105959 pretreat_mean=2.102447086016719 TRUTH_control_mean=2.282459182105959
  [new]   control==pretreat: NO (reldif 7.9e-02);  control==TRUTH: YES (reldif 1.9e-16)
  [new] insample_iqr=1.568425595760346 pre_iqr=1.334739089012146 post_iqr=1.341322660446167 control_iqr=1.372958660125733
  [new]   pre_iqr==TRUTH: YES (0.0e+00);  post_iqr==TRUTH: YES (0.0e+00);  control_iqr==TRUTH: YES (4.9e-16);  e(ypre_iqr)==SS[2,4]: YES
  byte-tie [old] summarystats[4,1]: repro +1.0d1cfc7201934X+001 vs 0911 dump +1.0d1cfc7201934X+001 -> SAME
  byte-tie [new] summarystats[4,1]: repro +1.24279f5acccceX+001 vs 0911 dump +1.24279f5acccceX+001 -> SAME
  dump-tie [old] basic: 177 cells in common with the 0911 dump, 0 differ at rtol 0
  dump-tie [new] basic: 278 cells in common with the 0911 dump, 0 differ at rtol 0
  ```
  Old: the control row IS the pre-treatment row (bit-identical), and all four IQR cells
  are the in-sample IQR; new: every cell equals the hand truth (the 1.9e-16 / 4.9e-16
  residues on the control cells, and the 2-ulp in-sample IQR residue that appears on
  BOTH builds, are machine-precision residues whose mechanism was not diagnosed — the
  cells pass through locals before entering the matrix — and sit five orders of
  magnitude inside the 1e-10 indifference threshold).
- **Classification: old-bug-fixed.**

### L.2 cfx matrices truncated to width 1 on every estimand after the first (B7) — ledger `A5-cfx-m`, `A5-cfx-n`

- **Symptom / impact:** with `makecfxplots` + multiple estimands, old overwrote the
  `cfxplotnose` option flag with the `marginsdof` contents mid-loop, so every estimand
  after the first returned cfx matrices truncated to the estimate column — SEs/CIs
  silently lost — and plot titles carried a tempvar name. The 2026-07 assessment proved
  the truncation cluster-independent and positional.
- **Mechanism (v1.3, re-located):** old lines 1869-1871 —
  `if "\`cfxplotnose'" != "" local cfxmatendpoint 1` /
  `if "\`cfxplotnose'" == "" local cfxplotnose \`marginsdof'` /
  `if "\`cfxplotnose'" == "cfxplotnose" local cfxplotnose nose` — the second line
  makes the first line true on the next pass of the estimand loop.
- **Fix + proof:** `3df74db` (2026-07-20), per-iteration copies; verified 3/3 with
  estimand-2+ matrices restored full-width and bit-identical to single-estimand
  references (`wooldid_solved_issues.md`, B7). **At HEAD** lines 2338-2339:
  `local cfxseopt \`marginsdof'` / `if "\`cfxplotnose'" == "cfxplotnose" local cfxseopt nose`,
  with the user flag never modified. Two later commits touched the *neighbouring*
  `cfxmatendpoint` line — `1561386` (2026-09-05, fastnose/permute) and `82934b3` (D-1,
  2026-09-11): it now reads `if "\`cfxplotnose'" != "" & "\`permute'" == "" local
  cfxmatendpoint 1`, a permute()-only carve-out. `trig_cfx` has neither `cfxplotnose`
  nor `permute()`, so the width stays 6 and the class is unchanged.
- **This run (2026-09-11), `trig_cfx`:** `A5-cfx-n` — four matnames rows,
  `cfxatt_it_cfxpost`, `cfxatt_it_relyr2`, `cfxame_it_cfxpost`, `cfxame_it_relyr2`:
  `3x1 |b` on old → `3x6 |b se t pvalue ll ul` on new (the stripe rename to
  `lattice_#` is NFFAST-02, row `A6-cfxstripe`). `A5-cfx-m` — `cfxatt_it_cfxpost` and
  `cfxatt_it_relyr2`: 15 cells each, all `ABSENT -> .`/value (the five inference
  columns that old never returned; the att point cells [r,1] are unchanged);
  `cfxame_it_cfxpost` and `cfxame_it_relyr2`: 18 cells each (the 15 absent inference
  cells plus the 3 point cells, see the caveat). The B7 signature is exactly as on
  2026-08-23.
- **Attribution caveat (open question, not papered over — **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E1]**; the
  `cfxame*`-scoped row it asks for is now `SI20-L-trig_cfx-cfxame` and `A5-cfx-m` is narrowed
  to `cfx*_it_*`):** `A5-cfx-m` this run ALSO
  lists `cfxame_cfxpost` and `cfxame_relyr2` — 3x6 on both builds — with 3 value cells
  each, worst reldiff 6.70e-01 ([1,1] 3.0091 → 1.0071), and the [r,1] point cells of the
  `cfxame_it_*` matrices move likewise (3.0091 → 1.0071; 3.1533 → 1.1513). Those are not
  B7: the att point cells do not move and old's 3x6 first-estimand matrices were never
  truncated. They are the SI-20 class (`1df8613`, 2026-09-06; ledger row
  `SI20-L-trig_cfx-b`: "lattice + contreatcontrols(basic) now absorbs the baseline dose
  control ... ame-family cells move by the baseline dose slope (O(1))"), which on this
  config also moves `e(b)`'s ame cells 3.0026 → 0.9980. The ledger judges first-match
  and `A5-cfx-m`'s `cfx*` pattern (ledger line 70) precedes the SI-20 rows (lines
  313-315), so the SI-20 movement on the ame cfx matrices is being explained by a row
  whose reason text is B7-only — ledger reason and solved-issues record disagree for
  those nine cells. The repo is read-only for this agent; the recommended repair (per
  the project rule "name a new cell class in the comparer, never blanket the matrix")
  is a `cfxame*`-scoped SI-20 row for `trig_cfx` placed ahead of `A5-cfx-m`, or an
  amendment of `A5-cfx-m`'s reason naming both classes. The B7 part is unaffected.
- **Reproduction:** not re-run this pass (the task's two reproductions are L.1 and
L.5). **No two-build do-file for this family** (see L.9.3 for the full coverage statement). What
  stands in its place is a **re-tie of the decisive cells to a dump that exists today**: the
  2026-09-13 dual-build run in `test/suite/out/run_old|run_new` is on the *same pinned pair* this
  section is about — `prov|old|301b706` and `prov|new|6a05127`, whose `wooldid.ado` is
  byte-identical to the pin `08ec0b4` (md5 `382c1e0ee3e1be37fa2e9b677a65fdcc`). Read 2026-09-14;
  artifact `critleft/retie_0914.txt`. On `trig_cfx`:
  ```
  matnames cfxatt_it_cfxpost  old 3x1|…|b                         new 3x6|lattice_1 lattice_2 lattice_3|b se t pvalue ll ul
  matnames cfxame_cfxpost     old 3x6|intensity:1bn.__000009 …|…  new 3x6|lattice_1 lattice_2 lattice_3|…
  e(b) main_ame               old 3.002574284                     new 0.9980211786
  cfxame_cfxpost[1,1]         old 3.009137279                     new 1.007105902
  ```
  The B7 width signature (`3x1|b` → `3x6|b se t pvalue ll ul` on the `_it_` matrices) and the
  SI-20 ame movement (3.0026 → 0.9980, 3.0091 → 1.0071) both reproduce to the digits quoted above.
- **Classification: old-bug-fixed** (the width-1 truncation). The co-fired ame value
  cells are a different, later, old-bug-fixed family (SI-20) misattributed by row order.

### L.3 Extended missing values (`.a`–`.z`) escaped the missingness screens (B9) — ledger `A6-ext-m`

- **Symptom / impact:** old tested `== .`, so `.a`-style missings — routine in survey
  data — stayed in the screened sample: they inflated the `w_per_it`/`w_per_i`/`t_per_i`
  weight denominators (diluting real observations' cell weights), and `t = .a` on a
  treated unit polluted the post-treatment classification (`t >= ttre` is true).
- **Mechanism (v1.3, re-located):** old line 421
  `qui replace \`ifvar' = 0 if \`i' == . | \`t' == . | \`y' == . | \`w' == . | \`w' <= 0`
  and old line 443 `qui replace \`ifvar' = 0 if \`mustbenonmissingvar' == .`.
- **Fix + proof:** `2adf447` (2026-07-20), `>= .` at both sites; verified 3+3 —
  corrupted-data runs bit-match cleaned-reference runs at every injection point,
  estimand, weight and estimator; clean data provably bit-identical
  (`wooldid_solved_issues.md`, B9). **At HEAD** new lines 674 and 696 carry `>= .`
  at both sites; the same tolerance is now used by the F13 onset screens (755-756).
- **This run (2026-09-11), `trig_extmiss`** (accepted blind-spot config, `dta_s_extmiss`
  with extended missings in y/t/w): `A6-ext-m` — `b` 2 cells, worst **5.07e-02**
  ([1,3] 0.70571 → 0.70423; [1,6] −0.08801 → −0.08355), `V` 2 cells 1.94e-03,
  `wdidmainresults` 12 cells 1.36e-01. Identical to 2026-08-23. Old's weights were
  diluted by the phantom rows; new matches the cleaned-data answer. Co-firing
  `PS62-esblock` (24 cells) is the PS-62 family.
- **Reproduction:** not re-run this pass. **No two-build do-file for this family** (see L.9.3
  for the full coverage statement). What stands in its place is a **re-tie of the decisive cells
  to a dump that exists today**: the 2026-09-13 dual-build run in `test/suite/out/run_old|run_new`
  is on the *same pinned pair* this section is about — `prov|old|301b706` and `prov|new|6a05127`,
  whose `wooldid.ado` is byte-identical to the pin `08ec0b4` (md5
  `382c1e0ee3e1be37fa2e9b677a65fdcc`). Read 2026-09-14; artifact `critleft/retie_0914.txt`. On
  `trig_extmiss`:
  ```
  e(b)[1,3]   old 0.705705705      new 0.7042335867
  e(b)[1,6]   old -0.08800843031   new -0.0835492375
  ```
  Both reproduce the magnitudes quoted above (worst 5.07e-02 on `b`).
- **Classification: old-bug-fixed.**

### L.4 Estimand weights built on the pre-regression sample, not `e(sample)` (B10) — ledger `A7-drop-m`, `A7-kit-s/m`

- **Symptom / impact:** reghdfe singleton drops (and B9-class rows) made the regression
  sample strictly smaller than the screened sample, but the att_it/att_i/att_itime cell
  weights normalized over the larger set — cells containing dropped observations were
  silently under-weighted, and old's att_it collapsed exactly onto plain att whenever a
  whole cell dropped.
- **Mechanism (v1.3, re-located):** the weight inputs at old lines 1063-1070
  (`egen \`w_per_it' = total(\`ifvar' * \`w')` and the `w_per_i`/`t_per_i` twins, both
  code paths) run BEFORE estimation; `e(sample)` is only materialized at old lines
  1168-1169 (`gen byte \`esamp' = e(sample)`).
- **Fix + proof:** `40e1913` (2026-07-20): the block moved to immediately after
  `esamp` exists, with `esamp` substituted in the three normalizers, the `tag()` cell
  counter and the `posttreatment` classifier; verified 3+3 — under singleton/separation
  drops full-data runs match cleaned-reference runs bit-for-bit, plain att provably
  invariant, bit-identical on no-drop data across 18+ configs
  (`wooldid_solved_issues.md`, B10). **At HEAD** `gen byte \`esamp' = e(sample)` is at
  new line 1430 and the normalizers at 1538-1540 read `total(\`esamp' * \`w')`,
  `total(\`w'*\`esamp')`, `total(\`tag')` with `tag`/`posttreatment` built on
  `\`esamp' == 1` (1520-1529).
- **This run (2026-09-11):** `A7-drop-m` on `trig_allest_drop` (blind-spot config with
  forced fe-singleton drops): `b` 4 cells, worst **8.22e-02**. Named by the dump's
  `matnames|b` row (`main_att main_att_it main_att_i main_att_itime pretreatment_att
  pretreatment_att_it pretreatment_att_i pretreatment_att_itime ES_att_* ES_att_it_*`):
  [1,1] `main_att` and [1,2] `main_att_it` are **bit-identical** across builds
  (`+1.12b0e168eb563X+001` on both); [1,3] `main_att_i` and [1,4] `main_att_itime` move
  1.8e-09 relative (2.2806937434 → 2.2806937476); [1,7] `pretreatment_att_i` and [1,8]
  `pretreatment_att_itime` move 8.22e-02 (0.28393 → 0.30935). That is the B10 signature
  on this design, read in the honest direction of inference: att_it's bit-identity
  IMPLIES the dropped observations shared no (i,t) cell with survivors (so surviving
  cells' `w_per_it` was unchanged on either sample base), while att_i/att_itime move
  because their per-unit normalizers `w_per_i`/`t_per_i` are what the fix moved onto
  `e(sample)`. (The cell structure of the drops was inferred from the dump, not
  observed directly.) `V` 4 cells 1.54e-02; `wdidmainresults` 24 cells 2.50e-01. `A7-kit-s/m` on
  `cov_kitchen` (controls with missings + B1/B2): `b` 11 cells worst 1.41e-07, `V`
  2.11e-09, `wdidmainresults` 62 cells 1.42e-07, `summarystats` 1.48e-01, `joint_horiz`
  0.11412498 → 0.11412499, `joint_horiz_i` 0.23377417 → 0.23377418, `ypre_iqr` 3.7133 →
  3.7775. All identical to 2026-08-23.
- **Reproduction:** not re-run this pass. **No two-build do-file for this family** (see L.9.3
  for the full coverage statement). What stands in its place is a **re-tie of the decisive cells
  to a dump that exists today**: the 2026-09-13 dual-build run in `test/suite/out/run_old|run_new`
  is on the *same pinned pair* this section is about — `prov|old|301b706` and `prov|new|6a05127`,
  whose `wooldid.ado` is byte-identical to the pin `08ec0b4` (md5
  `382c1e0ee3e1be37fa2e9b677a65fdcc`). Read 2026-09-14; artifact `critleft/retie_0914.txt`. On
  `trig_allest_drop` and
  `cov_kitchen`:
  ```
  trig_allest_drop  e(b) main_att            old 2.146022965   new 2.146022965   BIT-IDENTICAL +1.12b0e168eb563X+001
                    e(b) main_att_it         old 2.146022965   new 2.146022965   BIT-IDENTICAL +1.12b0e168eb563X+001
                    e(b) main_att_i          old 2.280693743   new 2.280693748
                    e(b) main_att_itime      old 2.280693743   new 2.280693748
                    e(b) pretreatment_att_i  old 0.2839344819  new 0.3093541662
                    e(b) pretreatment_att_itime old 0.2839344819 new 0.3093541662
  cov_kitchen       joint_horiz              old 0.1141249833  new 0.1141249916
                    joint_horiz_i            old 0.2337741748  new 0.2337741842
                    ypre_iqr                 old 3.713320289   new 3.777455563
  ```
  The att/att_it bit-identity that carries this item's inference — and the `%21x` string
  `+1.12b0e168eb563X+001` it rests on — is re-read directly out of both dumps.
- **Classification: old-bug-fixed.**

### L.5 `timetrends` numerically inert at large calendar codes (tcen) — ledger `A12-ttb-s/m/p/vn`, `A12-ttp-s`

- **Symptom / impact:** with raw large time codes (measured as low as t ~ 5e5) the
  demeaned `i.i#c.t` trend columns lost all precision and the requested unit trends went
  numerically inert — rc 0, no caution, results exactly equal to the no-timetrends run.
  Users believed they had controlled for trends when they had not.
- **Mechanism (v1.3, re-located):** old line 939
  `if "\`timetrends'" == "timetrends" local trendinjector i.\`i'##c.\`t'` — `t` enters
  raw.
- **Fix + proof:** `f7019c9` (2026-07-27, tranche 13): a mean-centered `gen double`
  trend tempvar, exactly shift-invariant; suite F8/F9, battery 28/28, two verifiers
  (48/48 non-timetrends configs BYTE-identical; timetrends configs numeric-only within
  the pre-registered footprint) (`wooldid_solved_issues.md`, tranche 13 "Centered
  trend"). **At HEAD** new lines 1166-1171: `tempvar tcen` / `gstats sum \`t' if
  \`ifstatement', meanonly` / `gen double \`tcen' = \`t' - r(mean)` /
  `local trendinjector i.\`i'##c.\`tcen'`. Certified on HEAD by
  `test/test_ttrefloat.do` F8 (large-grid timetrends == its recode) and F9 (timetrends
  vs no-timetrends must differ materially at large t).
- **This run (2026-09-11), `trig_tt_big`** (blind-spot config, `dta_s_ttbig`: t =
  500000..500012, onsets 500006/500009): `A12-ttb-s` — ll −116.40 → −70.05, r2 0.8887 →
  0.9069, r2_within 0.7359 → 0.7790, rmse 0.4067 → 0.4001, rss 47.639 → 39.861;
  `A12-ttb-m` — `b` 6 cells worst **2.26e-01** ([1,1] 1.60167 → 1.38545),
  `wdidmainresults` 36 cells (the se/t/p/CI columns become missing on new), `V` REMOVED;
  `A12-ttb-p` `properties` [b V] → [b]; `A12-ttb-vn` the matnames side of the V drop.
  `A12-ttp-s` on `trig_tt_pois` bounds the poisson IRLS reparameterization noise (this
  run rss/rmse 2.06e-09, bound 1e-4, run-varying class). Identical to 2026-08-23.
- **Runnable reproduction — `repro_A12TTB.do`, config `trig_tt_big` on
  `test/suite/data/dta_s_ttbig.dta`.** Three legs per build: A = the suite config
  (timetrends, raw t); B = the same without `timetrends`; C = timetrends on the recoded
  panel `t - 499999` (t = 1..13, onsets 7/10; shifted by 499999, not 500000, so no
  period lands on t = 0), the small-code benchmark both builds fit correctly. Decisive
  lines (`repro_cross.out`):

  ```
  [old] att TT_RAW=1.601672573691604  NOTT_RAW=1.601672573691604  TT_RECODED=1.385452648276623  props TT_RAW=b V
  [old] rss TT_RAW=47.63904846392923  NOTT_RAW=47.63904846392921  TT_RECODED=39.8607067737382
  [old]   reldif(att TT_RAW, NOTT_RAW)=0.000e+00  reldif(rss TT_RAW, NOTT_RAW)=4.475e-16   [0 = trends INERT]
  [old]   reldif(att TT_RAW, TT_RECODED)=1.350e-01  reldif(rss TT_RAW, TT_RECODED)=1.633e-01
  [new] att TT_RAW=1.38545264827662  NOTT_RAW=1.601672573691604  TT_RECODED=1.38545264827662  props TT_RAW=b
  [new] rss TT_RAW=39.8607067737382  NOTT_RAW=47.63904846392921  TT_RECODED=39.8607067737382
  [new]   reldif(att TT_RAW, NOTT_RAW)=1.350e-01  reldif(rss TT_RAW, NOTT_RAW)=1.633e-01
  [new]   reldif(att TT_RAW, TT_RECODED)=0.000e+00  reldif(rss TT_RAW, TT_RECODED)=0.000e+00   [shift-invariant]
  byte-tie [old] b[1,1] (leg A): repro +1.9a0736bae7d4aX+000 vs 0911 dump +1.9a0736bae7d4aX+000 -> SAME
  byte-tie [new] b[1,1] (leg A): repro +1.62ad065681323X+000 vs 0911 dump +1.62ad065681323X+000 -> SAME
  dump-tie [old] trig_tt_big: 96 cells in common with the 0911 dump, 0 differ at rtol 0
  dump-tie [new] trig_tt_big: 106 cells in common with the 0911 dump, 0 differ at rtol 0
  [cross] reldif(att old TT_RECODED benchmark, new TT_RAW)=2.244e-15  rss: 0.000e+00
  [cross] reldif(att old TT_RAW, new TT_RAW)=1.350e-01   (the A12-ttb-m ledger movement)
  ```
  The three-value triangle: **old's timetrends run at t ~ 5e5 is bit-identical to its
  own no-timetrends run** (att reldif exactly 0; rss 4.5e-16) — the trends are inert, rc
  0, `e(V)` posted as if nothing were wrong; **old's own code on the recoded panel gives
  1.38545264828**, which is **new's answer on the raw large codes to 2.2e-15** (rss
  bit-identical); and new's raw and recoded runs are bit-identical (shift invariance).
  Two side observations: (i) old on the recoded panel also posts `b` only
  (`PROPS_TT_RECODED=b`), so the missing `e(V)` under *working* unit trends on this
  panel is a property of the design on both builds — the "documented signature" of rows
  `A12-ttb-p/vn`, not a new-build regression; (ii) the [1,1] reldif 1.35e-01 is smaller
  than the ledger's 2.26e-01 because the ledger's worst is over all six `b` cells.
- **HEAD-side mechanism of the `e(V)` drop (rows `A12-ttb-p/vn`), checked in pass 2:**
  with WORKING unit trends `i.i##c.tcen` and `cluster(i)`, every regressor is
  within-cluster, so the cluster-robust VCE is identically zero; HEAD's vzero screen
  (`wooldid_amescreen , mode(vzero)`, HEAD ~1484-1488) names it — the 2026-09-11 shard
  log for `trig_tt_big` on the new build carries "Caution: the estimation command
  returned an all-zero variance-covariance matrix, so no standard error in this run is
  identified", "Caution: standard errors were not identified for 6 of 6 reported
  estimates, so no variance-covariance matrix e(V) was posted" and "Caution: no
  horizon of the att event study has an identified variance" — and the b-only posting
  falls through the `ereturn post` fallback at HEAD 3898-3931. Old's inert trends escaped
  the zero VCE for the wrong reason (its trend columns were numerically garbage, not
  within-cluster exact), which is why old posts `b V` and the old shard log
  (`wd_old/shard_trig_00/shard_trig_00.log`, the whole block from the previous config's
  `dumpres` at line 5145 to `dumpres trig_tt_big` at 5213) carries no caution — the only
  regex hit for caution/error/warning there is the banner word "Errors" in "Standard
  Errors: cluster(i)"; the new-build block (24227-24312) carries eight. Commit `b801282` (2026-09-08, refit-rescue removal) sits on this
  posting rule but its message records "No timetrends config moved"; nothing in it
  changes the class.
- **Classification: old-bug-fixed.**

### L.6 Sentinel-coded onsets silently biased estimates under `safety(off)` (mint staleness) — ledger `A15-mint-s/m/e`

- **Symptom / impact:** with onset variables carrying sentinel codes (0/9999 for
  never-treated) and `safety(off)` — exactly what the F13 error message tells such users
  to specify — old's row-level no-pre-period drop could strand surviving units without
  their reference period; stale `mint`/`mintp1` then made a treated dummy silently
  drop: 18-29% biased att with fabricated pre-trends up to sign flips, rc 0, the alarm
  disabled by the same `safety(off)`.
- **Mechanism (v1.3, re-located):** `mint`/`mintp1` are computed at old lines 481-500,
  the drop at 506-507 is ROW-level (`replace \`ifvar' = 0 if \`ttre' <= \`mint'` /
  `... <= \`mintp1'` under timetrends), nothing recomputes `mint` afterwards, and the
  reference-forcing clauses at 520-546 (`replace \`_t' = . if ... \`t' == \`mint'`)
  then target a period the surviving unit may no longer have.
- **Fix + proof:** `0c03010` (2026-07-27, tranche 12): unit-level drop (`gegen max by
  i`) plus a post-drop recompute of `mint`/`mintp1` as an invariant guard;
  `test/test_mint.do` pins sentinel == unit-deleted benchmark (plain + timetrends);
  parent probe att 1.4234 vs truth 1.7436; verifier cm_v1 10/10 (partial contamination
  bench-exact where the parent was biased mreldif .29; user-if twin closed, parent
  mreldif .43) (`wooldid_solved_issues.md`, "mint staleness"). **At HEAD** the code is
  in `wooldid_reltime` (moved verbatim by `308472c`), new lines 11553-11575:
  `gegen byte \`nopre' = max(\`ttre' <= \`mint') if \`ifstatement', by(\`i')` (and the
  `mintp1` twin), `replace \`ifvar' = 0 if \`nopre' == 1`, then `mint`/`mintp1`
  recomputed on the post-drop sample.
- **This run (2026-09-11), `trig_mintsent_off`** (blind-spot config, `dta_s_mintsent`
  with sentinel-0 onsets, `safety(off)`): `A15-mint-s` — N 296 → 290, N_dropped 4 → 10,
  Nclus 30 → 29, df_r 29 → 28, ll −327.42 → −323.75, r2 0.5892 → 0.5786;
  `A15-mint-e` — every e(sample) digest shifts (`__N` 296 → 290, sums of i/t/ttre/w/y);
  `A15-mint-m` — `b` 3 cells worst **1.84e-01** ([1,1] att 0.90455 → 1.06361; [1,5]
  0.96680 → 1.18453; [1,6] 1.61506 → 1.78800), `V` 6 cells 1.59e-01, `wdidmainresults`
  33 cells 9.07e-01. The default-safety companion `trig_mintsent_def` is the F13 refusal
  (rc-by-design 0 → 459). Identical to 2026-08-23 (which recorded N 296 → 290).
- **Reproduction:** not re-run this pass. **No two-build do-file for this family** (see L.9.3
  for the full coverage statement). What stands in its place is a **re-tie of the decisive cells
  to a dump that exists today**: the 2026-09-13 dual-build run in `test/suite/out/run_old|run_new`
  is on the *same pinned pair* this section is about — `prov|old|301b706` and `prov|new|6a05127`,
  whose `wooldid.ado` is byte-identical to the pin `08ec0b4` (md5
  `382c1e0ee3e1be37fa2e9b677a65fdcc`). Read 2026-09-14; artifact `critleft/retie_0914.txt`. On
  `trig_mintsent_off`:
  ```
  N          old 296            new 290          N_dropped  old 4    new 10
  Nclus      old 30             new 29           df_r       old 29   new 28
  ll         old -327.41807     new -323.7506043 r2         old 0.589158811  new 0.5786458842
  e(b)[1,1]  old 0.9045468954   new 1.063611069
  e(b)[1,5]  old 0.9667957393   new 1.184531679
  e(b)[1,6]  old 1.61505938     new 1.787992707
  ```
  Every figure quoted above reproduces. One refinement the re-tie adds: the cluster count is
  `e(Nclus)` on **both** builds (30 → 29); new additionally posts `e(N_clust)` = 29, which old does
  not post at all — an addition alongside the move, not the move itself.
- **Classification: old-bug-fixed.**

### L.7 The esf calendar-support drop ran AFTER the reference checks (G10B) — ledger `G10B-esf-s/m/e/n`

- **Symptom / impact:** in the `esfixedbaseperiod` branch old ran the per-group
  reference checks BEFORE the calendar-support drop; when the support drop deleted a
  period containing a cohort's fixed reference row, that cohort survived without a
  reference — an opaque r(322) under default safety and silently wrong estimates at rc 0
  under `safety(off)`, the switch the error message steers users toward.
- **Mechanism (v1.3, re-located):** old esf branch lines 630-640 build
  `refperiodobserved` and drop refless cohorts; the timetrends reference block follows
  (643-653); the support drop (`controlcheck2`, `replace \`ifvar' = 0 if
  \`controlcheck2' == 0`) only comes at 661-676 — after both reference checks.
- **Fix + proof:** `7ba1e16` (VAL-G10B, 2026-08-19): minimal support-first reorder with a
  fixed-point proof (support is furnished only by never-treated rows, which the
  reference screens never drop), two review rounds, benign esf configs bit-identical;
  HEAD contract pinned by `test/test_esfsupport.do`; parent probe att −0.151 vs
  hand-deleted-cohort benchmark 2.005 (`wooldid_solved_issues.md`, external campaign
  entry, and the 2026-08-23 coverage-closure record). **At HEAD** new lines 800-830 run
  the `controlcheck`/`controlcheck2` support drop first (with the VAL-G10B comment
  block), then `refperiodobserved` at 835-845 and the timetrends block after it.
- **This run (2026-09-11), `trig_esfsup_ref`** (deterministic 14-unit fixture,
  controls unobserved at t = 3, cohorts onset 4 — reference t = 3 unsupported — and
  onset 5; `safety(off)`): `G10B-esf-s` — N 70 → 50, N_dropped 8 → 28, Nclus 14 → 10,
  df_r 13 → 9, ll 48.55 → 26.27, r2 0.9618 → 0.9298; `G10B-esf-e` — `__N` 70 → 50 and the
  digests; `G10B-esf-n` — `b`/`V`/`wdidmainresults` 1x6 → 1x5 (the `ES_att_pre2` cell,
  whose rel −2 row sat only in unsupported t = 3, is gone with the cohort); `G10B-esf-m`
  — `b` [1,1] **0.61487 → 1.41062** (the 2026-08-23 three-value triangle: old full
  fixture 0.6149 silently wrong; old on the hand-benchmark sample 1.4106; new 1.4106 to
  16 digits), `V` 22 cells, `wdidmainresults` 39 cells. `trig_esfsup_ref322` (default
  safety) is pinned 322 → 0. Identical to 2026-08-23.
- **Reproduction:** not re-run this pass; the committed `test/test_esfsupport.do` pins the HEAD
contract. **No two-build do-file for this family** (see L.9.3 for the full coverage statement).
  What stands in its place is a **re-tie of the decisive cells to a dump that exists today**: the
  2026-09-13 dual-build run in `test/suite/out/run_old|run_new` is on the *same pinned pair* this
  section is about — `prov|old|301b706` and `prov|new|6a05127`, whose `wooldid.ado` is
  byte-identical to the pin `08ec0b4` (md5 `382c1e0ee3e1be37fa2e9b677a65fdcc`). Read 2026-09-14;
  artifact `critleft/retie_0914.txt`. On `trig_esfsup_ref`:
  ```
  N          old 70            new 50           N_dropped  old 8   new 28
  Nclus      old 14            new 10           df_r       old 13  new 9
  ll         old 48.55284255   new 26.26879591  r2         old 0.9618496683  new 0.929840854
  e(b)[1,1]  old 0.614874351   new 1.410624802
  matnames b old 1x6|y1|main_att pretreatment_att ES_att_pre3 ES_att_pre2 ES_att_post0 ES_att_post1
             new 1x5|y1|main_att pretreatment_att ES_att_pre3 ES_att_post0 ES_att_post1
  ```
  Including the decisive one: `ES_att_pre2` — the rel −2 row that sat only in the unsupported
  t = 3 — is present on old and gone on new, and `b[1,1]` 0.61487 → 1.41062 is the
  three-value triangle's old-full-fixture and new values to the digits quoted. Same `Nclus` /
  `N_clust` refinement as L.6 (14 → 10 on `Nclus`; `N_clust` added on new only).
- **Classification: old-bug-fixed.**

---

### L.8 Section 2 re-checked at HEAD `08ec0b4` — "old-build issues NOT visible as old-vs-new mismatches"

Status of each of the eight 2026-08-23 items, in order. Three RETIRE (their defect was
fixed after 2026-08-23 and now fires as an ordinary old-vs-new difference), one closes
its open clause, one is now visible as an rc change but stays a defect, three are
unchanged.

1. **`ccc_absorb()` — still a both-builds defect, but now VISIBLE.** `cov_ccc_absorb`
   is pinned `198 → 459` (`RC_EXPLAINED[rc-by-design]` this run) instead of the 08-23
   both-r198 parity: since the BS-33 guard (`f8fe14e`, 2026-09-06) plus the
   cause-naming message of `5c78560` (2026-09-07; HEAD line 1614: "the cohort-by-period
   treatment terms were absorbed or dropped by the specification rather than estimated:
   check absorbed fixed effects (ccc_absorb() in particular ...)"), the new build stops
   with a named r(459) where old died with a bare "invalid syntax" r(198) in the margins
   phase. The underlying defect — `ccc_absorb` absorbs the treatment interactions
   themselves, so the fitted model carries no treatment coefficient at any horizon — is
   unchanged and still an open ticket (`configs.csv` note: "PARITY DEFECT TICKET ...
   unchanged and still open"). Kind: both-builds-defect.
2. **`o.`-marker naming fragility on `ext_contract10_rcs_factor_fe` — RETIRED.** SI-7
   was fixed 2026-08-27 (`configs_ext.csv` note: "was 0->198 (bare '5o. invalid name')
   ... New reads the marker, drops that horizon from sg1 with a caution, and completes").
   Both builds now rc 0 and the config fires the value/name rows `SI7-eb-n`, `SI7-eb-m`,
   `SI7-sg-n`, `SI7-sg-m`. **CORRECTED 2026-09-13:** an earlier draft of this item said
   old reported sg1's non-estimable pre2 horizon as a machine-zero 8.9e-16 estimate and
   new drops it. On the run's default path it does not. `contsg` S4 reproduced the config
   on both builds: the 8.9e-16 is **retained bit-for-bit on both** (`%21x`
   `+1.0000000000000X-032` on each side), and `e(b)` goes **11 names → 12**, an
   *insertion* of `pretreatment_SG1_att`, not a drop to 10; `wdidsgresults_1` goes 3x7 →
   4x7 with rownames `[att pre2 post0]` → `[att att pre2 post0]`. The drop the retired
   note describes happens only on the `fullprecisionslow` leg, because route P's cell
   collapse (`4b224f0`, default since `111ecbb`) removes the omission marker the SI-7
   read keys on. The parenthetical traced to `configs_ext.csv`'s then-current note and to
   `7f704d3`'s commit message, both of which describe pre-route-P behaviour; the ledger
   itself now agrees with contsg — at `4280f3b` `SI7-eb-n` reads *"e(b) gains the
   pretreatment_SG1_att entry (default cell-collapse path); 11 names -> 12 on 2026-09-11
   (the reason once said 11 -> 10: a drop that the collapse path superseded)"*
   (see `LEDGER_DELTA.md` §B). The colleagues' Stata-17 old build hit
   the same bare error, which is why this was a real hazard on old too. It belongs to
   the post-08-23 families (another agent's proof), not to this group — and that hand-off
   goes to **`contsg` S4's verdict of *unclear***, not to a clean fix.
3. **Unit-level `i` at realistic scale exits r(907) — the "pre-flight size hint remains
   ticketed" clause CLOSES.** SI-8 landed in `083ea89` (2026-08-25, author-approved):
   under default safety, more than 100 distinct `i` values draws a caution steering to
   `fe()` + coarser `i`, and a variable-demand estimate exceeding `c(maxvar)` halts EARLY
   with an informative r(907) naming the demand (HEAD lines 1263-1295); both skipped
   under `safety(off)`; cert `test/test_si89.do`. The limit itself is real and
   unchanged — `adm_unit_i` still runs the ~280-county subsample by design — what
   changed is that the failure is now early and named rather than cryptic and late.
4. **Stata-17 environment-dependent rc pins — unchanged.** Pins remain calibrated to
   this environment (Stata 19 MP); e.g. `ext_path10_ppml_separation` runs rc 0 on both
   builds here where the colleagues' Stata-17 old build recorded r(503) — Stata's
   conformability error, which is what its `msg_old` column in `configs_ext.csv` says.
5. **Corrected T1 reading (old-build missing pairwise-SG cells = the %21x encoding of
   Stata missing, not maxdouble artifacts) — unchanged, and now re-measured.** The
   *encoding* reading stands and is part of the report of record. Re-measured
   2026-09-14 on the 2026-09-13 dual-build dumps: across the four DDG configs old carries
   **60 missing `pwsgs` cells, every one of them exactly `+1.0000000000000X+3ff`**
   (= 2^1023 = Stata system missing), and **zero** maxdouble-like artifacts — 120 `pwsgs`
   cells scanned in `run_old/ext_cov10_rcs_{pweight_ddd,twoway_cluster}_c{1,2}.csv`.
   **The hand-off that followed it does NOT stand, and this item now defers.** It used
   to say the DDG family "is re-verified by another agent"; the other agent did not
   confirm it, and on the primary sources that section is right. `inference` §5 classifies DDG
   **unclear**, and that classification is adopted here. **The disagreement the earlier
   draft left open is settled 2026-09-14 in `inference` §5's favour**, from the two
   builds, the dumps and the ledger row (`CRITIQUE_CLOSEOUT.md` item 1; working dir
   `critleft/`):
   - the real → missing move is real, on **two** configs: `joint_pre_att_allsgssame`
     `+1.ffffb8906f474X-001` (= `.999997871047198`) old → `+1.0000000000000X+3ff` (= `.`)
     new, identically on `ext_cov10_rcs_pweight_ddd_c1` and
     `ext_cov10_rcs_twoway_cluster_c1`;
   - the "recovery" on the `pwsgs` surface of `..._pweight_ddd_c1` is a **relabelling,
     bit-for-bit**: old's `att_i:2vs1.group` row (`+1.ffffffffffff8X-001`;
     `+1.0000000000000X-020` = 2^-32; `+1.ffffffffffff8X+01f` = 2^32;
     `+1.55fa32d00a453X-2af` = 2.08e-207; and its two CI cells) is **byte-identical to
     new's `att:2vs1.group` row and to new's `att_itime:2vs1.group` row**, while new's own
     `att_i` inference cells go missing. Census of that matrix: 48 cells, **15
     missing → value, 5 value → missing**, 20 bit-identical;
   - the ledger has already moved the same way. At `26ce421` `DDG-jt-1` and `DDG-jt-3`
     both read "same"; at `4280f3b` = `HEAD`, `DDG-jt-1` reads *"joint_att_it_allsgssame is
     recovered (. -> value) AND joint_pre_att_allsgssame is withheld (0.999998 -> .) -
     this one row blesses both directions; the withheld move is recorded here not
     re-diagnosed (hygiene 2026-09-13)"* and `DDG-jt-3` *"here ONLY the withheld direction
     fires (joint_pre_att_allsgssame 0.999998 -> .)"*. That wording is suite hygiene, not an
     author ruling on DDG — none was found — which is exactly why the class is `unclear`
     rather than either settled verdict;
   - the suite's own degeneracy component corroborates "degenerate on both builds":
     `test/suite/degeneracy_expected.csv` line 173 classes `ext_cov10_rcs_pweight_ddd_c1`
     **RULED / `PWSGS_DEAD`, citation `DDG-pw-1`**, and the 2026-09-13 scan of the **new**
     dump records `PWSGS_DEAD(0v+2i/8) SGDEGEN(8:report)`
     (`test/suite/out/degeneracy/degeneracy_table.csv`) — two of eight pairwise rows are
     still inference-dead *after* the "recovery".

   **One thing is still wrong in the repo and is not this report's to fix.** `DDG-pw-1`'s
   reason is **byte-unchanged** from `26ce421` and still reads *"v1.3 returns MISSING
   pairwise-SG cells … where HEAD returns valid values"* — which the five value → missing
   cells and the `PWSGS_DEAD(0v+2i/8)` scan contradict on that same config. `4280f3b` gave
   the two-direction wording to the `-jt` rows only, and expanded `DDG-pw-2` / `-pw-3` from
   "same" to *"as DDG-pw-1 … missing-cells-recovered class"* without changing the framing.
   Recorded here, not edited: `test/` and `wooldidinstall/` were read-only for this agent.
6. **T-M2 (continuous-SG joint pre-trend p stored as literal missing on BOTH builds) —
   RETIRED; the adjacent allsgssame-ame observation is superseded on `cont_sg*` and
   otherwise UNCLEAR.** SI-10/T-M2 was fixed 2026-08-25 (author-ordered; the
   `r(p)[1,2]` read switched to `[1,1]` at the two sg-continuous sites) and now fires as
   `TM2-sg-p` on `trig_tm2_sgjoint` and `TM2-ext-p` on the 24-config `ext_c2_*` grid —
   an ordinary old-vs-new difference, handled by the TM2 agent. Precisely what the
   2026-09-11 run shows on `trig_tm2_sgjoint`: old stored the MACROS `joint_Choriz_s1..s3`
   = "." ; new REMOVES those macros (CS-60, `605d6d7`: per-subgroup joint p's became
   scalars under the same names) and posts the SCALARS `joint_Choriz_s1` = 0.8140,
   `joint_Choriz_s2` = 1.41e-13, `joint_Choriz_s3` = 1.44e-36 (`run_new/trig_tm2_sgjoint.csv`,
   decoded from `%21x`) — the real per-subgroup p-values the `[1,2]` read used to
   discard, delivered as scalars rather than as the macros the ledger row's reason text
   anticipated. The "adjacent open observation" (the continuous
   SG-equality `allsgssame` ame joint tests missing-valued in both builds, undiagnosed)
   is now, on the `cont_sg*` configs, skipped by design under AMESE R3 (`7054878`,
   cell containment; ledger `AMESE-R1-ja`: "SG-equality ame joint test skipped under
   cell containment (old stored missing here - the allsgssame-ame open observation)"),
   and BS-15 (`cfac3f5`, 2026-09-10) blanks both equality tests under a rank-deficient
   contrast covariance with a caution — a different mechanism. Whether the old
   missing-valued continuous `allsgssame` on a NON-containment, full-rank design was
   ever diagnosed is not recorded anywhere I could find (`allsgssame` appears in
   `wooldid_solved_issues.md` only in the 08-23 final-record line) — left as an open
   question rather than asserted closed.
7. **`poisexpresults` without `poisson` silently accepted on both builds — RETIRED.**
   SI-9 (`083ea89`, 2026-08-25): now an option-contract r(198) (HEAD line 5305: "Error:
   poisexpresults requires poisson to have any effect; specify poisson, or drop
   poisexpresults."); `ext_contract10_poisexp_without_poisson` is pinned `0 → 198` and
   judged `RC_EXPLAINED[rc-by-design]` this run. With it the whole REPORT_2026-08-18
   ticket-4 validation-gap family (elasticity type, `safety()`, `poisexpresults`) is
   closed.
8. **Both-error parity pins — count unchanged at 42 configs** (identical rcs on both
   builds, pinned by return code only; list in `parity_pins.txt`: the `ext_val_*`
   validation deck, `scr_strclus` r109, `scr_impsg` r198, `trig_pois_sep` r322,
   `patho_nocontrols`, `patho_nonever`, ...). Count, not membership: `fcol_gate` is a
   post-08-23 config. The designed-rc-change count rose from 32 to 66: of the 42
   configs at `198 → 0`, two are the 08-23 bucket-A crash cures (`cont_uz_ss` B5,
   `adv_gapa` B16) and the other 40 are the new-only option surface (old r198 unknown
   option → new rc 0: the `boot_*` and `band_*` decks, `perm_*` sub-options, `fcol_*`,
   `level_analytic`, `patho_boot/gt/perm/ptp`, `scr_antic`, `scr_fastnose`, `scr_ptp`,
   `gt_basic`, ...), the section-3.4 class of the 08-23 report. The README's
   accepted-blind-spot list is now 7 configs (`cov_icbad2` joined the six).

---

### L.9 Open questions (consolidated)

1. **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E1]** **A5-cfx-m attribution (L.2)** — the repair proposed here **landed essentially as
   proposed**: `4280f3b` added `SI20-L-trig_cfx-cfxame` (config `trig_cfx`, kind `mat`, pattern
   `cfxame_*`) and narrowed `A5-cfx-m`'s pattern `cfx*` → `cfx*_it_*`, so the row no longer
   absorbs the SI-20 cells by first-match. The question as written was:
   the ledger row explains, by first-match, nine SI-20
   value cells on `trig_cfx`'s `cfxame_*` matrices that its B7-only reason does not
   cover. Not a wrong verdict — the cells ARE a documented fixed class — but a
   reason/record mismatch. Suggested repair: a `cfxame*`-scoped SI-20 row ahead of
   `A5-cfx-m` (repo read-only for this agent; not done).
2. **Section 2 item 6, second half:** the continuous `allsgssame` ame equality tests
   being missing-valued on both builds outside the cell-containment regime is
   unrecorded as diagnosed or as closed; the `cont_sg*` configs no longer reach it.
3. **Reproduction coverage — stated explicitly, because it is not one per family
   (critique S-6).** The report's template promises every family a runnable two-build
   reproduction. **This section delivers two, for seven families.** That is disclosed here
   rather than papered over. `CRITIQUE.md` S-6 reports it as the only such gap in the six
   sections — it checked the other five and found every family there covered by a named
   reproduction block; that check is the critique's, not re-done here.

   | family | ledger rows | two-build do-file? | what stands in its place |
   |---|---|---|---|
   | L.1 A1/A2 (B1, B2) | `A1-B1a..e`, `A2-B2a..e` | **yes** — `repro_A1A2.do` | — |
   | L.5 A12 (tcen) | `A12-ttb-*`, `A12-ttp-s` | **yes** — `repro_A12TTB.do` | — |
   | L.2 A5-cfx (B7) | `A5-cfx-m`, `A5-cfx-n` | **no** | decisive cells re-tied 2026-09-14 (above); `3df74db`'s 3/3 verification record |
   | L.3 A6-ext (B9) | `A6-ext-m` | **no** | decisive cells re-tied (above); `2adf447`'s 3+3 record |
   | L.4 A7-drop / A7-kit (B10) | `A7-drop-m`, `A7-kit-s/m` | **no** | decisive cells re-tied (above), incl. the `%21x` bit-identity the item's inference rests on; `40e1913`'s 3+3 record |
   | L.6 A15-mint | `A15-mint-s/m/e` | **no** | decisive cells re-tied (above); `test/test_mint.do` |
   | L.7 G10B-esf | `G10B-esf-s/m/e/n` | **no** | decisive cells re-tied (above); `test/test_esfsupport.do` |

   **Why only two.** This agent's brief asked for two reproductions and named which two; the
   remaining five were left on the run's dumps. That was a defensible allocation when written
   and is a real gap in the assembled report, because the dumps a reader was pointed at
   (`test/suite/out/` — gitignored) are regenerated in place by every run.

   **What closes the hole without five more Stata families.** The 2026-09-13 dual-build run
   is on the **same pinned build pair** this section is about — `run_old` carries
   `prov|old|301b706` and `run_new` carries `prov|new|6a05127`, whose `wooldid.ado` is
   byte-identical to the pin `08ec0b4` (md5 `382c1e0ee3e1be37fa2e9b677a65fdcc`) — so the five
   families' decisive cells were re-read out of it on 2026-09-14 and printed in each family's
   Reproduction bullet above. **Every figure reproduced**, including the two the items' own
   reasoning leans on: `trig_allest_drop`'s `main_att` / `main_att_it` bit-identity at
   `+1.12b0e168eb563X+001`, and `trig_esfsup_ref`'s `ES_att_pre2` present on old and absent on
   new. Artifact: `critleft/retie_0914.txt` (the script that produced it, `retie_final.py`,
   sits beside it). This is a re-tie to a live dump, **not** a substitute for a two-build
   reproduction: it re-runs nothing, so it cannot catch anything the suite harness itself would
   get wrong. The honest statement for the assembled report's front matter is: *legacy ships 2
   reproductions for 7 families; the other 5 are evidenced by the dual-build dumps, re-tied
   2026-09-14 to the 2026-09-13 run on the same pinned builds, plus the committed cert suites.*

### L.10 Files

The runnable do-files and helpers below are at
`review_process/assessment_2026-09-13/repro/legacy/` (parameterised). The rest are artifacts of
the drafting session and are regenerable by running them; see `repro/README.md`.

- `SECTION.md` — this report. `NOTES.md` — running notes.
- `old/wooldid.ado`, `new/wooldid.ado` — the pinned builds (`301b706`, `08ec0b4`).
- `repro_A1A2.do`, `repro_A12TTB.do` — the reproductions (now `args build repo WD_`: one build
  per
  Stata process); `run_repro.sh` — sequential driver with cap check and watchdog;
  `repro_cross.py` — decodes the `%21x` lines, prints the triangles, ties every common
  cell to `test/suite/out/run_old|run_new`.
- `run_repro.out` (driver console incl. every `REPRO_*` line; its tail is the
  regenerated `repro_cross.py` summary — the first pass of that script had imported
  `test/compare.py`, whose own main ran instead and printed a stray
  `TOTAL_DIFF_GROUPS=0`; the parser is now inlined and the file repaired), `repro_cross.out`
  (the decisive lines quoted above), `repro_A1A2_{old,new}.log`,
  `repro_A12TTB_{old,new}.log` (full Stata logs), `dumps_{old,new}/{basic,trig_tt_big}.csv`
  (`%21x` e() dumps from the reproductions).
- `parity_pins.txt` — the 42 both-error parity configs of the 2026-09-11 run.
- `tie_{old,new}/{basic,trig_tt_big}/{a,b}/` — pass 2's offline re-tie inputs (`a` = the
  0911 run dump minus prov/varlist/esamp/rc rows, `b` = the reproduction dump minus rc),
  compared with `test/compare.py a b --rtol 0`.


---


# Part 2 - SI-20: the two model-construction fixes


> ## Headline of record and evidence provenance (added 2026-09-13; identical in all six sections)
>
> **Superseding note added 2026-09-14.** Two things below have since moved and are corrected
> here rather than rewritten in place, because the block is accurate about the run it describes.
> (i) `4280f3b` is no longer `HEAD`; the report's own header names the run of record by its
> `prov[new]` sha, and **that** is the build the report's claims attach to. (ii) The closing run
> of this assessment is the **2026-09-14** one over a **559-config** universe against the same
> v1.3 old side, reported in Part 7.5: **559 of 559 compared, zero DIFFERS, unexplained 0,
> stale_ledger 0**. The 2026-09-11 run described below is where this section's EVIDENCE was
> gathered and its magnitudes measured; the 2026-09-14 run is what verifies that nothing still
> differs. Neither supersedes the other - they answer different questions.
>
> **Run of record.** The 2026-09-11 dual-build run: **old = `301b706`** (published v1.3),
> **new = `14c30d3`**, whose `wooldid.ado` is byte-identical to this section's pin `08ec0b4`
> and to current `HEAD` apart from a nine-line comment block and a shortened `dis as text`
> permutation caution (10 insertions, 1 deletion; no estimation code).
>
> **The report judges those dumps against the ledger and comparer at `4280f3b` = `HEAD`.** That
> suite defines a **538-config universe**; the 2026-09-11 dumps cover **507** of them:
> **507 compared — EXPLAINED = 374, MATCH = 25, RC_EXPLAINED = 108 — unexplained 0,
> stale_ledger 0, unscoped 0, prov_fail 0, COVERAGE_OK.** The other 31 configs are new at
> `4280f3b`, have no 2026-09-11 dumps, and are outside every claim in this report.
>
> An earlier headline — **507/507 compared, EXPLAINED = 399, RC_EXPLAINED = 108**, at ledger
> `26ce421` — is quoted by three of the six drafts (`legacy` and `intended` give the full
> 399/108 line; `inference` gives the 507/507). It is the **same dumps and the same 507 configs**
> judged by the pre-`4280f3b` comparer. The two reconcile exactly: `399 = 374 + 25`, the 25 being
> the configs on which MATCH became reachable once the comparer began comparing across batch 3's
> macro→scalar kind change (`4280f3b`, verbatim: *"It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)"*; the same paragraph
> puts the cells affected at *"361 values on 103 configs"*). **Both numbers are measured and
> neither is withdrawn.** See `LEDGER_DELTA.md` §0.
>
> **Evidence provenance.** Every `test/suite/out/run_old|run_new/<config>.csv` cited below is a
> **regenerable** artifact, not a stored one: `test/suite/out/` is gitignored (`test/.gitignore`
> line 5) and the 2026-09-11 dumps were overwritten in place by the 2026-09-13 run. Regenerate
> them with `./run_full.sh oldnew` from `test/suite`, which stages old from
> `git show main:wooldidinstall/wooldid.ado` and new from `git show HEAD:…`. `compare_verbose_0911.txt`
> — cited by four sections, twice as `scratchpad/compare_verbose_0911.txt`, which is **not**
> repo-relative and never was — is likewise a recipe: rebuild it with `compare_oldnew.py … --verbose`,
> passing `--ledger`/`--configs` from whichever commit's judgement you want. The command, the build
> table and the one harness caveat (`4280f3b` changed `mkrun.py`'s `dumpres`, so a regenerated dump
> may differ in its `files|` line or gain a `matnames|…|UNREADABLE` line — **value cells are
> unaffected**) are in **`review_process/assessment_2026-09-13/repro/README.md`**.
>
> **Reproductions.** The do-files this section cites "in this directory" are now parameterised —
> repo root and scratch workdir as positional arguments — and live at
> **`review_process/assessment_2026-09-13/repro/<section>/`**. Their exact invocations are in each
> file's header and in that README. Bulk artifacts (dumps, logs, copied `.ado` builds) are not
> copied: they are regenerable, and the README says how.
>
> **Per-edit log:** `EDITS.md`. **Ledger delta and closed open questions:** `LEDGER_DELTA.md`.


Group: **SI-20** (ledger families `SI20-L-*`, `SI20-I-*`, `SI20-P-*`), the largest explained class
of the run. Builds: **old** = published wooldid v1.3 (`git main`, dump prov `old|301b706`) vs
**new** = git HEAD `08ec0b4`. The 2026-09-11 dumps in `test/suite/out/run_new` carry prov
`new|14c30d3`; `git diff 14c30d3 08ec0b4 -- wooldidinstall/wooldid.ado` is a comment block plus a
shortening of the permute off-centre caution's display text (`08ec0b4`, audit #57; no estimation
code), so the dumps are the HEAD ado numerically. Sources of record: `wooldid_solved_issues.md`
"SI-20" (l.517-601), `review_process/si20_ship_2026-09-05.md`, `wooldid_ado_review.md` "Findings
(Claude, 2026-09-05)", the ledger comment block (`test/suite/ledger.csv` l.128) and its 206 `SI20-*`
rows, the verbose re-judge `compare_verbose_0911.txt`, fix commit `1df8613` (2026-09-05).

Old-code line numbers below refer to v1.3 (`301b706`) and were located in that build directly.

**This run, in one line:** 206 ledger rows (149 `SI20-L`, 38 `SI20-I`, 19 `SI20-P`) **as the
ledger stood at `26ce421`** — at `4280f3b`, the ledger of record (see the headline block), it is
**207 / 150 `SI20-L`**: `4280f3b` split `SI20-L-trig_cfx-cfxame` out of `A5-cfx-m`, and that row
fired, since the HEAD re-judge of these same dumps reports `stale_ledger=0` (`LEDGER_DELTA.md` §1,
§D). Of the 206 measured here, **all 206
fired**, 532 firing lines on **53 configs** (43 with an L row, 9 with an I row, 13 with a P row); 0
stale, 0 unscoped. The solved-issues entry and the ship record speak of "208 ledger rows": the
ledger held 208 at the ship commits (`42ff816`, `1661415`); `8379d29` (the R4-refit withdrawal,
author ruling 2026-09-08) retired the two `si20_lat_cfx_pw` rows `-cfxatt_` and `-properties` that
existed only for the refit's cfxatt drift, leaving the 206 that stand at HEAD — every one of which
fired.

Reproduction: `repro_si20.do` (now at `review_process/assessment_2026-09-13/repro/si20/`, with
`instrument.py`, `check_inert.py`, `summarize.py`, `named_cells.py`; one build per Stata process, run from `run_old/`
and `run_new/`; sentinels `REPRO_SI20F1_<cfg>_DONE`, `REPRO_SI20F2_<cfg>_DONE`, `REPRO_SI20_DONE`).
The two builds are **instrumented copies** (`old_i/`, `new_i/`, made by `instrument.py`; `diff`
against the pristine extractions shows 52 inserted lines and 0 deletions each): immediately before
the two pooled joint-test `margins ... over(_t) ... contrast(overjoint)` calls — old lines 1773
(binary att) and 2148 (continuous ame); new 2170 and 2662 — the same call is re-run without the
contrast under `cap qui`, never posted, and its `over(_t)` level margin for the reference cell (the
column whose name starts `0b`) is printed as `REFM|<F2|F1>|<eq>:<col>|<value>`. The probe is inert:
`check_inert.py` shows every e(b)/e(V)/e(wdidmainresults) cell and every e() scalar of the
instrumented runs **bit-identical (`%21x` string equality, rtol 0)** to the suite's 2026-09-11
dumps on both builds (old: 254/254, 254/254, 119/119 cells; new: 265/265, 265/265, 130/130).
Wall time 39-52 s per build.

---

### S.1 `lattice()` + `contreatcontrols(basic)`: v1.3 absorbed no baseline dose control on the reghdfe path, so every ame-family row carried the baseline dose slope (SI-20 fix 1)

- **Symptom / impact:** with a continuous treatment, `contreatcontrols(basic)` and any `lattice()`,
  v1.3's default (reghdfe) path inserted **no** baseline dose control at all. The base level of the
  dose therefore stayed estimated inside the cell block instead of being partialled out, the
  event-study **reference cell carried a non-zero average marginal effect**, and every ame-family
  row the user reads — `main:ame`, `pretreatment:ame`, each `ES_ame:` horizon, the per-subgroup
  ame rows and their pairwise contrasts, the `cfxame_*` bins — was off by the baseline dose slope,
  O(1). On the suite's fixture the pre-treatment ame rows printed 0.81-0.90 with intervals
  excluding zero where the corrected values are -0.10 to -0.25; on `dta_cont` (DGP base slope 2)
  they printed ~2.0 where the truth is 0 (ship record §3 A: on a null outcome v1.3 reports
  `main_ame` 1.99 for a truth of 0). Anyone who published a continuous-treatment ame table with a
  lattice published the dose slope added to every effect. The att rows were never affected.
- **Mechanism (v1.3):** the guard at old line 975,
  `if (strpos("\`contreatcontrols'","basic") != 0 & "\`latticepart2hash'" == "" ) | "\`unconditionalse'" == "unconditionalse" | "\`customvce'" != "" | "\`poisson'" == "poisson" {`,
  is the only place `basic` initialises the absorbed dose control (`local contreatcontrolinserterFE
  i.\`null'##\`ringerpoly'` at 978, plus `i.\`latticevar'##\`ringerpoly'` at 979 when a lattice is
  present). With a lattice (`latticepart2hash` = `##i.<lattice>`) the first disjunct is false, and
  on the default path the other three are false too, so the block is skipped and nothing else
  inserts a dose control; the `i.\`_i'#i.\`_t'#i.\`targetflag'##i.lattice##c.contreat` cell block
  (reghdfe at old 1147) then keeps its base-level dose columns. The three escape disjuncts are
  exactly the cures the ticket measured (`unconditionalse`, `customvce`, `poisson` restore a zero
  anchor), and the author's own comment at old 970-972 says that specifying `null##ringer` and
  `lattice##ringer` explicitly "ensures that the base level of contreat gets partialled out. if not
  specified, it does not always get partialled out". The pooled ame joint test (old 2148,
  `margins, dydx(contreat) over(_t) subpop(... (_t == 0 | _t >= horizonstarter)) contrast(overjoint)`)
  contrasts every pre horizon against that non-zero reference margin.
- **Fix + proof:** `1df8613` (2026-09-05, branch `si20-fix`), the `& "\`latticepart2hash'" == ""`
  disjunct deleted (the guard at new line 1210, block 1210-1216): `basic` always absorbs `i.null##c.(ringer)` and, with a
  lattice, its `i.lattice##c.(ringer)` twin; help `lattice()` (sthlp 909-912) and `basic` (964-965)
  document the change and date it. Ship record §3 A: 183 lattice configs on both builds — on all 59
  plain-basic + lattice configs that report ame rows the pristine reference margin is non-zero
  (0.27..2.00; flagship 1.019826251639038 = the ticket's value) and the fixed build's is exactly 0
  with the base-level dose columns omitted (`o.dose_common 2o.lat#co.dose_common`); every
  ame-family row moves by exactly minus the pristine anchor (<= 2e-14) where the lattice
  composition of the reference cell equals each ES cell's, row-specific (2e-9..9e-3) under
  person-varying weights or a time-varying dose; att rows unchanged (<= 1e-10); byTime-only,
  untreatedZero, unconditionalse, customvce and poisson configs bit-identical (61 configs, 14,541
  cells); rc identical on all 183. Solver dependence removed (pristine: `regtol` 9 -> 12 moves the
  anchor 0 -> 0.1313; fixed: tolerance-invariant). Battery signature exactly as pre-registered:
  only `cont_lat`, its six ame estimates by -2.0043 (§2). Cert `test/test_refmargin.do` (49c9fc3):
  ALL_PASS on the fix, 52 failures on the pristine build, calibration cases passing on both.
- **This run (ledger `SI20-L-*`, 149 rows, 43 configs, 428 firing lines):** scoped by config name
  and unbounded by design (the ledger block: "so a future regression there still fails the run" is
  reserved for the parity pins, S.3). The signature on every plain-basic + lattice reghdfe-path
  config is identical in kind: **ame-family cells move by the baseline dose slope; att point
  estimates unchanged (<= 1e-9); df_r -1** (the absorbed slope is counted) so att-side se/t/p/ci
  move by sqrt(df_old/df_new) (3.9e-4 rel on se; p and near-zero ci cells tail-amplify); `F ll_0
  r2_within* r2adj rmse` follow reghdfe's absorbed-model accounting; `ll rss r2 N` unchanged.
  Measured from the dumps (`summarize.py`; `main:ame` old -> new; shift = old-new on every ame
  entry of e(b)):
  - `si20_lat_pw` (pw): 1.24375 -> 0.24316, shift 1.000594 on all six ame rows, df_r 1268 -> 1267,
    `joint_Choriz` 0.3113 -> 0.3116, V worst 9.06e-01, wdidmainresults 66 cells worst 1.40e+00.
  - `si20_lat_aw_bt` (aw, baseTreated 2): 1.26644 -> 0.24316, shift 1.023278, df_r 1414 -> 1413.
  - `si20_lat_sg` (subgroups): 1.26735 -> 0.27118, shift 0.996171; `wdidsgresults_1..3` 42 cells
    each (worst 1.83e+00), `pwsgs` 60 cells (1.66e-02), `joint_*_allsgssame` at the third digit.
  - `si20_lat_esf` (esfixedbaseperiod, byTime lattice): 1.26735 -> 0.36292, shift 0.904425.
  - `si20_lat_poly` (contreatpoly(2)): shift 0.996171, df_r 1322 -> 1320 (two absorbed slopes).
  - `si20_lat_ccnoab`: shift 0.996171 with **no** df change (the ringer terms are explicit
    regressors under ccnoabsorb) — only the ame cells move, as the row says.
  - the 24-config class `ext_lat_*` (`SI20-L-ext_lat_glob-*`, 8 rows, 240 firing lines): the n2
    lattices 1.26735 -> 0.27118 (shift 0.996171, df_r 1450 -> 1449), the n3 lattices 1.24399 ->
    0.25570 (shift 0.988298, df_r 1322 -> 1321) — one number per lattice size, every lattice type.
  - `ext_cont_ext_custom_lattice`, `ext_cont_ext_cfx`, `ext_cont_ext_cfx_subgroups`,
    `ext_cont_ext_weighted_lattice`, `ext_cov10_cfx_inference`: shifts 0.988-1.000, plus the
    `cfxame_*` bin matrices (18 cells each, worst 1.95e+00; `cfxame_cfxpre` 0.755 -> -0.277).
    `ext_char_cfx_nonatt` (att_i): e(b) has no `main_ame` row; its ame entries shift 1.889935,
    df_r 272 -> 271.
  - the `dta_cont` class under `cluster(i)` (`cont_lat`, `plot_cont`, `trig_cfx`, `si20_lat_cfx_pw`;
    e(properties) `b` on **both** builds — v1.3 posted b-only silently, HEAD posts b-only with the
    all-zero-covariance caution under the 2026-09-08 ruling `b801282`, which removed the
    explicit-dummies refit before these dumps were made): `main:ame` 3.00432 -> 1.00001, shifts
    2.0040..2.0051 = the DGP base slope 2; df_r 9 -> 9 (cluster-count df, so no df factor here);
    att entries agree to 5.2e-12 (`cont_lat`, `si20_lat_cfx_pw`) / 5.45e-9 (`plot_cont`,
    `trig_cfx`). The ledger reasons on these rows (written at ship, `42ff816`) still attribute
    that residual to "the retired AMESE-R4-b regress-vs-reghdfe drift" and, on
    `si20_lat_cfx_pw-b`, say "properties [b] -> [b V]; att inference cells recovered"; neither
    holds on HEAD — both builds fit through reghdfe and post b only, so the 5e-9 is reghdfe
    convergence drift under a different absorbed set (the class `SI20-P-si20_lat_byt-b` bounds at
    7.07e-8). Stale reason text, not a stale row — open question (4). **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E5]** —
    `4280f3b` rewrote all three: `SI20-L-si20_lat_cfx_pw-b`, `SI20-L-cont_lat-b` and
    `SI20-L-trig_cfx-b` now state `e(properties)` is `[b]` on both builds and that no inference
    cell is recovered.
  - `si20_ic_cont` (both fixes at once): ame shift 0.944453..0.944475 and the att rows move too
    (worst 8.92e-02, `ES_ame_post0` 0.626 -> -0.319) — the interactivecontrols half is S.2.
  - `ext_san_input_missing` (pathological: 9 planted missings, 36 person clusters vs ~100
    parameters): ame shift 0.918..0.930 **and** properties `b V` -> `b`: the fixed build's extra
    absorbed slopes push this knife-edge fixture into the all-zero-VCE state, so no e(V) is posted
    and every se/t/p/ci is missing with the loud caution; att points agree to 3.2e-9. The b-only
    posting with the caution is HEAD's ruled contract for an entirely non-identified covariance
    (`b801282`, 2026-09-08: the estimates stand, every se is missing, no refit); the ledger
    reason's "the R4 refit cannot fire" predates that ruling. Observable failure, not wrong
    numbers — ship record §4.2.
  - `si20_lat_byt` (`basic byTreated` + lattice): the anchor was **already zero** (byTreated absorbs
    the base slope), the ame shift is 0 and the point estimates agree to 7.6e-13; what fires is
    df_r 1263 -> 1258 (the added ringer terms lie in the span of the absorbed byTreated x lattice
    slopes but are counted) with se 2e-3 rel, p/t tail cells up to 0.15. S.4 / open question.
- **Reproduction (`repro_si20.do`, sections 1a-1c; decisive lines pasted from `run_*/repro_si20.log`):**
  `si20_lat_pw` — `wooldid y_cont_common cohort t ttre [pw=probability_weight], att contreat(dose_common)
  contreatcontrols(basic) lattice(all 3) esprelength(2) espostlength(2) jointtests safety(off)` on
  `test/suite/data/continuous_control_dense.dta`:
  ```
  old:  REFM|F1|dose_common:0bn.__000003|1.00059351834833      REFMV|F1|...|.00452054133856487
        B|si20_lat_pw|main_ame|1.24375359087951    B|si20_lat_pw|pretreatment_ame|.883390245761555
        B|si20_lat_pw|ES_ame_pre2|.896798112890836   SCAL|si20_lat_pw|df_r|1268
        W|si20_lat_pw|pretreatment:ame|p|3.01434848490726e-88
  new:  REFM|F1|dose_common:0bn.__000003|0                     REFMV|F1|...|0
        B|si20_lat_pw|main_ame|.243160072531196    B|si20_lat_pw|pretreatment_ame|-.117203272586756
        B|si20_lat_pw|ES_ame_pre2|-.103795405457479  SCAL|si20_lat_pw|df_r|1267
        W|si20_lat_pw|pretreatment:ame|p|.136891542918375
  both: B|si20_lat_pw|main_att|1.3808962103 (old) / 1.3808962102 (new); REFM|F2 (binary probe) = 0
  ```
  Arithmetic on the dumps: `old - new` on all six ame entries = 1.000593518348 (each within
  2.6e-14 of the old anchor 1.00059351834833); att entries agree to <= 7.1e-10; `main:att` se
  new/old = 1.000394555 = sqrt(1268/1267) to nine digits. The binary probe (`F2`) at the same run
  reads `0bno.` = 0 on **both** builds — the att side never had an anchor (calibration).
  `si20_lat_aw_bt` (aweights): old anchor 1.02327770247951 (V 1.8026e-3), new 0; shifts 1.0232777
  on every ame row, within 2.7e-9 of the anchor (the person-varying-weights class the ship record
  names); df_r 1414 -> 1413, se ratio 1.000353795 = sqrt(1414/1413).
  Section 1c is the ticket's own fixture of record (unweighted `lattice(baseTreated 2)`, `robust`):
  ```
  old:  REFM|F1|dose_common:0bn.__000004|1.01982625163904      REFMV|F1|...|.00178439875953114
        B|ticket_bt2|main_ame|1.26734836006062    B|ticket_bt2|ES_ame_pre2|.806931487285221
  new:  REFM|F1|dose_common:0bn.__000004|0                     REFMV|F1|...|0
        B|ticket_bt2|main_ame|.24752210842158    B|ticket_bt2|ES_ame_pre2|-.212894764353825
  ```
  — the ticket's 1.0198262516390375 / V 1.784e-3 to every printed digit; `old - new` on all six
  ame entries = 1.019826251639036..046 (within 7e-15 of the anchor), att entries agree to <= 5e-14,
  df_r 1414 -> 1413, se ratio 1.000353794 = sqrt(1414/1413); `joint_Choriz` 0.00360 -> 0.00361 (the
  uniform shift cancels inside `contrast(overjoint)`, as the Findings block says).
- **Classification: old-bug-fixed.** A specification defect in the model construction (the author's
  own comment describes the intended partialling-out; the guard prevented it whenever a lattice was
  present), author-ruled and shipped as such. Not a convention change: the att side, the poisson
  path and every non-lattice `basic` run are unchanged.
- **Open questions:** (1) the df_a over-count (ship record §4.1): reghdfe counts the added ringer
  slopes, so **every** row's se/t/p/ci — att rows included — scales by sqrt(df_old/df_new) (+1 df
  on the plain path, +4..6 under `basic byX` + lattice where the terms are redundant with the
  absorbed byX x lattice slopes, V x1.0028..1.0050). No ruling on it was found in
  `qr_ship_2026-09-06.md`, `held_for_ruling_2026-09-09.md` or the rulings memory; the ledger rows
  document it as-built. (2) `ext_san_input_missing`'s loss of e(V): the ship record's §4.2 candidate remedy
  (widen the R4 refit trigger) was withdrawn with the refit itself by the 2026-09-08 ruling
  `b801282`; the b-only-with-caution posting seen here is that ruling's contract, so nothing is
  left to rule except whether a pathological fixture crossing the all-zero-VCE line needs more
  than the caution (the effort-cap doctrine says no). (3) The primary ppmlhdfe path with the absorbed
  ringer terms and `eyex`/`dyex` under a lattice are listed untested in the ship record §6;
  `si20_lat_pois` (S.3) shows the poisson path does not move, which is the guard's fourth disjunct
  being true there, not a test of the new absorbed terms under ppmlhdfe. (4) Stale reason
  text in the ledger: `SI20-L-{cont_lat,plot_cont,trig_cfx,si20_lat_cfx_pw}-b` still describe
  "this R4-refit config" / "properties [b] -> [b V]" and the `SI20-L-ext_san_input_missing-*`
  rows "the R4 refit cannot fire" — mechanisms removed by `b801282` on 2026-09-08, three days
  after the rows were written. Every one of these rows fires correctly on the fix-1 ame shift
  (S.1's signature), so this is a wording refresh for the ledger maintainer, not a comparer
  defect; `8379d29` already retired the two `si20_lat_cfx_pw` rows whose *firing* the refit had
  caused.

### S.2 `interactivecontrols()`: the pooled reference cell carried its own slope column, exactly collinear with `i.t#c.x` plus the treated-cell slopes, so the engine's collinearity drop decided the reported horizons (SI-20 fix 2)

- **Symptom / impact:** under `interactivecontrols()` v1.3 built the cell interaction with
  `i.targetflag`, which includes its base level, so the reference cell (`_i = _t = targetflag = 0`,
  the pooled untreated observations) got a slope column `0._i#0._t#0.targetflag#c.x_dm`. That column
  is exactly collinear with `i.t#c.x` plus the treated-cell slopes (modulo the absorbed cohort
  means), so reghdfe/reg/ppmlhdfe dropped **one** of the set — which one depended on regressor
  order, solver and engine — and the surviving reference slope `s000` leaked into every reported
  horizon: the reference-cell margin was `-s000 x E_w[x_dm | treated reference rows]`, non-zero
  wherever the control varies within cohort over time, under weights, or on an unbalanced panel.
  Consequences a user saw: every displayed horizon, `main:att` and `pretreatment:att` shifted by a
  horizon-specific amount (on the flagship -0.08..+0.14 against se 0.08-0.12), standard errors
  moved 18-40%, and the pre-trend joint test contrasted the pre horizons against the non-zero
  reference margin rather than zero (flagship p 0.0510 where the pre = 0 Wald on the same output is
  0.70). The Findings block's trending-covariate panel (true pre-effects 0, true ATT 1.0-1.5): v1.3
  printed pre rows -2.96/-2.10/-0.81, att 3.92, joint p 1.9e-38 — a wholly spurious pre-trend and a
  2.6x ATT — where the corrected construction gives ~0, 1.52, p 0.91. Every downstream consumer of
  `e(b_es_*)`/`e(V_es_*)` (HonestDiD, `pretrends`, `pretrendpower()`, the adversarial line,
  gteffects, the plots) inherited it. The reported numbers were not estimable functions of the model.
- **Mechanism (v1.3):** old line 933,
  `local interinject i.\`_i'#i.\`_t'#i.\`targetflag'#c.(\`interlist')  i.\`t'#c.(\`interactivecontrols')`,
  with `interlist` the within-cohort de-meaned controls (old 927, `gstats transform (demean)`); the
  comment at old 932 hopes "the null untreated dummy X interactive controls object included should
  be harmless". The pooled binary joint test (old 1773, `margins, dydx(targetflag) over(_t)
  subpop(... (_t == 0 | _t >= horizonstarter)) contrast(overjoint)`) includes `_t == 0`, so with a
  non-zero reference margin it tests "pre horizons = the reference-cell effect".
- **Fix + proof:** `1df8613` (new line 1152): `1.\`targetflag'` instead of `i.\`targetflag'` — the
  reference cell gets no slope column, its de-meaned-control variation is partialled out through
  `i.t#c.x`, and the reference margin is a structural zero (`0bno.`); help `interactivecontrols()`
  (sthlp 1051-1055) documents and dates it. Ship record §3 B: 116 config-build rows, 81 configs —
  fixed anchor exactly 0 on 92/92 running rows, pristine non-zero on 84/92 (the 8 zeros: no
  surviving reference-cell slope column); model fit identical across builds on every engine
  (rss/F/chi2 <= 1e-12 rel); on the fixed build regressor-order invariance <= 1.4e-13 (pristine
  moves up to 0.10, joint p .051 -> .639), reghdfe == reg <= 2.5e-13 (pristine 1.0-2.0), posted
  joint pre-trend p == the pre = 0 Wald from `e(b_es_*)`/`e(V_es_*)` <= 2e-13 across 23 configs
  incl. ppmlhdfe and reg (pristine fails wherever the anchor != 0; the ticket's .615 vs .146
  reproduced); replications under bootstrap/permute carry the fixed construction (per-rep anchors
  21/21 zero). rc identical on 115/116, the one mismatch pristine-only (r(459) VCV blackout under
  permtstat that the fixed build completes). `test_refmargin.do` detectors (a) joint-Wald identity
  and (c) the unconditionalse twin fail on pristine, pass on the fix.
- **This run (ledger `SI20-I-*`, 38 rows, 9 configs, 41 firing lines; plus `si20_ic_cont` under
  L rows):** the rows fire on b, V, wdidmainresults and, where present, the joint tests and
  subgroup objects. Measured from the dumps (`main:att` old -> new; worst e(b) cell):
  - `si20_ic_coh` (flagship: cohort-level i via prep `covars_coh`, `fe(i)`, `cluster(i)`,
    time-varying `x1`): 1.45625 -> 1.30051; `ES_att_pre2` -0.0473 -> +0.0352 (rel 1.74e+00);
    `pretreatment:att` -0.0567 (p 0.43) -> -0.0031 (p 0.97); **`joint_horiz` 0.05098 -> 0.73068**;
    V worst 4.22e-01; wdidmainresults 42 cells worst 1.78e+00. df_r 149 both.
  - `si20_ic_aw` (aw): 1.53256 -> 1.31547; `joint_horiz` 0.00200 -> 0.17038.
  - `si20_ic_pw` (pw, x1 x2, att att_it): 1.50836 -> 1.30216; `pretreatment_att` -0.036 -> +0.043.
  - `si20_ic_sg` (subgroups): 1.46821 -> 1.31169; `joint_horiz` 0.0489 -> 0.4608; per-subgroup
    joint p macros `joint_horiz_s1/s2` 0.031/0.075 on old are REMOVED on new (posted as scalars
    instead: 0.372/0.548, the row says); `joint_pre_att_allsgssame` 0.347 -> 0.020; `pwsgs` 12
    cells, `wdidsgresults_1/2` 36 cells each (worst 1.62e+00).
  - `si20_ic_esf` (esfixedbaseperiod, adversarial): 1.40124 -> 1.26548; `joint_horiz` 0.581 ->
    0.924; `adversarial_pval_att` 9.7e-11 -> 1.1e-11 and `adversarial_line_att` [t _cons] ->
    [slope intercept_rt0] — these two ride on the NFADV-01/02/05 recalibration and the BS-12 ruling
    (template §3.5), an **intended change** inside this family; row `advln` is pure NFADV-05.
  - `si20_ic_admin` (realistic admin panel, pw, x1 x2, df_r 28642): 0.67129 -> 0.67080;
    `ES_att_pre2` 0.0632 -> 0.0753 (rel 1.60e-01) — small here because the control's within-cohort
    trend is weak, exactly as the mechanism predicts.
  - `si20_ic_pois` (ppmlhdfe, cohort-level i): 0.20568 -> 0.19903 (4.0e-02); V [1,1] 1.30e-3 ->
    2.42e-4 (8.1e-01); wdidmainresults 24 cells worst 1.00e+00 (t 5.70 -> 12.79).
  - `ext_bin_ext_interactive` (binary_balanced): 1.23766 -> 1.16794, `ES_att_pre1` 0.0069 ->
    0.0179 (6.14e-01); these rows absorb the retired bounded XIC-b/V/w B21 IC-copy rows (1.33e-4),
    now dominated by the fix. `ext_cont_ext_interactive` (contreat basic, no lattice): att rows
    move (8.59e-02) while the ame rows are untouched (shift <= 1e-6) — fix 2 without fix 1.
  - Independent corroboration inside the dumps (detector (c) of `test_refmargin.do`): on
    `dta_s_covars` with `interactivecontrols(x1)` the plain-reg twin `si20_ic_uncond`
    (`unconditionalse`; reg dropped the reference column itself) posts `main_att`
    1.3116904373277194 on **both** builds; the new build's default reghdfe path on the same model
    (`si20_ic_sg`) posts 1.3116904373277247 — agreement to 4e-15 — while v1.3's default path posted
    1.4682140200278555, 0.157 away from its own reg twin.
- **Reproduction (`repro_si20.do`, section 2; `dta_s_covars.dta` + `gen coh =
  cond(ttre==5,1,cond(ttre==9,2,3))`; `wooldid y coh t ttre, att cluster(i) fe(i) esprelength(2)
  espostlength(3) interactivecontrols(x1) jointtests`):**
  ```
  old:  REFM|F2|0b.__000008:0bno.__000004|0            (the vacuous 0.targetflag equation block)
        REFM|F2|1.__000008:0bn.__000004|.206686266294091   REFMV|F2|1.__000008:0bn.__000004|.00607539362946496
        SCAL|si20_ic_coh|joint_horiz|.0509757696006915
        B|si20_ic_coh|main_att|1.45624524730235   B|si20_ic_coh|pretreatment_att|-.0566563939724496
        B|si20_ic_coh|ES_att_pre2|-.0473364140370975
  new:  REFM|F2|0b.__000007:0bno.__000006|0
        REFM|F2|1.__000007:0bno.__000006|0                 REFMV|F2|1.__000007:0bno.__000006|0
        SCAL|si20_ic_coh|joint_horiz|.730675099835656
        B|si20_ic_coh|main_att|1.30051262684996   B|si20_ic_coh|pretreatment_att|-.00309508161032479
        B|si20_ic_coh|ES_att_pre2|.0351589506667569
  ```
  The old reference-cell margin 0.2067 with variance 0.006075 (se 0.0779) is the Findings block's
  "0.2067 (se 0.078, t 2.65)"; on the new build the column comes back omitted (`0bno.`) at exactly
  0, and `joint_horiz` 0.0510 -> 0.7307 is the number the ledger row and the config note carry. The
  `dydx(i.targetflag) over(_t)` margins return two equation blocks; the `0.targetflag` block is
  all-omitted on both builds and is printed only for completeness.
- **Classification: old-bug-fixed.** The reported quantities were not identified (solver- and
  order-dependent); the fix makes them estimable functions of the model. Author-ruled as a design
  ruling against the "harmless" base-level term (Findings §Solutions 2), documented in the help. The
  `si20_ic_esf` `advl`/`advln`/`advp` cells additionally carry the NFADV/BS-12 intended change.
- **Open questions:** none that the records contradict. The mechanism predicts (and the ship record
  measured) that where the control has no within-cohort variation the two builds agree — `cov_ic`
  and `si20_ic_uncond` (S.3) are that case and carry no SI20 row by design.

### S.3 The parity pins: `si20_lat_byt`, `si20_lat_pois`, `si20_ic_uncond`, `cov_ic`, and the B21 fit-statistic pins on the `si20_ic_*` configs

The ledger block places four configs where the fixes must **not** move the numbers, and gives them
bounded or no rows so that a future regression there still fails the run:

- `cov_ic` (`interactivecontrols(icvar)`, icvar has no within-cohort variation) and `si20_ic_uncond`
  (the plain-reg twin, `unconditionalse`; reg's own drop discarded the reference column): **no SI20
  row at all**, and in this run they fired only `CORE10-cmdline` (**retired at `4280f3b`**;
  `LEDGER_DELTA.md` §A) and `PS62-esblock` — e(b) is
  bit-identical old vs new (`main_att` rel 0.00e+00 on both). The design held.
- `si20_lat_pois` (`poisson` + lattice + basic): fix 1 does not reach the ppmlhdfe path (the guard's
  fourth disjunct was already true), so the rows are the v1.3-vs-HEAD IRLS class only:
  `SI20-P-si20_lat_pois-b` (rtol 1e-6) measured 2.47e-08, `-V` (1e-4) 4.97e-06, `-w` (1e-3)
  1.81e-05 (t/p amplify), `-rss` (1e-8) 1.9e-10. Bit-for-bit parity of the model on this path was
  established pristine-vs-fixed (156 cells) at ship time. (v1.3 needs 856 s for this config.)
- `si20_lat_byt` (`basic byTreated` + lattice): the anchor was already zero, so the point estimates
  are pinned — `SI20-P-si20_lat_byt-b` (1e-6) measured 7.07e-08 (the reghdfe tol-1e-9 class with
  more absorbed slopes), `-ll` (1e-8) 3.7e-10, `-rss` 9.5e-10 — while V is bounded at 1e-2 and
  measured 3.96e-03, the df factor 1263/1258. The unbounded L rows on this config are the df_r,
  F, r2*, rmse, wdidmainresults (63 cells, worst 1.52e-01 on p/t tails) and the two joint p's.
  **Classification: unclear (a category gap, not a missing proof)** — the b/ll/rss pins are
  parity assertions (nothing moved; no defect on either build); the V/df rows are the unruled
  df_a over-count of S.1's open question (1), a mechanical consequence of counting redundant
  absorbed slopes, not a change in what is estimated.
- `SI20-P-si20_ic_*-fit` (8 configs, rtol 1e-5): the B21 IC-copy class (v1.3 stored float copies of
  the de-meaned control, moving F ll rss r2* rmse, chi2/r2pseud under poisson) on the new
  interactivecontrols configs; measured ceiling 1.4e-7 (`si20_ic_admin` F). Fix 2 leaves the fit
  statistics identical (the model space is the same; only the parameterisation changed), so these
  rows are old-bug-fixed of a different, already-documented kind (B21), bounded so a fit-statistic
  regression cannot hide behind SI-20.

### S.4 The df-factor class: `SI20-P-*-cfxatt` (bounded, rtol 1e-1)

On `ext_cont_ext_cfx`, `ext_cont_ext_cfx_subgroups` and `ext_cov10_cfx_inference` the att-side cfx
bin matrices (`cfxatt_cfxpost`, `cfxatt_cfxpre`, `cfxatt_relyr*`, `cfxatt_cfxpost_sg1..3`) fire a
bounded row: the estimate column is unchanged, se and t move by sqrt(df_old/df_new) (measured
3.94e-4 / 3.78e-4, i.e. df_r 1268 -> 1267 and 1322 -> 1321), and only p and near-zero ci cells
tail-amplify (worst 2.36e-02 on `cfxatt_relyrn2` where p is 2.6e-14; 8.94e-03 on `cfxatt_cfxpre`).
The bound of 1e-1 exists "so a real att-bin regression cannot hide". These are the same df -1 that
S.1 records on every plain-basic + lattice config (`main:att` se ratio 1.000394555 =
sqrt(1268/1267) in the reproduction); the ame-side twins `cfxame_*` are the O(1) fix-1 movers under
unbounded L rows. **Classification: unclear** — a bounded, mechanical side effect of the old-bug-fixed fix 1
(the estimate column is untouched), but whether the extra absorbed slope should count against df_r
is the open judgment call (S.1, open question 1), so these rows cannot be called right or wrong
until it is ruled.

---

#### Relevant but untested here
The reproduction runs the reghdfe path only (no poisson, no bootstrap/permute, no subgroups, no
cfx plots); those surfaces are covered by the suite dumps cited above and by the ship record's
grids, not re-run. The `0bn.`-vs-`0bno.` column *name* is not used as the detector anywhere in this
section (the Findings block shows it can read `0bn.` with b = 0); the *value* is.


---


# Part 3 - the continuous-treatment and subgroup degeneracy family


> ## Headline of record and evidence provenance (added 2026-09-13; identical in all six sections)
>
> **Superseding note added 2026-09-14.** Two things below have since moved and are corrected
> here rather than rewritten in place, because the block is accurate about the run it describes.
> (i) `4280f3b` is no longer `HEAD`; the report's own header names the run of record by its
> `prov[new]` sha, and **that** is the build the report's claims attach to. (ii) The closing run
> of this assessment is the **2026-09-14** one over a **559-config** universe against the same
> v1.3 old side, reported in Part 7.5: **559 of 559 compared, zero DIFFERS, unexplained 0,
> stale_ledger 0**. The 2026-09-11 run described below is where this section's EVIDENCE was
> gathered and its magnitudes measured; the 2026-09-14 run is what verifies that nothing still
> differs. Neither supersedes the other - they answer different questions.
>
> **Run of record.** The 2026-09-11 dual-build run: **old = `301b706`** (published v1.3),
> **new = `14c30d3`**, whose `wooldid.ado` is byte-identical to this section's pin `08ec0b4`
> and to current `HEAD` apart from a nine-line comment block and a shortened `dis as text`
> permutation caution (10 insertions, 1 deletion; no estimation code).
>
> **The report judges those dumps against the ledger and comparer at `4280f3b` = `HEAD`.** That
> suite defines a **538-config universe**; the 2026-09-11 dumps cover **507** of them:
> **507 compared — EXPLAINED = 374, MATCH = 25, RC_EXPLAINED = 108 — unexplained 0,
> stale_ledger 0, unscoped 0, prov_fail 0, COVERAGE_OK.** The other 31 configs are new at
> `4280f3b`, have no 2026-09-11 dumps, and are outside every claim in this report.
>
> An earlier headline — **507/507 compared, EXPLAINED = 399, RC_EXPLAINED = 108**, at ledger
> `26ce421` — is quoted by three of the six drafts (`legacy` and `intended` give the full
> 399/108 line; `inference` gives the 507/507). It is the **same dumps and the same 507 configs**
> judged by the pre-`4280f3b` comparer. The two reconcile exactly: `399 = 374 + 25`, the 25 being
> the configs on which MATCH became reachable once the comparer began comparing across batch 3's
> macro→scalar kind change (`4280f3b`, verbatim: *"It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)"*; the same paragraph
> puts the cells affected at *"361 values on 103 configs"*). **Both numbers are measured and
> neither is withdrawn.** See `LEDGER_DELTA.md` §0.
>
> **Evidence provenance.** Every `test/suite/out/run_old|run_new/<config>.csv` cited below is a
> **regenerable** artifact, not a stored one: `test/suite/out/` is gitignored (`test/.gitignore`
> line 5) and the 2026-09-11 dumps were overwritten in place by the 2026-09-13 run. Regenerate
> them with `./run_full.sh oldnew` from `test/suite`, which stages old from
> `git show main:wooldidinstall/wooldid.ado` and new from `git show HEAD:…`. `compare_verbose_0911.txt`
> — cited by four sections, twice as `scratchpad/compare_verbose_0911.txt`, which is **not**
> repo-relative and never was — is likewise a recipe: rebuild it with `compare_oldnew.py … --verbose`,
> passing `--ledger`/`--configs` from whichever commit's judgement you want. The command, the build
> table and the one harness caveat (`4280f3b` changed `mkrun.py`'s `dumpres`, so a regenerated dump
> may differ in its `files|` line or gain a `matnames|…|UNREADABLE` line — **value cells are
> unaffected**) are in **`review_process/assessment_2026-09-13/repro/README.md`**.
>
> **Reproductions.** The do-files this section cites "in this directory" are now parameterised —
> repo root and scratch workdir as positional arguments — and live at
> **`review_process/assessment_2026-09-13/repro/<section>/`**. Their exact invocations are in each
> file's header and in that README. Bulk artifacts (dumps, logs, copied `.ado` builds) are not
> copied: they are regenerable, and the README says how.
>
> **Per-edit log:** `EDITS.md`. **Ledger delta and closed open questions:** `LEDGER_DELTA.md`.


Drafted 2026-09-12/13 for the old-wooldid bug attribution report of the 2026-09 full assessment
run. Builds: **old** = published v1.3 (`git main`, `301b706`, `old/wooldid.ado`, 3,450 lines)
vs **new** = HEAD `08ec0b4` (`new/wooldid.ado`). Sources of record: ledger
`test/suite/ledger.csv` (row ids cited), the verbose re-judge of the 2026-09-11 dumps
(regenerable; the filename `scratchpad/compare_verbose_0911.txt` is a session path, **not**
repo-relative — see the headline block above), the per-config `%21x` dumps in
`test/suite/out/run_old` / `run_new`, `wooldid_solved_issues.md`, `wooldid_standing_issues.md`
section A, `review_process/bug_sweep_2026-09-06.md`, `review_process/rulings_2026-09-07.md`
(section J carries the 2026-09-08 rulings), `review_process/si_reports/si7_omission_marker_report.md`,
`review_process/cell_collapse_subgroups_2026-09-01.md`, `review_process/qr_ship_2026-09-06.md`.
Old-code line numbers below were located in `old/wooldid.ado` (v1.3) for this report; new-code
lines in `new/wooldid.ado` (HEAD).

**Reproduction design (all sections).** `repro_run.sh`, `repro_old.do`, `repro_new.do`,
`repro_helpers.doh`, `repro_summarize.py`, `dumpparse.py` and `gen_repro.py` are at
`review_process/assessment_2026-09-13/repro/contsg/`, parameterised on repo root and scratch
workdir. `repro_run.sh` runs `repro_old.do` and `repro_new.do`
in two SEPARATE Stata processes (each `adopath ++` its own build directory), mirroring
`test/suite/mkrun.py`'s per-config block (seed/sortseed 98765, `ereturn clear`, the config's
prep lines); each block prints tagged values (`<ID>_<quantity>=` in `%21x` and decimal) and a
`REPRO_<ID>_DONE` sentinel, and dumps the run with the tested harness. `repro_summarize.py`
re-joins Stata's 80-column wraps and writes `repro_result.txt`; its last section runs
`test/compare.py` at `--rtol 0` between the reproduction's dumps and the suite's 2026-09-11
dumps: **every value cell is bit-identical on all five reproduced configs on both builds**
(the only difference groups are `matnames`, because the harness's `dumpres` writes rownames
without their equation prefix). All nine sentinels fired on both legs. The new-build do-file
also runs two extra `fullprecisionslow` legs (`QRFPS` on cont_sg, `SI7FPS` on
ext_contract10_rcs_factor_fe) that turn out to be decisive for S4 and S6. Data classes
(`test/suite/data_classes.csv`): dta_cont wellbehaved; survey_rcs wellbehaved;
dta_path_relyrgap stress; cont_sg_tidose stress; trig_d3void pathological (purpose-built d3
pin).

Family map: S1 AMESE-R1 · S2 AMESE-R4-l/-r + R4W-fitstats · S3 SIA4 (SI-A4/R6) · S4 SI7 ·
S5 BS8 · S6 QR-sgjt · S7 C-SGlab (dupcol) · S8 A17-sgv (e(V)-preserve).

---

### S1. Cohort-clustered continuous-treatment ame inference was a structurally invalid plug-in (AMESE R1/R3, family A)

Ledger rows: `AMESE-R1-w`, `AMESE-R1-v`, `AMESE-R1-sgp`, `AMESE-R1-sgr`, `AMESE-R1-jc`,
`AMESE-R1-ja`, `AMESE-R1-jp`, `AMESE-R1-jm` (and the siblings `AMESE-R1-w2`/`-w3` on
plot_cont/trig_cfx, outside this group's configs).

- **Symptom / impact:** with `contreat()` and a cluster-type VCE in which a treated
  cohort-by-period cell lies wholly inside one cluster (`cluster(i)` is the canonical case),
  old reported the delta-method ame standard errors that `margins` computes from the posted
  cluster covariance — a covariance that cannot see within-cell dose variation at all. The
  standing-issues A diagnosis (26-case x 256-rep MC grid, two locked theory reviews, a blind
  theorist, a referee) measured the G1 signature: `e(V_es_ame)` of exact rank 1, **every ame SE
  equal to the SE of one pooled slope**, SE/truth ratios ~0.20 (ES horizons) to ~0.52 (main),
  silent. Old also posted ame joint pre-trend p's and pairwise-subgroup ame contrasts on that
  covariance — the contrast SE being cancellation residue. On the suite's own well-behaved
  continuous panel the user saw a main ame with **t = 2285.6, p = 3.0e-27**, a subgroup ame
  contrast with **SE 1.8e-10, t = 1.3e7, p = 4.2e-61**, and a joint ame pre-trend p of .0031.
  Old wooldid shares every pathology of that grid (same spec/VCE paths); the standing record
  says so explicitly (`wooldid_standing_issues.md` A.1, last paragraph).
- **Mechanism (v1.3):** line 1809, inside `if "`contreat'" != ""`: `qui margins [w],
  dydx(contreat) subpop(... & t >= ttre) noestimcheck ...` followed by
  `matrix ___C`item'' = r(table)[1..6,1]'` — margins' se/t/p/ci from the posted cluster VCE are
  copied verbatim into the ame row; the pre-treatment twin at 1839; the per-subgroup twins at
  2496/2502 (`r(table)[1..6,1..colsof(r(b))]'`); the subgroup horizon twin at 2516; the ame
  joint pre-trend test at 2149 (`local joint_Choriz = r(p)[1,1]`); the pairwise ame contrasts
  named at ~2921. There is no containment screen anywhere in v1.3 (`grep -c` of "within a
  single cluster" and of "all-zero variance" in `old/wooldid.ado`: 0 and 0).
- **Fix + proof:** `7054878` (2026-08-24; R1 + R3, the `wooldid_amescreen` subprogram — new
  line 12363 — sthlp, battery re-bless, cert suite), `ea69544` (08-25, refit narrowing),
  `294513a` (08-25, V1-wave closing). Author-ruled (R2, the replacement estimator, skipped);
  design contract `review_process/amese/R1_design_note.md` committed before implementation.
  Behaviour: if ANY treated (i, relyr) cell is wholly inside one cluster of `cluster()`
  ("any", not "every" — the referee demonstrated a partial-containment false negative with a
  spurious p = 1.2e-08), print the mechanism + three valid alternatives, scrub every ame-family
  se/t/p/ci on all surfaces (incl. cfx and plot bands), withhold `e(V_es_ame*)` while posting
  `e(b_es_ame*)`, skip the ame V_es tests; never under fastnose/impvalues/unconditionalse;
  att side untouched by design. Verification (solved-issues record): 10-scenario cert
  `test_amescreen` (incl. the referee's partial-containment acid test), battery 34/34 with
  exactly the 8 predicted contreat+cluster(i) configs re-blessed, THREE full 450-config
  oldnew cycles to 0 unexplained, a 3-verifier wave (V2 hunk audit, V1 adversarial incl. the
  cfx-leak find, V3 contract with all 26 ledger rows fired).
- **This run:** the genuine R1 signature among this group's configs is **cont_sg** only:
  `wdidmainresults` 30 cells (`AMESE-R1-w`: every ame se/t/p/ci, e.g. [3,2] .000437716 → .,
  [3,3] 2285.57 → ., [3,4] 2.99e-27 → .), `V` 18 cells (`AMESE-R1-v`: ame diagonals
  1.9159545e-07 → the −1e100 sentinel), `pwsgs` 6 cells (`-sgp`: ame contrast se 1.76e-10, t
  1.32e7, p 4.18e-61 → .), `wdidsgresults_1/_2` 30 cells each (`-sgr`: per-subgroup ame rows;
  t 2282.39 / 2287.70 → .), `joint_Choriz` .00308 REMOVED (`-jc`), `joint_ame_allsgssame` and
  `joint_pre_ame_allsgssame` REMOVED (`-ja`/`-jp`; old stored them as missing —
  **cause not diagnosed**. *Corrected 2026-09-13:* this used to read "the T-M2 `[1,2]`
  read defect". It is not that. The two v1.3 sites that produce these scalars are 2651
  and 2663, and both read `r(p)[1,1]` — verbatim, `local joint_C`item''_sdiff = r(p)[1,1]`
  and `cap local prejoint_C`item''_sdiff = r(p)[1,1]`, each off the
  `contrast(overjoint)` margins call immediately above it (2650 / 2662). T-M2 is a
  different pair of sites — `local joint_Choriz_s`sgcase'' = r(p)[1,2]` at 2620 and 2631,
  located precisely by `inference` §4 and fixed at new 3235/3246. Following `legacy`
  L.8 item 6, and the ledger row `AMESE-R1-ja` whose own text calls this *"the
  allsgssame-ame open observation"*, the old missing-valued continuous `allsgssame` on a
  NON-containment, full-rank design is **left as an open question rather than asserted
  closed**), `joint_Choriz_s1/s2` macros REMOVED (`-jm`). Att-side cells agree to
  ≤1e-12 (main att se .001425182992 on both; pwsgs att se .000027529347 on both).
  `AMESE-R1-w`'s `cont_*` pattern also matched `wdidmainresults` on **cont_lat** (6 cells: the
  SI-20 ame point shift 3.0043 → 1.0000064, SEs missing on BOTH builds — not an R1 scrub) and
  on **cont_sg_gap / cont_sg_tidose** (6 cells: 0 → ., the R6 void of S3); the siblings
  `-w2` (plot_cont, 7 cells) and `-w3` (trig_cfx, 14 cells) are likewise all SI-20 point
  shifts (3.0026 → 0.9980) with **no** "old-missing → recovered att cells" — consistent with
  the R4 refit having been withdrawn (S2).
- **Reproduction:** `repro_run.sh`, block `AMESE_R1` (config cont_sg: `wooldid y i t ttre,
  contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2)
  subgroups(sg) jointtests safety(off)` on `test/data/dta_cont.dta`, `sg = mod(i,2)+1`).
  Decisive lines from `repro_result.txt`:
  ```
  AMESE_R1_main_ame_est      OLD 1.00043254097          NEW 1.00043254097
  AMESE_R1_main_ame_se       OLD .00043771617           NEW .
  AMESE_R1_main_ame_t        OLD 2285.57364476          NEW .
  AMESE_R1_main_ame_p        OLD 2.9916501e-27          NEW .
  AMESE_R1_ESame_pre2_se     OLD .00043771617           NEW .      <- identical to the main
  AMESE_R1_ESame_post1_se    OLD .00043771617           NEW .         ame SE: the rank-1 signature
  AMESE_R1_V_3_3             OLD 1.9159545e-07          NEW -1.000000e+100
  AMESE_R1_V_1_1 (att)       OLD 2.0311466e-06          NEW 2.0311466e-06
  AMESE_R1_main_att_se       OLD .001425182992          NEW .001425182992
  AMESE_R1_pwsgs_ame_se      OLD 1.7613228e-10          NEW .
  AMESE_R1_pwsgs_ame_t       OLD 13202020.7819          NEW .
  AMESE_R1_pwsgs_ame_p       OLD 4.1793344e-61          NEW .
  AMESE_R1_joint_Choriz      OLD .003080477038          NEW ABSENT
  NEW log: "Caution: with standard errors clustered on i, 40 treated i-t cell(s) lie completely
  within a single cluster, preventing valid analytic estimation of SEs for AME estimates ..."
  and "standard errors were not identified for 18 of 36 reported estimates."
  ```
- **Classification: old-bug-fixed.**
- **Open questions:** (i) the `cont_*` pattern of `AMESE-R1-w` absorbs two other causes on
  this group's configs (SI-20's point shift on cont_lat; the R6 void on cont_sg_gap/tidose) —
  attribution there is by pattern, not by cause; a per-config split or a comment on the row
  would make the ledger self-explaining. (ii) The att-side G1a finding (cohort-clustered att
  SEs ~25-30% understated) is unaddressed on both builds by design and produces no diff.

---

### S2. The AMESE R4 refit and its withdrawal: what the `AMESE-R4-l/-r` and `R4W-fitstats` rows actually record

Ledger rows: `AMESE-R4-l-cont_lat`, `AMESE-R4-r-cont_lat`, `R4W-fitstats-lat`,
`R4W-fitstats-lat2` (siblings on plot_cont, trig_cfx, si20_lat_cfx_pw outside this group's
configs).

- **What the rows record (not an old bug):** on cont_lat five model-fit statistics differ
  old → new while every estimate, standard error and test agrees: `ll_0` −2419.2233 →
  −2315.5552 (4.3e-2), `r2_within` .9999901132 → .9999878353 (2.3e-6), `r2_withinadj`
  (2.9e-6), `r2adj` .99999836045 → .99999835827 (2.2e-9), `rmse` .0098460017 → .0098525374
  (**6.6e-4**). `ll` agrees to 1e-14 (3342.93274585), `rss` to 5e-14 (.073095587313), `df_r`
  = 9 on both, `e(properties)` = `b` on both (all-zero VCE on this design on both builds).
- **Mechanism — reghdfe's residual-dof accounting on the SI-20 model, not the refit:** the
  identity `rss / rmse²` recovers reghdfe's residual degrees of freedom exactly from the dumps:
  **754.000000 (old) vs 753.000000 (new)**. New consumes one more absorbed dof because SI-20
  (`1df8613`, 2026-09-05: "basic always absorbs the baseline dose control and, with a lattice,
  its lattice twin") moved the ringer term into the absorb set; `r2adj = 1 − (1−r2)(N−1)/df`
  then moves by exactly the observed 2.2e-9, and `r2_within`/`ll_0` change because the
  within-transformation and the null model now include the absorbed term. The SI-20 record
  itself lists "df_a over-count from the absorbed ringer terms" as a judgment call left to the
  author. Old's statistics are reghdfe's on the v1.3 model (where that control entered as an
  explicit regressor — the SI-20 point defect that shifted every ame-family row by the baseline
  dose slope, 3.0043 → 1.0000, row `SI20-L-cont_lat-b`, another group).
- **Why the ledger attributes them to R4, and why that is stale:** R4 (`7054878`, 08-24) was
  the "blackout repair": detect Stata's ereturn-post V-zeroing after reghdfe and re-estimate
  via `regress` on explicit dummies to salvage the att standard errors. Between `7054878` and
  `b801282` the new build's stored fit statistics on cont_lat therefore came from that
  `regress` refit, whose N − k dof is invariant to absorbed-vs-explicit: the baseline history
  of `test/baseline/cont_lat.csv` shows `7054878` changing `r2adj` by one ulp, dropping
  `r2_within` (regress does not report it) and moving `ll_0`, with **`rmse` untouched** (the
  refit reproduced old's rmse bit-exactly); SI-20's re-bless `5aabe1e` did not touch these
  lines; and **`b801282` (2026-09-08, the withdrawal) is the only commit that moves `rmse`
  (+1.42a23fadd0474X-007 → +1.42d912f236402X-007), `r2adj`, restores `r2_within` and flips
  `properties` `b V` → `b`**. The withdrawal thus EXPOSED reghdfe's own accounting on the SI-20
  model; the `R4W-fitstats` rows were written from that battery diff with "r2adj and rmse move
  at 1e-8 relative" — true for r2adj (2.2e-9), wrong for rmse (6.6e-4). The `AMESE-R4-l/-r`
  reasons ("ll/ll_0 are regress-vs-reghdfe accounting"; "within-R2 not reported by regress
  (removed)") describe a refit that no longer exists: `r2_within` is present on both builds
  this run and differs at 2.3e-6.
- **The withdrawal itself (author ruling 2026-09-08, `rulings_2026-09-07.md` §J; commits
  `b801282`, `3743aec`; record `review_process/research_2026-09-08/fastnose_precision.md`):**
  the refit "salvages a subset of standard errors from a specification whose containment has
  already invalidated the ame inference, at the cost of a solve over hundreds of collinear
  dummy columns, which delays the user's discovery that the specification is broken." Ruled
  consequences: estimates move at 1e-12 relative (the partialled fit is better conditioned —
  here `rss` differs at 5e-14 and `main_att` at 1e-15), and a family whose covariance is
  entirely non-identified now posts its point estimates and pre-horizon count under a caution
  instead of withholding them with the covariance. Old-vs-new on the posting: old posted b-only
  **silently** (line 2987 `cap ereturn post eb ev` fails on the missing diagonal, 2989 posts
  `eb` alone; v1.3 has no all-zero-VCV screen at all); new posts b-only under three cautions
  ("the estimation command returned an all-zero variance-covariance matrix ..."; "standard
  errors were not identified for 12 of 12 reported estimates, so no ... e(V) was posted";
  "no horizon of the att event study has an identified variance, so e(b_es_att) posts without
  e(V_es_att)"). Same `e(properties)`, so no ledger row fires for it; the difference is console
  only.
- **This run:** cont_lat (values above); plot_cont and trig_cfx: `ll_0` −2419.4328 →
  −2316.7228, `r2_within` .9999898913 → .9999875862, `r2adj` 2.2e-9, `rmse` .0099579556 →
  .0099645656 (6.6e-4); si20_lat_cfx_pw: `rmse` .0098460017 → .0098525374.
- **Reproduction:** `repro_run.sh`, block `R4W` (config cont_lat: `wooldid y i t ttre [aw=w],
  contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2)
  espostlength(2)` on dta_cont). Decisive lines:
  ```
  R4W_rss          OLD +1.2b66474fd929aX-004 (.073095587313)  NEW +1.2b66474fd93a5X-004 (.073095587313)
  R4W_rmse         OLD .009846001714                          NEW .009852537395      -> rss/rmse^2 = 754.000000 vs 753.000000
  R4W_r2adj        OLD .999998360451                          NEW .999998358273
  R4W_r2_within    OLD .99999011322                           NEW .999987835345       (present on both)
  R4W_ll_0         OLD -2419.22327162                         NEW -2315.55522277
  R4W_ll           OLD 3342.93274585                          NEW 3342.93274585
  R4W_df_r         OLD 9                                      NEW 9
  R4W_properties   OLD [b]                                    NEW [b]
  R4W_main_att     OLD 10.0097918111                          NEW 10.0097918111
  R4W_main_att_se  OLD .                                      NEW .
  R4W_main_ame     OLD 3.00432338576                          NEW 1.00000636895       (SI-20, other group)
  ```
- **Classification: unclear** (revised 2026-09-13; see `RECONCILIATION.md` §S-2). This section
  first read the family as *intended-change* on the R4 withdrawal, then wrote, correctly, that
  "old-vs-new here is not an old bug and not the withdrawal of one: v1.3 never had the refit."
  A ruling that withdraws something v1.3 never had cannot be what an old-vs-new difference
  records. What the difference records is reghdfe's residual-dof accounting on the (fixed) SI-20
  model — 754 old against 753 new, the new build consuming one more absorbed dof — and THAT is
  unruled. It is item 1 of `review_process/held_for_ruling_2026-09-13.md` and `si20` S.1's open
  question 1 confirms no ruling exists in `qr_ship_2026-09-06.md`,
  `held_for_ruling_2026-09-09.md` or the rulings memory. `si20` S.3 and S.4 call the identical
  mechanism *unclear*; all three now agree and will move together on one ruling. No estimate,
  standard error or test moves either way.
- **Open questions:** (i) **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E2]** — all four were rewritten by `4280f3b`;
  `R4W-fitstats-lat2` now reads *"~6.6e-4 relative - not 1e-8 - measured 2026-09-11"* and
  `AMESE-R4-l/-r-cont_lat` now say the refit was removed by `b801282`. The question as written
  was: the four reason texts are stale or inaccurate (`R4W-fitstats-lat2`
  understates rmse by five orders of magnitude; `AMESE-R4-l/-r` describe the removed refit);
  (ii) whether reghdfe's extra absorbed dof on the SI-20 model is a df_a over-count is the
  judgment call the SI-20 record already leaves to the author — if it is, new's `rmse` is
  6.6e-4 too large on lattice+basic runs (a fit statistic only).

---

### S3. The d3 class: old presented fabricated structural zeros as ame point estimates (SI-A4 / R6)

Ledger rows: `SIA4-d3`, `SIA4-d3p` (trig_d3void), `SIA4-gp-p`, `SIA4-gp-s` (cont_sg_gap),
`SIA4-td-p`, `SIA4-td-s` (cont_sg_tidose).

- **Symptom / impact:** when the dose is constant within every treated cohort-by-period cell
  (standing-issues A group G7, "d3": e.g. a dose that is a function of the unit id, or
  time-invariant per unit on the subgroup fixtures), every treated-cell dose-slope term is
  dropped for collinearity and the within-cell ame is genuinely unidentified — bucket (iv), no
  option rescues it. Old showed the ame rows of `e(wdidmainresults)`, `e(wdidsgresults_*)` and
  the ame contrasts of `e(pwsgs)` with **point estimate exactly 0 and missing SE**, behind the
  generic collinearity-count cautions, and posted b-only; a user reading the table sees an
  average marginal effect of zero, not "not identified".
- **Mechanism (v1.3):** the same line 1809 as S1: `margins ... noestimcheck` computes the
  margin anyway when every slope column is omitted — `r(table)` then carries b = 0 with se
  missing — and `matrix ___C`item'' = r(table)[1..6,1]'` files the 0 as the estimate; the
  pre-treatment twin 1839, the per-subgroup twins 2496/2502, the subgroup horizons 2516, the
  pairwise contrasts. Then 2976 `ev = diag(hadamard(se,se))` carries missing, 2987 fails and
  2989 posts b-only (S8).
- **Fix + proof:** R6, author-ruled ("named refusal of the ame block ... instead of
  structural-zero rows behind generic cautions", standing A.3), landed in the SI-A4/R6 closure
  wave, final `38c49ae` (2026-08-25; full suite 455/455, 0 unexplained, cert 24/24, battery
  identical). The landing took one recorded fix-up: voiding the ame POINT column reached
  `e(b)`'s construction and `ereturn post` refuses missing coefficients, so `e2a8d74` lets the
  voided entries re-enter `e(b)` as structural zeros — the one raw surface out of reach, stated
  in the caution — while every wooldid-native surface keeps them missing (new lines 1628/3872,
  `amevoid`). The same run found cont_sg_gap and cont_sg_tidose to be LIVE d3 members; their
  battery baselines were re-blessed on evidence after the diffs were confirmed to be exactly the
  void class.
- **This run:** trig_d3void: `wdidmainresults` 6 cells ([3,1], [4,1], [9,1]..[12,1]: 0 → .),
  `properties` `b` → `b V` (`SIA4-d3p`, the S8 sentinel machinery; att rows bit-identical).
  cont_sg_gap / cont_sg_tidose: `wdidsgresults_0/1/2` 29 cells each (row 3 `main:ame` 0 → .;
  rows 8-11, the four ES_ame rows, ABSENT — that part is BS-8, S5); `pwsgs` 36 cells each, of
  which **6** are R6 (the ame-contrast estimates 0 → .) and **30** are the att-row se/t/p/ci
  cells blanked by BS-15 (`cfac3f5`, see open question).
- **Reproduction:** `repro_run.sh`, block `SIA4` (config trig_d3void, the inline 12-unit panel
  with `dose = 0.5 + 0.1*i`: `wooldid y i t ttre, att robust contreat(dose)
  contreatcontrols(basic) esprelength(2) espostlength(2)`), plus block `SGLAB` for the
  subgroup twins. Decisive lines:
  ```
  SIA4_main_att_est   OLD .275664260704   NEW .275664260704     (att bit-identical)
  SIA4_main_att_se    OLD .046901528538   NEW .046901528538
  SIA4_main_ame_est   OLD 0               NEW .                 <- the fabricated zero
  SIA4_main_ame_se    OLD .               NEW .
  SIA4_ESame_pre2_est OLD 0               NEW .
  SIA4_ESame_post1_est OLD 0              NEW .
  SIA4_b_main_ame     OLD 0               NEW 0                 (e(b): the e2a8d74 seam, on both)
  SIA4_properties     OLD [b]             NEW [b V]
  SIA4_V_1_1          OLD (no V)          NEW .002199753379     (= .046901528538^2)
  SIA4_V_3_3          OLD (no V)          NEW -1.000000e+100
  NEW log: "Caution: every treated cohort-period dose-slope term for the continuous treatment
  variable dose was dropped for multicollinearity ...; this occurs, for example, when dose is
  constant within each cohort. The within-cell average marginal effects (ame) that wooldid
  reports for continuous treatments are NOT IDENTIFIED in this design ..." and "standard
  errors were not identified for 6 of 12 reported estimates."
  SGLAB_main_ame_est  OLD 0               NEW .        SGLAB_sg0_r3_est OLD 0  NEW .
  SGLAB_pwsgs_ame_est OLD 0               NEW .        SGLAB_pwsgs_ame_se OLD . NEW .
  ```
- **Classification: old-bug-fixed.**
- **Open question:** on cont_sg_gap/cont_sg_tidose the `SIA4-gp-p`/`SIA4-td-p` rows absorbed
  the 30 att-row `pwsgs` cells that BS-15 blanks (repro: `SGLAB_pwsgs_att_se` .021559719779 → .,
  `_t` −14.98 → ., `_p` 3.47e-15 → .; new caution "rank-deficient in the contrast direction
  for 6 pairwise row(s) (largest implied correlation ... 0.9990)"; `e(wdid_sgdegen)` = 6). The
  ledger has `BS15-jt`/`BS15-jtp` rows for both configs but no `BS15-pw` row, so the SIA4
  pattern explains cells it does not describe; adding `BS15-pw-cont_sg_gap/tidose` rows would
  close the gap.

---

### S4. The SI-7 rows on ext_contract10_rcs_factor_fe: a real old defect, but this run's cells are not the SI-7 drop

Ledger rows: `SI7-sg-m`, `SI7-sg-n`, `SI7-eb-n`, `SI7-eb-m`.

- **Symptom / impact (old):** on the repeated-cross-section survey design with `fe(i.group##i.t)`
  absorbing the subgroup-by-time variation, `cluster(state)` (24 clusters) and `subgroups(group)`,
  old returned rc 0 with a table the SI-7 diagnosis calls "quietly degenerate": an all-zero
  clustered VCV (every se/t/p/ci in every returned matrix missing, with **no comment** — v1.3 has
  no all-zero-VCV screen), pre-trend "estimates" that are machine zeros (8.88e-16) computed
  under protest, and **subgroup 1's pretreatment row silently absent** (`e(wdidsgresults_1)` 3
  rows vs subgroup 2's 4; `e(b)` without `pretreatment_SG1_att`). Where Stata did flag the
  non-estimable horizon, the pre-fix new build died r(198) "5o. invalid name" after estimation
  had succeeded — the SI-7 ticket.
- **Mechanism (v1.3):** the subgroup event-study parse at 2285-2300 reads relative periods off
  `r(table)` colnames and strips `bn`, the `_t` tempvar and `1.targetflag` — an omission-marked
  `5bno.` becomes `5o.` and `local relyrinj = 5o. -6 -1` (2304) throws; an unmarked machine-zero
  column parses as a number and is reported. The pretreatment-row match at 2815-2829 compares
  the stripped `over()` names with the subgroup value and sets `sghaspre no` on any mismatch,
  after which 2833-2834 skip the row and 2867-2874 drop its label with no message. The all-zero
  VCV is reghdfe's cluster VCE zeroed by Stata's ereturn-post screen on this design — on both
  builds (si7 report §(b): `sum(e(V))` = +0.0 bit-exact on both).
- **Fix + proof:** `c014be7` (2026-08-27, author-ruled "proceed, do not refuse"): READ the
  marker, drop the non-estimable horizons from that subgroup's event study with a caution, skip
  a subgroup with none; `dd4e85f` (the shared `wooldid_sglbl` subprogram, byte-identical);
  `20c431d` (BS-26: fast mode keeps the drop); `49bfb52` (BS-8: the AME twin, S5). The new
  all-zero-VCV screen prints the caution and "12 of 12 not identified". Cert `test/test_si7.do`
  (T1 firing, T2 dormant; T1 needs `fullprecisionslow` since route P — `qr_ship_2026-09-06.md`
  §5). Suite pin rc_new 198 → 0 with the five SI7 rows (`7f704d3`).
- **This run — the rows fired on a different diff than they describe:** the reasons say the
  non-estimable sg1 pre2 horizon is dropped (e(b) 11 → 10 names). Observed: e(b) **1x11 → 1x12**
  — new ADDS `pretreatment_SG1_att` (= bitwise 0, relyr .), `wdidsgresults_1` 3x7 → 4x7 —
  while `ES_att_SG1_pre2` = 8.8817842e-16 is present on BOTH builds bit-identically. The 7
  `SI7-eb-m` cells are the positional shift from the inserted entry; the 11 `SI7-sg-m` cells in
  `wdidsgresults_1` are the inserted row, plus `wdidsgresults_2` [3,1] 8.88e-16 → 1.78e-15
  (machine-zero residue under the collapse's summation order). Why no drop on new: route P
  (`4b224f0`, 2026-09-02; collapse default since `111ecbb`) runs the cell collapse with
  `subgroups()` behind a per-read degeneracy screen, and under `safety(off)` "there is no
  fallback — the run proceeds on the collapsed rows and prints a caution"
  (`e(wdid_fcguard)` = caution); on the collapsed rows the sg1 pre2 margin is not
  omission-marked, so the machine zero survives as an estimate on new as it did on old. The
  `SI7FPS` leg (new + `fullprecisionslow`, same process) shows exactly the behaviour the rows
  describe: sg1 = [att att post0], `ES_att_SG1_pre2` ABSENT, "Caution: the relative period
  effect(s) pre2 are not estimable for subgroup 1 (margins reports them as omitted), so they
  are dropped from that subgroup's event study".
- **Reproduction:** `repro_run.sh`, block `SI7` (config ext_contract10_rcs_factor_fe:
  `wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) fe(i.group##i.t)
  esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) safety(off)` on
  `test/suite/data/survey_rcs.dta` after `replace y_rcs = y_rcs - rcs_cell_shock`), plus the
  new-only block `SI7FPS` (same call + `fullprecisionslow`). Decisive lines:
  ```
  SI7_b_pretreatment_SG1_att  OLD ABSENT (ncols 11)          NEW 0 (firstpos 6 of 12)
  SI7_sg1_shape               OLD 3x7 [att pre2 post0]       NEW 4x7 [att att pre2 post0]
  SI7_sg1_r2_est / relyr      OLD 8.8817842e-16 / -2         NEW 0 / .
  SI7_sg1_r3_est              OLD .6976081812 (post0)        NEW 8.8817842e-16 (pre2)
  SI7_b_ES_att_SG1_pre2       OLD 8.8817842e-16              NEW 8.8817842e-16   (retained on both)
  SI7_main_att_se (sg1)       OLD .                          NEW .
  SI7_properties              OLD [b]                        NEW [b]
  SI7_wdid_fcguard            OLD ABSENT                     NEW [caution]
  SI7FPS_sg1_shape            NEW 3x7 [att att post0];  SI7FPS_b_ES_att_SG1_pre2 ABSENT;  SI7FPS_wdid_fcguard ABSENT
  NEW log (default path): 7 cautions - all-zero VCV (x2), "wooldid found subgroup margins with
  zero variance at a subgroup margins call ...", "12 of 12 not identified", and e(b_es_att),
  e(b_es_att_s1), e(b_es_att_s2) posting without their covariances.  OLD log: none.
  ```
- **Classification: unclear.** The family is a genuine old defect (silent row loss, no
  all-zero-VCV screen, machine-zero pre-trends reported without comment), and the SI-7 drop is
  proven on the `fullprecisionslow` leg — but the cells that differ this run are not produced
  by the SI-7 drop; the four reasons are stale for the run's default path since route P. The
  pretreatment-row match site is byte-identical on new (3500-3514 vs old **2812-2826**
  — *range corrected 2026-09-13*; `diff` of those two spans is empty, both beginning
  `local truerow `sgnum'` and ending `if "`truerow'" == "" local sghaspre no`), so the
  row's return is a numerical knife-edge on the margins column name (the si7 report §(d)
  attributes the flip to B21's float→double `w_att`), not a fix at that site; and the value new
  reports there is a bitwise 0 under the fcguard caution.
  **A second commit did act just outside that span, and this section previously asserted
  byte-identity across it.** Old 2829 reads `local trurerow` — a typo clearing a local that
  does not exist — and new 3517 reads `local truerow`, repaired by `01a3b58` (2026-07-22,
  *"Author edit: fix trurerow typo in the SG horizon-matrix assembly"*, the only commit
  besides the 2023 upload that touches the string: `git log -S'trurerow'`). The `else` leg
  is reached only when `round` ∈ {`pre`, `contpre`} **and** `sghaspre` != `yes`; in that
  state both append lines are gated off on either build (old 2834/2835, new 3522/3523), so
  the stale `truerow` is never used as a row subscript — but old 2836 / new 3528
  (`if "`truerow'" != "" local firsttimethrough 0`) is reached, and there old still holds
  `sgnum` while new holds empty. **`01a3b58`'s reach on this config was UNDETERMINED in the
  earlier draft; it is now MEASURED, and the answer is that it has none.** Four one-line-twin
  builds were run on this config, 2026-09-14 (working dir `critleft/`):

  ```
  build      wooldid.ado md5                    what it is
  old        5f067940742613346d90fc65eb088590   v1.3, git main = 301b706; L2829 `local trurerow`
  old+fix    5fddf7298a87fc6822fbe2fa758bbcb8   old with L2829 -> `local truerow`
  new        382c1e0ee3e1be37fa2e9b677a65fdcc   08ec0b4 (= 6a05127, the run's new leg); L3517 `local truerow`
  new+typo   3b52a608fe31a6e25e169bea6022fa75   new with L3517 -> `local trurerow`
  each pair differs by exactly one line (diff: 1 changed line, nothing else).

  e() dumps (%21x, this section's own dumpres helper), default path AND fullprecisionslow:
    old vs old+fix    si7_default.csv BYTE-IDENTICAL   si7_fps.csv BYTE-IDENTICAL (rc 198 both: v1.3 has no fullprecisionslow)
    new vs new+typo   si7_default.csv BYTE-IDENTICAL   si7_fps.csv BYTE-IDENTICAL (rc 0 both)
  ```

  **Mechanism, measured rather than reasoned.** A counter inserted in the `else` leg, and a second
  in the enclosing per-object loop body, show the leg is entered **0 times** on this config against
  **4** loop iterations — on both builds and on both paths (`TR_ELSE_COUNT=0 TR_ITER_COUNT=4`;
  `TRFPS_ELSE_COUNT=0 TRFPS_ITER_COUNT=4` on new). The instrumented builds' dumps are byte-identical
  to the un-instrumented ones, so the probe is inert. The line `01a3b58` changed **never executes on
  this config**: whenever the `pre`/`contpre` branch is taken here, `sghaspre` is `yes` and the
  `if` leg runs, not the `else`. So `01a3b58` is **not** part of this config's old-vs-new
  difference, and S7's attribution of the dupcol repair to `9d2e4df` alone stands *for this
  config*. What survives from the critique is the narrower point, which this section accepts: the
  site is **not** byte-identical across the builds, and a second commit acted inside the range the
  earlier draft named.

  **Tie to the run of record.** The probe's dumps match the suite's 2026-09-13 dual-build run cell
  for cell: **124 of 124** value cells on old and **218 of 218** on new are identical to
  `test/suite/out/run_old|run_new/ext_contract10_rcs_factor_fe.csv` (old `301b706`; new `6a05127`,
  whose ado is byte-identical to this section's pin `08ec0b4`, md5 `382c1e0e…`). The only four
  non-matching lines are `matnames` rows, where this section's helper prints bare rownames
  (`att att pre2 post0`) and the suite's prints the equation-prefixed stripe (`main:att
  pretreatment:att ES_att:pre2 ES_att:post0`) — a harness formatting difference, not a value. That
  also re-ties this section's SI-7 shape claims (old `wdidsgresults_1` 3x7, new 4x7) to a dump that
  exists today.
- **Open questions:** (i) **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E3]** — `4280f3b` rewrote all four SI-7 reasons for
  the default path; `SI7-eb-n` now reads *"e(b) gains the pretreatment_SG1_att entry (default
  cell-collapse path); 11 names -> 12 on 2026-09-11 (the reason once said 11 -> 10: a drop that
  the collapse path superseded)"*. The question as written was: rewrite the four reasons for the
  default path (or pin the config with
  `fullprecisionslow` so the rows describe what fires); (ii) should a bitwise-zero subgroup
  pretreatment estimate under `safety(off)` be reported with the caution (current) or dropped
  like the SI-7 horizons? — a consistency question for the author; (iii) the exact colname old
  saw when it lost the SG1 row was not bisected (the si7 report records the outcome, not the
  name).

---

### S5. Non-estimable continuous-treatment subgroup horizons: old carried twelve structural zeros in e(b) (BS-8)

Ledger rows: `BS8-cont-sg-gap-b`, `BS8-cont-sg-tidose-b`.

- **Symptom / impact:** on the two live d3 fixtures every continuous-treatment subgroup
  horizon is non-estimable (the subgroup's dose slopes are all dropped). Old reported them as
  four `ES_ame_SG#_{pre2,pre1,post0,post1}` rows per subgroup with estimate 0 and missing SE —
  **twelve structural zeros in `e(b)`** (`_b[ES_ame_SG0_post1]` = 0), the same fabricated-zero
  class as S3, sitting behind the dupcol label shift (S7) with no statement that they are not
  estimates. BS-8's headline was a new-build crash — a bare `r(503)` when the binary SI-7 drop
  shortened the relyr column the AME call borrowed — which had no old analogue (old never
  dropped anything); the old-vs-new content of these rows is the zeros.
- **Mechanism (v1.3):** 2516 `matrix ___`horizitem'_Choriz_s`sgcase'' = (r(table)[1..6,
  1..colsof(r(b))]', relyrcolumn_s`sgcase'')` copies the AME subgroup-horizon `r(table)` verbatim
  under `noestimcheck` (zeros for omitted margins); 2855/2859 stack them into the per-subgroup
  matrix; 2969-2970 carry the same rows into `e(b)`/`e(V)`.
- **Fix + proof:** `49bfb52` (2026-09-06; the per-call AME scan, introduced as `wooldid_sgamelbl` and at HEAD
  carried by the shared `wooldid_sglbl` — new line 12078, called for the continuous-treatment
  block at 3185 with `kind(continuous-treatment)`): the AME call is scanned, labelled and
  filtered on its OWN colnames, prints its own SI-7 caution, skips a
  subgroup whose AME horizons are all non-estimable; `3e03179` re-blessed cont_sg_gap and
  cont_sg_tidose ("every other cell is byte-identical; an isolation run of the pre-refactor
  BS-8 build produced the identical diff"); cert `test_si7` T8 (the sweep's fixture: rc 0 with
  both skip cautions) and T9/T10 (one AME row dropped, siblings untouched, no repeated label).
  The ship record flags "blank rows versus dropped rows is a display choice; SI-7 consistency
  argued for dropping" for the author (`qr_ship_2026-09-06.md` §6); no ruling on it was found in
  `rulings_2026-09-07.md`.
- **This run:** `e(b)` 45 → 33 columns on both configs. The 25 "differing cells" the comparer
  counts are positional (old's label-shifted 45-vector against new's 33-vector, e.g. [1,20] old
  0.0 — old's `ES_ame_SG0_pre1` slot — → new 1.3272 = `main_SG1_att`); the twelve dropped
  entries are the `ES_ame_SG{0,1,2}_*` zeros. The per-subgroup matrices lose rows 8-11 (the
  `SIA4-gp-s`/`-td-s` ABSENT cells of S3).
- **Reproduction:** `repro_run.sh`, block `SGLAB` (config cont_sg_gap: `wooldid y i t ttre,
  contreat(dose2) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2)
  subgroups(sgr)` on `test/data/dta_path_relyrgap.dta`, `dose2 = 0.5 + mod(i,5) + 0.2*mod(i,3)`,
  `sgr = mod(i,3)`). Decisive lines:
  ```
  SGLAB_b_ES_ame_SG2_post1  OLD 0 (count 4, firstpos 42: four repeats of the label)   NEW ABSENT (ncols 33)
  SGLAB_sg0_shape           OLD 11x7 [att ame pre2 pre1 post0 post1 pre2 pre1 post0 post1 post1]   NEW 7x7 [att att ame pre2 pre1 post0 post1]
  SGLAB_sg0_r8_est          OLD 0                                                       NEW (row absent)
  NEW log, per subgroup: "the continuous-treatment relative period effect(s) pre2 pre1 post0 post1
  are not estimable for subgroup 0 (margins reports them as omitted), so they are dropped ..." and
  "no continuous-treatment relative period effect is estimable for subgroup 0, so its
  continuous-treatment event study is skipped entirely."  (x3 subgroups)
  ```
- **Classification: old-bug-fixed** (the zeros are the R6 fabricated-zero class; the row drop is
  the ruled SI-7 convention).
- **Open question:** the blank-vs-dropped display choice awaits the author; if "blank" were
  chosen, `e(b)` would regain the twelve entries as missing-SE rows and these two rows would
  need re-pinning.

---

### S6. Per-subgroup joint pre-trend p under the cell collapse (QR-sgjt): last-place residue, proven by the fullprecisionslow leg

Ledger rows: `QR-sgjt-1` — **retired at `4280f3b` (2026-09-13)**, when the cross-kind compare put
its 4e-12 move inside `rtol`; only `QR-sgjt-2` survives, and its reason now records the retirement
(`LEDGER_DELTA.md` §A) — and `QR-sgjt-2` (cont_sg).

- **What moved:** the per-subgroup joint pre-trend p's `joint_horiz_s1`/`s2` (old macros
  .1340523654409293 / .4122428908131096; new scalars .134052365441 / .412242890813) at 4.0e-12
  and 1.4e-12 relative, and `joint_pre_att_allsgssame` (1.57586260706e-19 → 1.57586260678e-19,
  1.8e-10 relative on a ~1e-19 chi-square tail). Old line 2430 `local joint_horiz_s`sgcase' =
  r(p)[1,2]` stored these as ereturn LOCALS; CS-60 (`605d6d7`, 2026-09-09) makes them ereturn
  scalars like the overall families (new 4177) — which is why the verbose shows the macro
  REMOVED and a scalar ADDED.
- **Mechanism:** the cell collapse (`111ecbb` default 2026-09-01; route P `4b224f0` 2026-09-02
  extends it to `subgroups()`/`contreat()`) changes the order in which margins sums cell
  contributions; a quantity that is mathematically identical comes out with last-place residue
  (`cell_collapse_subgroups_2026-09-01.md` §1). Deterministic across reruns; blessed in the
  baseline 2026-09-05.
- **Proof (decisive):** the `QRFPS` leg (new + `fullprecisionslow`, same Stata process as the
  collapsed leg) returns `joint_horiz` = `+1.b0445c643364aX-001` and `joint_pre_att_allsgssame`
  = `+1.74170cb89d771X-03f` — **bit-identical to old** — where the collapsed leg gives
  `+1.b0445c6430cc9X-001` and `+1.74170cb78182fX-03f`; `joint_horiz_s1` fps
  `+1.128a0bec296f6X-003` decodes to 0.1340523654409293 = old's macro EXACTLY (difference
  0.0; likewise `s2` fps = .4122428908131096 exactly), while the collapsed
  `+1.128a0bec2e265X-003` = 0.1340523654414653 differs from old by 5.4e-13 (s2: −6.0e-13). Nothing but the collapse separates the two new legs, so the whole
  movement is its summation-order residue.
- **This run:** cont_sg only (the two rows). No other config in this group carries
  `jointtests` with subgroups.
- **Reproduction:** `repro_run.sh`, block `AMESE_R1` (cont_sg; sentinel `REPRO_QRSGJT_DONE`)
  and the new-only block `QRFPS`. Decisive lines:
  ```
  AMESE_R1_joint_horiz             OLD +1.b0445c643364aX-001   NEW(collapsed) +1.b0445c6430cc9X-001   QRFPS +1.b0445c643364aX-001
  AMESE_R1_joint_pre_att_allsgssame OLD +1.74170cb89d771X-03f  NEW(collapsed) +1.74170cb78182fX-03f   QRFPS +1.74170cb89d771X-03f
  AMESE_R1_joint_horiz_s1          OLD [.1340523654409293] macro  NEW .1340523654414653 scalar (+1.128a0bec2e265X-003)   QRFPS +1.128a0bec296f6X-003 = .1340523654409293 (= OLD)
  AMESE_R1_joint_horiz_s2          OLD [.4122428908131096] macro  NEW .412242890813 scalar
  ```
- **Classification: intended-change.**
- **Open questions:** none on substance. **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E4]** — `QR-sgjt-2`'s reason at
  `4280f3b` now names both mechanisms and records that `QR-sgjt-1` was retired *"when the
  cross-kind compare put its 4e-12 move inside rtol"*. The observation as written was: the row
  reason mentions only the collapse, while the
  macro-to-scalar change (CS-60) rides on the same object and is what turns the verbose line
  into REMOVED + ADDITION.

---

### S7. Per-subgroup stacked results: silent label shift via Stata's rowname padding (dupcol)

Ledger rows: `C-SGlab-b`, `C-SGlab-s` (cont_sg_tidose), `C-SGlab-b2`, `C-SGlab-s2`
(cont_sg_gap).

- **Symptom / impact:** when an early per-subgroup "pre" round appended a row but a later one
  failed (here: a dose structure that makes the continuous pre-treatment margins all-omitted),
  old's namelist ended up shorter than the matrix and **Stata silently pads a short namelist by
  repeating the last name**: every later entry wore the next quantity's label. On cont_sg_gap
  `_b[main_SG1_att]` returned 0 (the true SG1 att, 1.3272, sat under `main_SG1_ame`),
  `_b[ES_att_SG0_post1]` returned the post0 value (1.4769; the true post1 estimate, 1.6820, sat
  under `ES_ame_SG0_pre2`), the pretreatment att (−0.8028) was displayed as `main:ame` with a
  missing relyr, and the last four `e(b)` columns all read `ES_ame_SG2_post1`. Values were
  positionally correct; the labels lied.
- **Mechanism (v1.3):** 2815-2829 — in the `pre`/`contpre` rounds the `truerow` match on the
  `over()` names sets `sghaspre no` when the subgroup's row is not found (the `contpre` round's
  `o.`-marked names defeat the match); 2833-2834 append each non-ES row per the RUNNING flag
  (round order main, pre, cont, contpre: the binary pre row is already in); 2867-2874 build the
  namelist from the flag's FINAL value (both `pretreatment:` labels dropped); 2886
  `matrix rownames wdid_sg = rownamessg` with a list one short → padded; 2969-2970 carry the
  same short list into `e(b)`/`e(V)`, so the shift accumulates across subgroups (three missing
  names → three pads → four trailing `ES_ame_SG2_post1`).
- **Fix + proof:** `9d2e4df` (2026-07-26, tranche 12 Wave C; pins `0b12ee1`): a per-(round x
  object) append record `sgrowkeep` (new 3526) keeps the k-th pooled non-ES name iff the k-th
  row exists — names can never exceed rows. Verifier wC_v1: parity strictly matnames-only on 7
  dumps; names invariant on 14 runs across the option surface; 130+ ES entries matched the
  independently-labelled `e(b_es_*_s<sg>)` at reldif 0; replay byte-parity.
- **This run:** cont_sg_gap and cont_sg_tidose: `e(b)` names 1x45 → 1x33 (the count change is
  BS-8, S5; the label repair is this row), `wdidsgresults_0/1/2` 11x7 with a repeated
  `ES_ame:post1` → 7x7 `[main:att pretreatment:att main:ame ES_att:pre2 .. post1]`; **zero value
  cells** differ in these groups.
- **Reproduction:** `repro_run.sh`, block `SGLAB` (cont_sg_gap, as in S5), plus the positional
  alignment read off the two suite dumps with `dumpparse.py` (old `e(b)` slots 13-23 are
  value-identical to new slots 13-19 followed by four zeros):
  ```
  old slot  old label            value               new slot  new label
  [14]      main_SG0_ame         -0.8028388411       [14]      pretreatment_SG0_att
  [15]      ES_att_SG0_pre2       0.0                [15]      main_SG0_ame        (the R6 zero)
  [18]      ES_att_SG0_post1      1.4769030794       [18]      ES_att_SG0_post0
  [19]      ES_ame_SG0_pre2       1.6819851986       [19]      ES_att_SG0_post1
  [23]      main_SG1_att          0.0                (old's own ES_ame_SG0_post1 structural zero)
  [24]      main_SG1_ame          1.3272440931       [20]      main_SG1_att
  [25]      ES_att_SG1_pre2      -1.1721099992       [21]      pretreatment_SG1_att
  old slots [42..45] all labelled ES_ame_SG2_post1
  SGLAB_b_main_SG1_att        OLD 0               NEW 1.3272440931
  SGLAB_b_ES_att_SG0_post1    OLD 1.47690307943   NEW 1.68198519856
  SGLAB_b_pretreatment_SG0_att OLD ABSENT         NEW -.802838841126
  SGLAB_sg0_r2_est / relyr    OLD -.802838841126 / .  (labelled main:ame)   NEW same value, labelled pretreatment:att
  ```
- **Classification: old-bug-fixed.**
- **Open questions:** none.

---

### S8. Silent b-only posting when some SEs are missing: the e(V)-preserve fix on the subgroup and d3 configs

Ledger rows: `A17-sgv-1` (cont_sg_tidose), `A17-sgv-2` (cont_sg_gap); the same machinery
fires as `SIA4-d3p` (trig_d3void, S3) and `A17-uz` (ext_cont_ext_untreatedzero, this group's
config list).

- **Symptom / impact:** whenever any reported estimate's SE was missing — here the voided
  ame rows of the d3 fixtures and the plain-continuous missing `*_ame` SEs of
  ext_cont_ext_untreatedzero — old silently posted **b-only** (`e(properties)` = `b`, no
  caution): lincom/test/esttab downstream lost the VCV for every healthy estimand too.
- **Mechanism (v1.3):** 2976 `matrix ev = diag(hadamard(se, se))` carries a missing diagonal;
  2987 `cap ereturn post eb ev` fails (Stata rejects a VCE with a missing diagonal); 2989
  `cap ereturn post eb` posts the coefficients alone; nothing is printed.
- **Fix + proof:** `47defa0` (2026-07-25) + doc `c53e25b`: a degenerate row gets a negative
  sentinel (−1e100) V diagonal that renders as missing on every consumer surface tested
  (ereturn display, esttab, lincom incl. mixed combinations, nlcom; `test` drops exactly the
  affected constraints); healthy rows keep their true VCV; a caution reports the count. The
  assessor wave refuted the zero-diagonal and identified-rows-only alternatives. Verified 3+3,
  32/32 clean configs bit-identical.
- **This run:** cont_sg_gap / cont_sg_tidose `properties` `[b]` → `[b V]` ("standard errors
  were not identified for 9 of 33 reported estimates"); ext_cont_ext_untreatedzero `[b]` →
  `[b V]` (`A17-uz`; new also adds `e(es_ndegen_ame)`); trig_d3void `[b]` → `[b V]`
  (`SIA4-d3p`, "6 of 12"). On cont_lat both builds post b-only (all-zero VCE), so no A17 row
  fires there — consistent.
- **Reproduction:** `repro_run.sh`, blocks `SGLAB` and `SIA4` (sentinel `REPRO_A17SGV_DONE`).
  Decisive lines:
  ```
  SGLAB_properties  OLD [b]        NEW [b V]
  SGLAB_V_1_1       OLD (no V)     NEW .202693722935     (= .450215196251^2, the SG-run main att se)
  SGLAB_V_3_3       OLD (no V)     NEW -1.000000e+100    (the voided main ame)
  SIA4_properties   OLD [b]        NEW [b V];  SIA4_V_1_1 NEW .002199753379;  SIA4_V_3_3 NEW -1e100
  OLD log: no caution.  NEW log: "Caution: standard errors were not identified for 9 of 33
  reported estimates. e(V) contains placeholder entries for the affected estimates that display
  as missing standard errors; all other estimates retain their variances in e(V)."
  ```
- **Classification: old-bug-fixed** (and, as the 2026-08-23 report's §3.6 records, simultaneously
  a deliberate contract change: partially-degenerate runs post e(V) with sentinel diagonals plus
  a caution instead of posting nothing).
- **Open questions:** none.

---

### Summary

| section | ledger family | kind | old-vs-new magnitude on this group's configs |
|---|---|---|---|
| S1 | AMESE-R1-* | old-bug-fixed | cont_sg: every ame se/t/p/ci → missing (t 2285 → ., contrast t 1.3e7 → .), 18 V diagonals → −1e100, joint_Choriz .0031 removed |
| S2 | AMESE-R4-l/-r, R4W-fitstats | intended-change | cont_lat: rmse 6.6e-4, r2adj 2.2e-9, r2_within 2.3e-6, ll_0 4.3e-2; rss 5e-14; dof 754 → 753; no estimate/SE/test moves |
| S3 | SIA4-* | old-bug-fixed | trig_d3void 6 ame cells 0 → .; sg fixtures: sg ame rows and 6 pwsgs contrasts 0 → . |
| S4 | SI7-* | unclear | e(b) 11 → 12 names (pretreatment_SG1_att added, = 0); the described drop fires only under fullprecisionslow |
| S5 | BS8-cont-sg-* | old-bug-fixed | e(b) 45 → 33 (12 structural zeros dropped); 25 positional cells |
| S6 | QR-sgjt-* | intended-change | 4.0e-12 / 1.4e-12 / 1.8e-10, bit-identical to old under fullprecisionslow |
| S7 | C-SGlab-* | old-bug-fixed | labels only; zero value cells |
| S8 | A17-sgv-* (A17-uz, SIA4-d3p) | old-bug-fixed | properties b → b V with sentinel diagonals + caution |

Relevant-but-untested here: the `A17-uz` config (ext_cont_ext_untreatedzero) and the sibling
configs plot_cont / trig_cfx / si20_lat_cfx_pw were read from the verbose re-judge and the
dumps only, not re-run; nothing in this group required a Stata run beyond the one two-process
reproduction.


---


# Part 4 - the adversarial line-test family


> ## Headline of record and evidence provenance (added 2026-09-13; identical in all six sections)
>
> **Superseding note added 2026-09-14.** Two things below have since moved and are corrected
> here rather than rewritten in place, because the block is accurate about the run it describes.
> (i) `4280f3b` is no longer `HEAD`; the report's own header names the run of record by its
> `prov[new]` sha, and **that** is the build the report's claims attach to. (ii) The closing run
> of this assessment is the **2026-09-14** one over a **559-config** universe against the same
> v1.3 old side, reported in Part 7.5: **559 of 559 compared, zero DIFFERS, unexplained 0,
> stale_ledger 0**. The 2026-09-11 run described below is where this section's EVIDENCE was
> gathered and its magnitudes measured; the 2026-09-14 run is what verifies that nothing still
> differs. Neither supersedes the other - they answer different questions.
>
> **Run of record.** The 2026-09-11 dual-build run: **old = `301b706`** (published v1.3),
> **new = `14c30d3`**, whose `wooldid.ado` is byte-identical to this section's pin `08ec0b4`
> and to current `HEAD` apart from a nine-line comment block and a shortened `dis as text`
> permutation caution (10 insertions, 1 deletion; no estimation code).
>
> **The report judges those dumps against the ledger and comparer at `4280f3b` = `HEAD`.** That
> suite defines a **538-config universe**; the 2026-09-11 dumps cover **507** of them:
> **507 compared — EXPLAINED = 374, MATCH = 25, RC_EXPLAINED = 108 — unexplained 0,
> stale_ledger 0, unscoped 0, prov_fail 0, COVERAGE_OK.** The other 31 configs are new at
> `4280f3b`, have no 2026-09-11 dumps, and are outside every claim in this report.
>
> An earlier headline — **507/507 compared, EXPLAINED = 399, RC_EXPLAINED = 108**, at ledger
> `26ce421` — is quoted by three of the six drafts (`legacy` and `intended` give the full
> 399/108 line; `inference` gives the 507/507). It is the **same dumps and the same 507 configs**
> judged by the pre-`4280f3b` comparer. The two reconcile exactly: `399 = 374 + 25`, the 25 being
> the configs on which MATCH became reachable once the comparer began comparing across batch 3's
> macro→scalar kind change (`4280f3b`, verbatim: *"It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)"*; the same paragraph
> puts the cells affected at *"361 values on 103 configs"*). **Both numbers are measured and
> neither is withdrawn.** See `LEDGER_DELTA.md` §0.
>
> **Evidence provenance.** Every `test/suite/out/run_old|run_new/<config>.csv` cited below is a
> **regenerable** artifact, not a stored one: `test/suite/out/` is gitignored (`test/.gitignore`
> line 5) and the 2026-09-11 dumps were overwritten in place by the 2026-09-13 run. Regenerate
> them with `./run_full.sh oldnew` from `test/suite`, which stages old from
> `git show main:wooldidinstall/wooldid.ado` and new from `git show HEAD:…`. `compare_verbose_0911.txt`
> — cited by four sections, twice as `scratchpad/compare_verbose_0911.txt`, which is **not**
> repo-relative and never was — is likewise a recipe: rebuild it with `compare_oldnew.py … --verbose`,
> passing `--ledger`/`--configs` from whichever commit's judgement you want. The command, the build
> table and the one harness caveat (`4280f3b` changed `mkrun.py`'s `dumpres`, so a regenerated dump
> may differ in its `files|` line or gain a `matnames|…|UNREADABLE` line — **value cells are
> unaffected**) are in **`review_process/assessment_2026-09-13/repro/README.md`**.
>
> **Reproductions.** The do-files this section cites "in this directory" are now parameterised —
> repo root and scratch workdir as positional arguments — and live at
> **`review_process/assessment_2026-09-13/repro/<section>/`**. Their exact invocations are in each
> file's header and in that README. Bulk artifacts (dumps, logs, copied `.ado` builds) are not
> copied: they are regenerable, and the README says how.
>
> **Per-edit log:** `EDITS.md`. **Ledger delta and closed open questions:** `LEDGER_DELTA.md`.


Group: ledger families `A13-esfadv-*`, `A13-altb-*`, `A13-gapb-*`, `NFADV-core-*`, `NFADV-adm-*`,
`NFADV-ext-*` (with the BS-12 clause on every primary-esf p row), plus the rc-by-design pins
`adv_gapa` 198 -> 0, `adv_r2001` 2001 -> 0, `trig_sgadv_deeppre` 2000 -> 0. Builds: old = v1.3
(`main`, 301b706); new = HEAD (08ec0b4); the isolation build for BS-12 is 7e30e85 (= 25727a6^, the
tree immediately before the BS-12 commit). Sources: `compare_verbose_0911.txt` (verbose re-judge of the
2026-09-11 dumps), `test/suite/ledger.csv`, `test/suite/out/run_old|run_new`, the committed battery
baselines at each fix commit (decoded in `baseline_history.txt` here), `wooldid_solved_issues.md`
(B16 main-block entry; the tranche-7 SG-crash entry; the 2026-08-18/19 twelve-commit entry),
`review_process/research_2026-09-08/adversarial_null.md`, `review_process/rulings_2026-09-07.md` J,
`review_process/bug_sweep_2026-09-06.md` BS-12, `test/test_advstat.do`, `test/test_advguard.do`.

This supersedes section 1.8 (and the pointer in 3.5) of the 2026-08-23 report. Figures that changed
since that report: `adv_gapb` 0.00274 -> 1.576e-4, `esf_adv` 8.29e-11 -> 1.84e-11, `adm_esf_adv`
2.98e-12 -> 1.085e-12 (all BS-12); `ext_b1_h11_f1` is no longer a removal but a posted value; and the
"145 ext configs" now decompose as 109 value-changed configs + 36 removal configs. Four layers stack on
the same e() objects (`e(adversarial_pval_*)`, `e(adversarial_line_*)`), so the ledger rows cite several
of them at once; the sections below separate the layers by commit, and the committed baselines plus the
three-build reproduction give the exact value at every step.

Verification state (2026-09-13, resumed from the 2026-09-12 draft on disk): the three build copies are
md5-identical to `git show main|7e30e85|HEAD:wooldidinstall/wooldid.ado`; every old- and new-build line
number below was re-located on those copies; the ext counts (595 cells per `NFADV-ext-*` row = 145 removal
cells on 36 pooled configs + 450 value cells on 109 esf configs; 145 configs on `-l`/`-ln`) were recounted
from `compare_verbose_0911.txt`; the committed baseline chain was re-decoded straight from git
(`baseline_chain_from_git.txt`); and the NEW-build leg of the reproduction was re-executed on 2026-09-13
(`new/launch_time.txt`: 00:41:23-00:41:57, exit 0, five sentinels) and matched the 2026-09-12 run
byte-for-byte on every dump line (`new/prev_run_0912/` holds the earlier run). Every rc, scalar and matrix
cell in `old/out` and `new/out` is bit-identical to the suite's `run_old`/`run_new` dump for the same
config; the only textual difference is the harness's `matnames` line, which omits the `y1` rowname the
suite's dump carries.

**The full chain on the three battery configs (every number is a committed baseline or a dump from this
run; `baseline_history.txt`, `threeway_table.txt`):**

| config | v1.3 (old dump) | c0e17f8 B16 | 4fcf127 NFADV (= 7e30e85, reproduced) | 25727a6 BS-12 (= HEAD, reproduced) |
|---|---|---|---|---|
| `esf_adv` p | 3.53039879927e-10 | unchanged (6267522 svmat: 3.53039968068e-10) | 8.29342036543e-11 | **1.8367718897e-11** |
| `adv_gapb` p | 0.0089204175420555 | **0.0018428792242799** | 0.0027377035367557 | **0.0001576357885506** |
| `adv_altb` p | 6.06579176253e-23 | **1.00214760803e-22** | **3.68856446693e-06** | unchanged |
| `adv_gapa` | rc 198 | rc 0, p 8.51034466062e-07 | 0.0182876111152825 | 4.59070323557e-06 |
| `adv_r2001` | rc 2001 | rc 0, p 9.47227267238e-21 | 1.20615084072e-05 | unchanged |

Bold = the step that moved the value. Reading across a row gives the decomposition the task asks for:
the c0e17f8 column is the bug-fix share (same plug-in statistic, corrected geometry), the 4fcf127 and
25727a6 columns are the two ruled convention changes.

---

#### S1. B16 + 4483758: the esfixedbaseperiod adversarial blocks fitted the wrong points, planted the reference at a data hole, and crashed on degenerate fits (old-bug-fixed)

- **Symptom / impact:** three coupled defects in the four adversarial blocks' `esfixedbaseperiod`
  branch of v1.3, plus a fourth in the two subgroup blocks. (1) *Altmethod cut off-by-one/two:* with the
  fixed reference below the reported window (`esrelativeto(-3)` + `esprelength(2)`), the relyr -1
  coefficient fell out of the pre-trend line fit and was tested as if post-treatment, so
  `altmethod_adversarial` tested the wrong hypothesis (`adv_altb`: p 6.07e-23 on a 1-point-short fit).
  (2) *Reference-gap hunt:* the branch located the omitted reference by scanning for jumps in the
  consecutive-t sequence of the svmat'd coefficient frame. On gap-year data with two jumps the local
  held two numbers and `replace t = 5 7 - 1`-style syntax aborted with r(198) `invalid '7'`
  (`adv_gapa`); with exactly one jump it silently planted the reference zero at the DATA hole instead
  of at the true reference period, so the reported p tested a line through a fabricated point
  (`adv_gapb`: p 0.00892 where the correctly planted reference gives 0.00184 on the same statistic -
  crossing 0.05 on aw/pw/sg/att_it variants per the B16 record). (3) *Degenerate fits:* a 1-row fit
  aborted the whole run with r(2001) `insufficient observations` (`adv_r2001`), and silent
  `finaladvtest` cancellations dropped the requested test with no message. (4) *SG skiplocation
  crash (4483758):* the two subgroup blocks built their altmethod cut from `skiplocation - 1`; with the
  reference below the reported pre-window that local is empty, the restriction degenerated to
  `if t <= -1`, selected zero rows, and `reg` aborted with r(2000) `no observations`
  (`trig_sgadv_deeppre`) - while the identical non-SG run returned rc 0. Who it matters for: anyone
  running `adversarial` with `esfixedbaseperiod` on a panel with calendar gaps (planting: silently
  wrong p), with `esrelativeto()` below the window (cut: wrong hypothesis; SG: crash), or with a short
  pre-window (r(2001) abort).
- **Mechanism (old build, v1.3; MAIN binary block, the other three blocks are copies at the lines in
  parentheses):** `qui svmat `pretrendregmat'` at 1652 (2036/2343/2535) builds the coefficient frame
  (float storage - see S3 for that storage class); the gap hunt is `gen missingfind = t - t[_n-1]`
  (1676), `glevelsof(t) if missingfind == bignjump + 1, ... local(skiplocation)` (1677) and
  `replace t = `skiplocation' - 1 if _n == _N` (1680) - a multi-valued `skiplocation` is the r(198),
  a single wrong jump is the hole planting; the altmethod cut is
  `local altadv_restrict if t <= `espre'+`altmethod_min'-1` (1699; 2078), which under-counts by the
  distance between the reference and the window; the SG blocks instead use
  `if t <= `skiplocation'-1` (2383; 2574), the empty-local r(2000). The fit `qui reg b t
  `altadv_restrict'` (1704; 2083/2387/2579) is uncaptured, so a 1-row or 0-row sample is the
  r(2001)/r(2000) abort. The silent cancellation is `if `finaladvtest' <= 2 { local adversarial ...}`
  at 1640 with no message.
- **Fix + proof:** `4483758` (2026-07-22, SG half: fall back to the MAIN formula only when
  `skiplocation` is empty, r(min) captured into a local first - a symbolic r(min) in the restriction
  resolves to missing at reg time and would silently select ALL rows, proven by the first assessor:
  regN=8 vs the correct 4); `c0e17f8` (2026-07-26: the cut becomes `t <= min_t`, exactly every
  pre-period coefficient plus the planted reference rows in both regimes; the hunt is deleted and the
  reference is planted analytically at `t == refperiod = esrelativeto + min_t + 1`);
  `05a243b` (closed-form guard for the min_t==1 1-row cell, a runtime count backstop before all four
  `reg b t` sites, and one-line Cautions for every cancellation); `44fcf47`/`291090c` verifier
  corrections (Cautions gated on a pending request; backstop threshold 1 + bignjump so the planted
  rows cannot satisfy it alone). Note that the r(2001) cure on `adv_r2001` is c0e17f8's cut fix
  alone - the wrong cut had selected a single row, and at c0e17f8 the config is already rc 0 with
  p 9.47227267238e-21; 05a243b's guards cover the residual min_t == 1 and gap-hollowed cells pinned
  by test_advguard T1/T2b, which no suite config reaches. Decisive evidence in the B16 record: verifier wA_v1 reconstructed
  the corrected `adv_gapb` p from e(b_es/V_es/df_r) alone to reldiff 1.7e-14 AND reproduced the
  parent's wrong p by planting at the hole; four pre-pinned battery configs (5a332ba) moved in
  exactly the pre-registered objects and nothing else moved at rtol 0 (21 -> 27 configs);
  `test/test_advguard.do` T1-T7 ALL_PASS with T1 aborting r(2001) on the parent (clean bisection);
  the SG half: both crash configs go r2000 -> rc 0 with per-subgroup p-values matching the non-SG
  reference (2.19e-13/1.99e-14), 18/18 battery identical, empty-skiplocation-but-no-crash
  counterexample hunt negative. The committed baselines pin the bug-fix share by themselves:
  `adv_gapb` 0.0089204 (5a332ba) -> 0.0018429 (c0e17f8) with the display slope moving 0.5706 ->
  0.3837 because the reference row moved from the hole to the true period; `adv_altb` 6.07e-23 ->
  1.002e-22 with slope -0.0611 -> -0.0913 because the fit gained the relyr -1 coefficient (the commit
  message says "gains the reference row"; the record corrects that to the relyr -1 coefficient).
- **This run:** the B16 share of `A13-gapb-p` (`adv_gapb` p 0.0089204175420555 -> 0.0001576357885506,
  of which 0.00892 -> 0.00184 is B16 and the rest is S3/S4) and `A13-gapb-l` (slope
  0.5706470523561749 -> 0.3837199010520373, entirely B16); the B16 share of `A13-altb-p` (`adv_altb` p
  6.06579176253e-23 -> 3.68856446693e-06, of which 6.07e-23 -> 1.00e-22 is the cut fix and the rest is
  NFADV-02) and `A13-altb-l` (slope -0.06111615896224976 -> -0.09131282072514192, entirely the cut
  fix); rc-by-design pins `adv_gapa` 198 -> 0, `adv_r2001` 2001 -> 0, `trig_sgadv_deeppre` 2000 -> 0
  (their post-fix p-values - 4.59e-06, 1.21e-05, 9.27e-10 with sg0 3.89e-08 / sg1 2.18e-11 - are never
  compared because old has no values). `esf_adv` carries "B16 family" in its `A13-esfadv-p` reason as
  a family label only: its reference sits inside the window on gap-free data, and its baseline did not
  move at c0e17f8 (the only pre-NFADV movement is the 6267522 svmat class, 3.53039880e-10 ->
  3.53039968e-10, 2.5e-7 relative).
- **Reproduction:** `repro_template.do` instantiated as `old/repro_adv_old.do` and `new/repro_adv_new.do`
  (separate Stata processes; sentinel `REPRO_B16_DONE`; the new leg re-run 2026-09-13, identical). The
  batch logs hard-wrap at 80 columns, so the decisive lines are quoted from the %21x dumps decoded in
  `threeway_table.txt` (and `threeway_new_rerun_0913.txt` for the fresh new leg); every rc, scalar and
  matrix cell in `old/out` and `new/out` is bit-identical to the suite's `run_old`/`run_new` dump for the
  same config:
  ```
  adv_gapb @ old :: scalar|adversarial_pval_att|0.0089204175420555      mat|adversarial_line_att|1|1|0.5706470523561749
  adv_gapb @ new :: scalar|adversarial_pval_att|0.0001576357885506      mat|adversarial_line_att|1|1|0.3837199010520373
  adv_altb @ old :: scalar|adversarial_pval_att|6.06579176253e-23       mat|adversarial_line_att|1|1|-0.06111615896224976
  adv_altb @ new :: scalar|adversarial_pval_att|3.68856446693e-06       mat|adversarial_line_att|1|1|-0.09131282072514192
  adv_gapa @ old :: rc|198      (log: invalid '7')          adv_gapa @ new :: rc|0  p 4.59070323557e-06
  adv_r2001 @ old :: rc|2001    (log: insufficient observations)   adv_r2001 @ new :: rc|0  p 1.20615084072e-05
  trig_sgadv_deeppre @ old :: rc|2000  (log: no observations)      new :: rc|0  p 9.26522980657e-10, s0 3.89318837775e-08, s1 2.1784808203e-11
  ```
  The B16-only intermediate (0.0018428792242799 / 1.00214760803e-22, with the line [1,1] already at its
  final 0.3837199010520373 / -0.09131282072514192) is NOT a fourth Stata process: it is the committed
  c0e17f8 battery baseline, decoded straight from git in `baseline_chain_from_git.txt` (`git show
  c0e17f8:test/baseline/adv_gapb.csv`, `.../adv_altb.csv`), bracketed by the 5a332ba pre-fix pins
  (0.0089204192224396 / 6.06579268565e-23, the v1.3 values up to the svmat-float class) and the 4fcf127
  re-bless. The whole chain old -> c0e17f8 -> 4fcf127 -> 25727a6 -> HEAD is therefore closed with a committed
  or reproduced value at every step.

#### S2. NFADV-05 / SI-14: a pooled-reference adversarial p-value was posted although no single omitted reference period exists (old-bug-fixed, by author ruling)

- **Symptom / impact:** without `esfixedbaseperiod` the comparison period is a cohort-specific pool of
  pre-periods, so no single relative period is omitted and there is no (reference, 0) point through
  which a line can be anchored. v1.3 nevertheless posted `e(adversarial_pval_*)` under the pooled
  reference by fabricating one: it appended a zero coefficient one horizon below the earliest real
  horizon, fitted the display line through it, and ran the frozen-line F test. The resulting
  p-value (`cov_advpool`: 2.18e-22; ext pooled cells down to 1e-165) answered no question about the
  estimator's normalisation, and readers took it as a pre-trend diagnostic.
- **Mechanism (old build):** in the pooled branch at 1669-1672, `sum t, meanonly` then
  `replace t = r(min) - 1 if _n == _N` plants the zero row; the same `reg b t` (1704), frozen
  constraints (1720) and `test` (1726) follow; and the posting at 3122/3173 has no esf gate
  (`if "`adversarialP...'" != ""` is the only condition), so every pooled run posted the scalar.
- **Fix + proof:** `4fcf127` (2026-08-19, author-approved: "p-values are now esf-scoped like
  pretrendpower: pooled-reference runs get one Caution and the display line without a p"). New build:
  the worker prints `Caution: the adversarial pre-trend p-value requires esfixedbaseperiod (the pooled
  reference has no single omitted period); the fitted line is reported without a p-value.` (969), and
  `wooldid_advstat` (12247) computes inference only inside `if "`esf'" != ""` (the block at
  12306-12325), returning the line with `status nop`; the posting at 4162/4190 is conditional on the
  local, so nothing is posted. `test/test_advstat.do` pins the pooled Caution / no-p / line-retained
  contract; SI-14 (author-approved 2026-08-25, solved-issues) added `cov_advpool` so a core config
  exercises the path - "3 scoped rows, all fired".
- **This run:** `NFADV-core-p` on `cov_advpool`: `adversarial_pval_att` REMOVED (was 2.1774111081e-22);
  `NFADV-core-l`/`-ln`: the line survives with renamed columns (slope 0.37710243947803984 ->
  0.3771024403538533, a 2.3e-9 svmat-float class; intercept re-expressed, S3). `NFADV-ext-p` removals:
  all **145 removed cells sit on exactly the 36 pooled-reference ext configs** (`_p` suffix; every
  one lacks `esfixedbaseperiod`): `ext_b1_h11_p ext_b1_h22_p ext_b1_h33_p ext_b1_h55_p`,
  `ext_b2_h{10,11,20,22,30,33,55}_p`, `ext_b3_h{10,11,20,22,30,33,55}_p`,
  `ext_b4_h{10,11,20,22,30,33,55}_p`, `ext_c1_h{11,22,33,55}_p`, `ext_c2_h{10,11,20,22,30,33,55}_p`.
  The removed cells are every pooled p the old build posted: `att`, `att_it` (h22/h55 multi-estimand
  rows), per-subgroup `att_s1..s3` / `att_it_s1..s3` (b2/b3/b4 = 1/3/3-level subgroups), and the
  continuous `ame`, `ame_it`, `ame_s1..s3` (c1/c2). Examples: `ext_b1_h22_p` att and att_it 1.13e-40;
  `ext_b3_h22_p` eight cells incl. att_it_s3 3.51e-56; `ext_c2_h33_p` ame_s2 0.0063 and att 2.97e-162.
  Each removal is right for the same single reason: the config has no omitted reference period, so the
  planted point the old test was anchored on was fabricated. **No degenerate-covariance removal fired
  this run** - the ledger reason's second clause ("degenerate-covariance p-values are no longer
  posted") describes a class that BS-12 emptied on this suite (see S4 and open questions); the only
  other unpost surface, the rank-1 `V_es_ame` under cohort clustering (NF-2 in the 08-23 report,
  old posted 0.6147 on a rank-1 covariance), is not reached by any config.
- **Reproduction:** same do-files (sentinel `REPRO_ADVPOOL_DONE`), config `cov_advpool` on
  `test/suite/data/dta_s_covars.dta`:
  ```
  cov_advpool @ old :: scalar|adversarial_pval_att|2.1774111081e-22    matnames|adversarial_line_att|1x2|t _cons
  cov_advpool @ new :: (no adversarial_pval_att line)                   matnames|adversarial_line_att|1x2|slope intercept_rt0
  new log (2026-09-13 run) line 364 [2026-09-12 run: line 376]: Caution: the adversarial pre-trend p-value requires esfixedbaseperiod (the pooled reference has no single omitted period); the fitted line is reported without a p-value.
  ```
  `compare_old_new.txt`: `[cov_advpool] scalar adversarial_pval_att: REMOVED (was 2.1774111081e-22)`;
  the pre-BS-12 build (`mid/out`) already lacks the scalar, so the removal is 4fcf127's, not BS-12's.

#### S3. NFADV-01/02/05: the primary and altmethod statistics recalibrated on the event-study covariance, and the line republished as (slope, intercept_rt0) (old-bug-fixed; the 2026-08-19 ruling governs the remedy)

> **Class revised 2026-09-13** from *intended-change* to *old-bug-fixed*; see `RECONCILIATION.md` §S-3. The rule the report applies is the one `inference` §1 derives from BS-15: **the defect is the wrong number, the remedy is the convention.** A primary statistic rejecting a true null at **0.90 against a nominal 0.05** (T7), ~1.4x nominal on a realistic V, with the altmethod statistic missing the forecast-variance and cross terms, is wrong output, not a convention anyone chose — and the author's ruling verb was "FIX BOTH". The 2026-08-19 ruling is real and governs the REMEDY (recalibration onto the event-study covariance, the line republished as (slope, intercept_rt0)). Classifying the defect by the name of its fix is what produced the contradiction this section's own open question 1 already flagged. **No measured number in this section changes**; the critic verified them. Item 10 of `review_process/held_for_ruling_2026-09-13.md` puts the same question to the author and recommends the same reading.

- **Symptom / impact (what moved and why the author ruled it):** v1.3's primary test froze the fitted
  line's predicted values as constants into `test` against the margins model - an F test of "each
  non-reference coefficient equals its fitted value" that ignores the sampling variability of the
  fitted line, so its null size depends on the margins VCV. The review measured it over-rejecting at
  ~1.4x nominal on a tame realistic esf covariance and at **0.90 against nominal 0.05** in the
  committed heteroskedastic-V null Monte Carlo (`test/test_advstat.do` T7, 200 draws; the
  recalibrated statistic 0.04). The altmethod forecast test compared post coefficients with a
  pre-period projection but omitted the projection's own variance and the pre/post cross-covariance
  (anti-conservative; `adv_altb` 1.0e-22 -> 3.7e-6, matching an independent F(2,79) anchor). And
  `e(adversarial_line_*)`'s second cell was the intercept in the internal nonnegative t coordinate
  (t = relyr + min_t + 1) while the help promised relative time (NFADV-05). The author's ruling
  (2026-08-19, AskUserQuestion, solved-issues): "FIX BOTH adversarial pre-trend statistics ... in a
  new subprogram, esf-scoped with a pooled-regime Caution, folding the NFADV-05 intercept-coordinate
  fix". The ledger and the 08-23 report record this as an author-approved recalibration; the same
  record calls the pre-campaign statistics "review-verified defects" - see open questions.
- **Mechanism (old build):** the frozen-line construction is `local testconstructor ... ([1.`targetflag']
  `horiztocheck'.`_t' = `xbline')` at 1720 with `xbline` read from `predict linearproject, xb` (1707),
  then `est restore __marginsmodel` and `qui test `testconstructor'` (1726): the constraints carry the
  fitted values as numbers, so the test statistic is a plain Wald form of b against fixed constants.
  The altmethod restriction (1699) selects the pre rows for the same `reg b t`, and the same frozen
  `test` is applied to the post rows only - no forecast variance anywhere. The returned line is
  `matrix `adversarialB_...' = e(b)` from `reg b t` (1705-1706): columns `t _cons` in the internal
  coordinate. **New build:** `wooldid_advstat` (12247) rebuilds the display frame (svmat double at
  12290, `reg` for the line), publishes `advL = (slope, _cons + (min_t+1)*slope)` with colnames
  `slope intercept_rt0` (12300-12301; slope bitwise-invariant by OLS affine invariance), then calls
  mata `wooldid_advx` (13705) -> `wooldid_advW` (13606): primary = minimum-distance composite-null
  statistic in the whitened eigenspace of `e(V_es_*)` (generalized inverse at tolerance 1e-12 of the
  largest eigenvalue, q = residual rank, `Ftail(q, df_r)` with a chi-squared fallback); altmethod =
  full forecast covariance `S = VQQ + Aop*VPP*Aop' - Aop*VPQ - VPQ'*Aop'` (13684) with the planted
  reference rows in the pre design; a degenerate S or q <= 0 returns missing and prints
  `Caution: adversarial: the composite-null covariance is degenerate; no p-value computed.` (12322).
- **Fix + proof:** `4fcf127` (2026-08-19). `test/test_advstat.do` ADVSTAT_ALL_PASS: independent
  hand-reconstruction of both statistics from `e(b_es_att)/e(V_es_att)` by a different route (pinv +
  GLS normal equations vs the ado's whitened eigenspace) at ~3-4e-13 relative agreement; slope bitwise
  vs the parent build (pins `+1.5808d87098235X-003` esf_adv, `-1.76046eab851a1X-004` adv_altb);
  intercept identity bitwise; pooled/timetrends Cautions; per-sg reconstruction; T7 null MC. Battery:
  exactly the 5 adversarial configs moved (pval + line colnames + [1,2] only; [1,1] bit-identical),
  24 others IDENTICAL at rtol 0. Blind post-commit verification AGREE (three-route reconstruction
  <= 2.5e-12, per the 08-23 report). The committed baselines isolate the NFADV share: `esf_adv`
  3.53039968e-10 (6267522) -> 8.29342036543e-11 (4fcf127); `adv_gapb` 0.0018429 (c0e17f8) ->
  0.0027377 (4fcf127); `adv_altb` 1.00214760803e-22 (c0e17f8) -> 3.68856446693e-06 (4fcf127), and
  the pre-BS-12 build 7e30e85 reproduces all three of those values bit-for-bit in this reproduction.
  Finite-sample calibration of the recalibrated statistic on realistic designs is the ADVBOUND campaign
  (`test/suite/mc/REPORT_ADVBOUND_2026-08-26.md`, harness committed at e478711 on 2026-08-25, i.e. on the
  4fcf127 two-parameter statistic, before BS-12): 25 configs, 17,560 reps, zero missing p-values; binary
  `adversarial` KEEP-GATED (flagship certified, sizes 4.80-7.20% at the boundary cells; gate m* = 6
  clusters per cell, n* >= 180 clusters); `altmethod_adversarial` KEEP and strictly dominating the
  primary at every boundary; continuous `adversarial` (ame) formally WITHDRAW, recommended GATE ON CELL
  SIZE (31.55% at 6 obs/cell -> 9.92% at 20 -> 5.95% PASS at 50 under robust). Two consequences for this
  report: ADVBOUND measured the pre-BS-12 statistic, so the one-parameter form's finite-sample size on
  that grid is unmeasured (the research record's level-validity statements are for known Gaussian V); and
  the ame p-values on the `ext_c*` value cells sit in the small-cell class ADVBOUND ruled oversized on the
  NEW statistic - a validity caveat on both builds, not part of the old-vs-new movement (open question 9).
- **This run:** the NFADV-01 share of `A13-esfadv-p` (`esf_adv` 3.53e-10 -> 8.29e-11 of the total
  3.53039879927e-10 -> 1.8367718897e-11), `A13-gapb-p` (0.00184 -> 0.00274 of the total) and
  `NFADV-adm-p` (`adm_esf_adv` 2.51529006878e-11 -> 1.08542355053e-12 in total; the NFADV-01 vs BS-12
  split is NOT isolated on the current tree - the 08-23 run's intermediate 2.98e-12 was measured at
  e53671f, before 111ecbb and the SI-20 batch, and no battery baseline exists for this config); the whole of `A13-altb-p` beyond the cut
  fix (1.00e-22 -> 3.69e-06, NFADV-02; BS-12 does not touch altmethod) and of `NFADV-ext-p` on
  `ext_bin_ext_altadv` (1.5698170341e-137 -> 1.7202424073e-110). NFADV-05 is the whole of every
  `-l` and `-ln` row: `A13-esfadv-l` [1,2] -1.024264054638999 -> 0.4876065376068326 (reldiff 1.48),
  `A13-altb-l` [1,2] 0.3666969537734986 -> -0.2638729079523001, `A13-gapb-l` [1,2] -2.3153004958516066
  -> 0.6363311527163849 (this cell also carries the B16 planting move), `NFADV-adm-l` [1,2]
  -1.0295867361128335 -> 0.34718056941453757, `NFADV-core-l` [1,2] -0.2550473727285864 ->
  0.8762599506407355, `NFADV-ext-l` on 595 cells; the `-ln` rows (`t _cons` -> `slope intercept_rt0`)
  fire on 1 + 1 + 1 + 1 + 1 + 595 matrices. The [1,1] slope cells are bit-identical wherever B16 did
  not move the fit and the old frame was not float-truncated; where they differ at 1e-9 (`esf_adv`
  0.16798562157366959 -> 0.16798562139139653; `adm_esf_adv` 2.0e-9; `cov_advpool` 2.3e-9;
  `ext_b1_h11_f1` 5.3e-8) that is the `6267522` svmat-float storage class the ledger folds in
  ("+ svmat"), not NFADV. The 450 value-changed `NFADV-ext-p` cells on the 109 esf ext configs
  (`_f1/_f2/_f4` + `ext_bin_ext_altadv`) carry NFADV-01 (or -02) and BS-12 together; typical
  magnitudes are whole orders in the tail (`ext_b2_h20_f1` att 5.68e-69 -> 7.67e-107), moderate on
  the continuous ame cells (`ext_c1_h11_f1` ame 0.9096 -> 0.1928; `ext_c2_h22_f2` ame_s3 0.8119 ->
  4.09e-05).
- **Reproduction:** same do-files plus `mid/repro_adv_mid.do` on the 7e30e85 build (sentinel
  `REPRO_NFADV_DONE`; 7e30e85's `wooldid.ado` md5-matches `git show 7e30e85:...`). Old -> pre-BS-12
  isolates NFADV-01/02/05 (with B16 and svmat riding along on gapb/altb):
  ```
  esf_adv  @ old :: p 3.53039879927e-10   line [t _cons]            [1,2] -1.024264054638999
  esf_adv  @ mid :: p 8.29342036543e-11   line [slope intercept_rt0] [1,2] 0.4876065376068326   (slope 0.16798562139139653 both mid and new)
  adv_altb @ old :: p 6.06579176253e-23   adv_altb @ mid :: p 3.68856446693e-06   adv_altb @ new :: p 3.68856446693e-06 (BS-12 inert on altmethod)
  adv_gapb @ mid :: p 0.0027377035367557  (= the 4fcf127 baseline; c0e17f8's 0.0018428792242799 is the B16-only value)
  ```

#### S4. BS-12: the primary adversarial null is the one-parameter line through the reference period (intended-change, author ruling 2026-09-08)

- **Symptom / impact (what moved and why the author ruled it):** under `esfixedbaseperiod` every
  reported coefficient is `b_h = theta_h - theta_r`, so a linear underlying path implies `b_h =
  c*(h - r)`: a line THROUGH the omitted reference period, one free parameter. The NFADV-01
  statistic (and, in its own way, v1.3's) minimised over the free-intercept pair `[1, relyr]`; because
  the constant vector lies in that column space the statistic was exactly invariant to adding the same
  number to every coefficient - which is precisely what an anomalous or anticipatory reference period
  produces - so its power against a reference-period level shift of ANY size was exactly the nominal
  level (research record: p .75 analytic / .83 boottstat / .90 permute on a panel with a five-unit
  anomaly at each cohort's own reference period, every event-study coefficient near -5 with |t| 41-68
  and joint pre-trend p 0.0000). The research record is explicit that this is NOT a size failure - the
  old test was exact-size for the null it labelled (measured .0508 against a pure level shift, 40,000
  draws) - but a different null; and because col(rel - esrel) is nested in col([1, rel]) the
  one-parameter statistic's noncentrality is at least as large against every alternative (a theorem
  of nested minimisation; 0 violations in 20,000 random smooth alternatives), at the price of one
  degree of freedom (q rises from rank(V) - 2 to rank(V) - 1). The ruling (rulings_2026-09-07.md J,
  2026-09-08): "Adopt the one-parameter null - the line through the reference period, matching the
  estimator's normalisation and what `wooldid_ptp` already does." Consequences the ruling accepted:
  every primary esf p moves again on all three live surfaces (analytic, `boottstat`, `permute()`);
  two-real-horizon cells that NFADV-01 refused (q = 0) become testable at q = 1; under `timetrends`
  two periods are normalised so the only line through both is the zero line and the null collapses to
  b = 0 (a strictly stronger test, documented in the help).
- **Mechanism (pre-BS-12 build 7e30e85 vs HEAD):** in mata `wooldid_advW`, `alt == 0` branch, the
  design was `X = J(rows(rel), 1, 1), rel` with `q = k1 - rankM` (= rank(V) - 2 when both columns are
  identified); HEAD has `X = rel :- esrel` (13649) with `q = k1 - rankM` (13653; = rank(V) - 1) and a
  `tt` branch `q = k1; W = z'z` (13644-13646). The printed label changed from "all (non-reference) att
  coefficients are on a single line" to "the att coefficients and the reference period are on a
  single line" (10937/11081). Altmethod (`alt == 1`, 13660-13703), the display line
  (`e(adversarial_line_*)`, deliberately still a free-intercept two-parameter fit) and the pooled unpost
  are untouched - which is exactly what the reproduction shows.
- **Fix + proof:** `25727a6` (2026-09-10) plus `60256b5` (re-pin test_advguard T2b to a posted p).
  `test/test_advstat.do` T1/T6 reconstruct the one-parameter statistic independently (pinv + normal
  equations against the ado's whitened eigenspace; `X = rel :- esrel`, `q = rV - 1`), re-pinned to
  q = 5 and p = 1.83677e-11 on `esf_adv` - the value the battery and this run carry; T7 re-pinned to
  q = 2 (4% rejection vs nominal 5%; the old frozen-line statistic 90% on the same draws); T2/T3
  deliberately not re-pinned (altmethod and display-line tripwires) and unchanged. Commit-verified:
  the battery moved on exactly `esf_adv`, `adv_gapa`, `adv_gapb` and only on `adversarial_pval_att`
  (line matrices byte-identical; `adv_altb` and `adv_r2001` unmoved). Research record
  `research_2026-09-08/adversarial_null.md`: the closed form `b'V^-1 b - (m'V^-1 b)^2/(m'V^-1 m)`
  with `m = rel - esrel` reproduces the new statistic to 4.2e-14 over 200 random full-rank designs and
  IS the eps -> 0 limit of planting the reference at variance eps, which the two-parameter statistic
  was not (W 80 -> 400 as eps -> 0 vs the shipped 3e-29).
- **This run:** the BS-12 share of every primary-esf p row: `A13-esfadv-p` 8.29342036543e-11 ->
  1.8367718897e-11 (of the total 3.53e-10 -> 1.84e-11), `A13-gapb-p` 0.0027377035367557 ->
  0.0001576357885506 (of 0.00892 -> 1.58e-4), `NFADV-adm-p` (total 2.51529006878e-11 -> 1.08542355053e-12;
  the split is not isolated on the current tree, see S3), `NFADV-ext-p` on the 108 primary esf ext configs (all `_f*` rows; every
  primary esf ext value cell carries both NFADV-01 and BS-12 - `ext_b1_h55_f4` att 1.17294108801e-85 ->
  1.98694557158e-86, `ext_b3_h22_f1` att 4.37569631777e-36 -> 8.3818954611e-111 under timetrends where
  the null is now b = 0), and the SI-20 rows `SI20-I-si20_ic_esf-advp` (documented with the SI-20
  group). Two cells deserve naming. (i) `ext_b1_h11_f1` (`b_es_att` is 1x2: `pre2 post0`): v1.3
  posted 0.0002864706283682, NFADV-01 refused it (q = 2 - 2 = 0; the 08-23 report lists it as a
  removal, and the 7e30e85 log here prints the degenerate-covariance Caution), BS-12 posts
  1.64433515556e-28 at q = 1 - a removal in the last report that is a value this run. Its value equals
  `ext_b1_h11_f2`'s (esrelativeto(-2); v1.3 gave 2.52524655326e-22 there): the one-parameter null
  "periods {-2, -1, 0} are collinear" does not depend on which of them is the reference, so the
  Wald statistic is reparameterisation-invariant across the two configs - a sanity property of the
  new null that v1.3's plug-in test lacked (2.86e-4 vs 2.52e-22 on the same three periods); the
  free-intercept form posted nothing on either cell, so it offers no comparison. (ii) The
  timetrends + esf ext cells (`ext_b3_*_f*`, `ext_c2_*_f*`) test b = 0 on the new build; v1.3 posted
  a frozen-line p there (`ext_b3_h22_f1` att 4.38e-36, s1-s3 4.5e-17..1.6e-21). Not moved by BS-12:
  `A13-altb-p`, `adv_r2001`, `trig_sgadv_deeppre` (altmethod), `NFADV-core-*` and the 36 pooled ext
  configs (no p either way), every `-l`/`-ln` row.
- **Reproduction:** the three-build run (sentinel `REPRO_BS12_DONE`); `compare_mid_new.txt` (7e30e85 vs
  HEAD at rtol 0) shows the adversarial objects differing ONLY as `adversarial_pval_att` on `esf_adv`,
  `adv_gapb`, `adv_gapa` and ADDED on `ext_b1_h11_f1`; every other difference in that compare is the
  PS-62 `e(V)` block splice (`wdid_vcvstate analytic-diagonal -> analytic-blocks`), which does not
  touch `e(V_es_*)` or any adversarial object:
  ```
  esf_adv       @ mid :: p 8.29342036543e-11    @ new :: p 1.8367718897e-11     line [1,1]/[1,2] identical mid vs new
  adv_gapb      @ mid :: p 0.0027377035367557   @ new :: p 0.0001576357885506   line identical
  adv_gapa      @ mid :: p 0.0182876111152825   @ new :: p 4.59070323557e-06    line identical
  ext_b1_h11_f1 @ mid :: (no p; log 454: Caution: adversarial: the composite-null covariance is degenerate; no p-value computed.)
  ext_b1_h11_f1 @ new :: p 1.64433515556e-28
  adv_altb / adv_r2001 / trig_sgadv_deeppre / cov_advpool: adversarial objects bit-identical mid vs new
  ```

---

#### Open questions and record discrepancies

1. **Classification of NFADV-01/02 (S3).** The ledger rows, the 08-23 report (1.8 and 3.5) and this
   section call the recalibration an author-approved intended change, and the movement on
   well-behaved cells is a replaced statistic rather than a coding slip. But the solved-issues record
   introduces the same items as "review-verified defects" with measured miscalibration (primary 0.90
   null rejection at nominal 0.05 in T7; ~1.4x nominal on a realistic V; altmethod missing the
   forecast-variance and cross terms), and the author's ruling verb was "FIX BOTH". A reader could
   defensibly file NFADV-01/02 under old-bug-fixed; the kind above follows the ledger's wording and the
   task's framing. NFADV-05 (line coordinates) and BS-12 are unambiguously convention changes; the
   pooled unpost (S2) is filed as a bug fix on the strength of the "no single omitted reference period"
   argument, although the ledger's own wording for it is also "author-approved".
2. **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E9]** — `4280f3b` rewrote it; the clause is now explicitly marked
   historical (*"the pre-BS-12 degenerate-covariance withholding clause is historical"*). The
   observation as written was:
   **`NFADV-ext-p` reason text is partly stale.** It still says "pooled-reference and
   degenerate-covariance p-values are no longer posted (removals)"; this run's 145 removals are all
   pooled-reference, and the degenerate-q class it had in mind (two real horizons, q = 0) was emptied
   by BS-12 (`ext_b1_h11_f1` is now a value). The row still explains every firing; only the clause is
   historical.
3. **`BS12-gapa-p`.** 25727a6's message says a new ledger row for `adv_gapa` was added; it is absent
   from HEAD's ledger. `f3f0c83` removed it ("added during batch 4, never fired: adv_gapa is
   rc-explained") - correct, since old has no `adv_gapa` values to compare. No action.
4. **The B16 record's residual "adversarial is a silent no-op under timetrends in every regime"**
   describes a campaign-era HEAD, not v1.3: the old dumps and this run's old build post p-values on
   the timetrends ext cells (`ext_b3_h22_f1` 4.38e-36). Do not import that residual as an old-build
   symptom.
5. **`esf_adv`'s "B16 family" citation** in `A13-esfadv-p` is a family label; that config's p has no B16
   share (baseline unmoved at c0e17f8). The row is still correct as an explanation of the firing.
6. **Three normalisations now coexist in one output block** (research record, open question 2-3): the
   primary p tests the line through the reference (one parameter), the altmethod p treats the planted
   reference as an ordinary observation in the forecast design, and `e(adversarial_line_*)` remains a
   free-intercept two-parameter display fit. Flagged in the record, not a defect finding.
7. **Untested here:** the `boottstat` and `permute()` surfaces of the adversarial p (the research
   record and 25727a6 measured them on the BS-12 design: .83/.90 -> the .0164 replication floor); the
   NF-2 rank-1 `V_es_ame` refusal (no suite config reaches it); `adm_esf_adv` and the ext grid were
   not re-run (their values are taken from the 2026-09-11 dumps, verified here to be the dumps the
   verbose re-judge read, and the 08-23 report). Only the NEW-build leg of the eight-config reproduction
   was re-executed on 2026-09-13; the old and 7e30e85 legs are the 2026-09-12 runs on disk, whose builds
   are md5-verified and whose adversarial values match the committed baselines and the suite dumps.
9. **A both-builds validity caveat on the continuous (ame) ext cells, not a movement.** ADVBOUND
   (2026-08-26) found the continuous `adversarial` test oversized at small cells on the recalibrated
   statistic (31.55% at 6 obs/cell under robust; 5.95% only at 50), and the v1.3 frozen-line statistic
   is worse-calibrated still (T7: 0.90 null rejection). The `ext_c1_*_f*` / `ext_c2_*_f*` ame value cells
   in `NFADV-ext-p` (e.g. `ext_c1_h11_f1` ame 0.9096 -> 0.1928; `ext_c2_h22_f2` ame_s3 0.8119 -> 4.09e-05)
   therefore move from one poorly calibrated p to a better-but-gated one; the ledger row explains the
   movement, but neither side's ame p on those small synthetic cells should be read as a calibrated
   diagnostic. No gate is implemented in HEAD for those cells (ADVBOUND's recommendation, not a ruling
   in this group's records) - flagged, not a finding this section can close.
10. **BS-12's finite-sample size is unmeasured on the ADVBOUND grid.** The research record's
   level-validity results are known-V Gaussian (40,000 draws) and T7 is a 200-draw known-covariance MC;
   ADVBOUND's 17,560-rep realistic-design calibration was run on the pre-BS-12 two-parameter statistic.
   The record's own open question 1 (whether the added constant direction interacts with small-cluster
   oversizing) is still open.
8. **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E10]** — `4280f3b` corrected the `configs.csv` note. (This section's
   account is also the one `intended` I.8b was corrected to on 2026-09-13; see `EDITS.md`, B-3.)
   The observation as written was:
   **`configs.csv`'s note for `adv_r2001` is stale.** It reads "B16: old degenerate 1-row fit rc2001;
   new clean skip + Caution", but the new build posts p 1.20615084072e-05 and a line for that config
   with no skip and no Caution (this run's `new/repro_adv_new.log` and dump): the cut fix made the
   cell estimable, and the skip-with-Caution path is the min_t == 1 / gap-hollowed class of
   test_advguard T1/T2b. The rc pin (2001 -> 0) is correct; only the annotation is wrong.

#### Artifacts

`repro_template.do` and the three instantiated `repro_adv_*.do` (plus `hexdec.py`) are at
`review_process/assessment_2026-09-13/repro/adversarial/`, parameterised on repo root and scratch
workdir. The rest below are drafting-session artifacts, regenerable by re-running them:

`repro_template.do`; `old/repro_adv_old.do` + `.log`, `mid/repro_adv_mid.do` + `.log` (build 7e30e85),
`new/repro_adv_new.do` + `.log` (the 2026-09-13 re-run; `new/prev_run_0912/` holds the 2026-09-12 log
and dumps, byte-identical; `new/launch_time.txt`); dumps `old/out`, `mid/out`, `new/out` (%21x);
`threeway_table.txt` and `threeway_new_rerun_0913.txt` (decoded); `compare_old_new.txt`,
`compare_mid_new.txt` (test/compare.py, rtol 0); `baseline_history.txt` and `baseline_chain_from_git.txt`
(decoded committed baselines, the latter re-derived on 2026-09-13); `nfadv_ext_p.txt` (the 595
`NFADV-ext-p` firings with config names); `hexdec.py`; `NOTES.md`; `SECTION_draft_0912.md.bak` (the
draft this supersedes).


---


# Part 5 - the inference-convention families of the 2026-09 batches


> ## Headline of record and evidence provenance (added 2026-09-13; identical in all six sections)
>
> **Superseding note added 2026-09-14.** Two things below have since moved and are corrected
> here rather than rewritten in place, because the block is accurate about the run it describes.
> (i) `4280f3b` is no longer `HEAD`; the report's own header names the run of record by its
> `prov[new]` sha, and **that** is the build the report's claims attach to. (ii) The closing run
> of this assessment is the **2026-09-14** one over a **559-config** universe against the same
> v1.3 old side, reported in Part 7.5: **559 of 559 compared, zero DIFFERS, unexplained 0,
> stale_ledger 0**. The 2026-09-11 run described below is where this section's EVIDENCE was
> gathered and its magnitudes measured; the 2026-09-14 run is what verifies that nothing still
> differs. Neither supersedes the other - they answer different questions.
>
> **Run of record.** The 2026-09-11 dual-build run: **old = `301b706`** (published v1.3),
> **new = `14c30d3`**, whose `wooldid.ado` is byte-identical to this section's pin `08ec0b4`
> and to current `HEAD` apart from a nine-line comment block and a shortened `dis as text`
> permutation caution (10 insertions, 1 deletion; no estimation code).
>
> **The report judges those dumps against the ledger and comparer at `4280f3b` = `HEAD`.** That
> suite defines a **538-config universe**; the 2026-09-11 dumps cover **507** of them:
> **507 compared — EXPLAINED = 374, MATCH = 25, RC_EXPLAINED = 108 — unexplained 0,
> stale_ledger 0, unscoped 0, prov_fail 0, COVERAGE_OK.** The other 31 configs are new at
> `4280f3b`, have no 2026-09-11 dumps, and are outside every claim in this report.
>
> An earlier headline — **507/507 compared, EXPLAINED = 399, RC_EXPLAINED = 108**, at ledger
> `26ce421` — is quoted by three of the six drafts (`legacy` and `intended` give the full
> 399/108 line; `inference` gives the 507/507). It is the **same dumps and the same 507 configs**
> judged by the pre-`4280f3b` comparer. The two reconcile exactly: `399 = 374 + 25`, the 25 being
> the configs on which MATCH became reachable once the comparer began comparing across batch 3's
> macro→scalar kind change (`4280f3b`, verbatim: *"It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)"*; the same paragraph
> puts the cells affected at *"361 values on 103 configs"*). **Both numbers are measured and
> neither is withdrawn.** See `LEDGER_DELTA.md` §0.
>
> **Evidence provenance.** Every `test/suite/out/run_old|run_new/<config>.csv` cited below is a
> **regenerable** artifact, not a stored one: `test/suite/out/` is gitignored (`test/.gitignore`
> line 5) and the 2026-09-11 dumps were overwritten in place by the 2026-09-13 run. Regenerate
> them with `./run_full.sh oldnew` from `test/suite`, which stages old from
> `git show main:wooldidinstall/wooldid.ado` and new from `git show HEAD:…`. `compare_verbose_0911.txt`
> — cited by four sections, twice as `scratchpad/compare_verbose_0911.txt`, which is **not**
> repo-relative and never was — is likewise a recipe: rebuild it with `compare_oldnew.py … --verbose`,
> passing `--ledger`/`--configs` from whichever commit's judgement you want. The command, the build
> table and the one harness caveat (`4280f3b` changed `mkrun.py`'s `dumpres`, so a regenerated dump
> may differ in its `files|` line or gain a `matnames|…|UNREADABLE` line — **value cells are
> unaffected**) are in **`review_process/assessment_2026-09-13/repro/README.md`**.
>
> **Reproductions.** The do-files this section cites "in this directory" are now parameterised —
> repo root and scratch workdir as positional arguments — and live at
> **`review_process/assessment_2026-09-13/repro/<section>/`**. Their exact invocations are in each
> file's header and in that README. Bulk artifacts (dumps, logs, copied `.ado` builds) are not
> copied: they are regenerable, and the README says how.
>
> **Per-edit log:** `EDITS.md`. **Ledger delta and closed open questions:** `LEDGER_DELTA.md`.


Group: `BS15-pw / BS15-jt / BS15-jtp`, `BS10-dfr`, `CS60-jhm` (the four `CS60-jhm-*` rows were
retired at `4280f3b`; §3 is now framed on the cells, not the rows), `TM2-sg / TM2-ext`,
`DDG-pw / DDG-jt`, `A18`, `A16-imd / A17 / A17-part`, `A3 / B4-ext`.

Builds: **old** = published wooldid v1.3 (`git main`, `301b706`); **new** = git HEAD
`08ec0b4`. The 2026-09-11 dual-build dumps carry `prov new|14c30d3`; the only ado change
between `14c30d3` and `08ec0b4` is one permutation-caution string (commit `08ec0b4`, plus a
comment), so those dumps stand for HEAD on every family below. Sources: per-config re-judge
`compare_verbose_0911.txt` (507/507 compared, `unexplained=0 stale_ledger=0` — **at ledger
`26ce421`**; see the headline block above), ledger
`test/suite/ledger.csv` (row ids cited), dumps `test/suite/out/run_old|run_new/<config>.csv`
(%21x), records `wooldid_solved_issues.md` and `review_process/*.md`. Old-code line numbers
refer to v1.3 as extracted into `old/wooldid.ado` in this directory (3,450 lines) and were
located there, not taken from any record; new-code lines refer to `new/wooldid.ado` (HEAD).

Reproduction: `repro_inference.do` (now at `review_process/assessment_2026-09-13/repro/inference/`),
run once per build in SEPARATE Stata
processes from the build's own subdirectory (`old/`, `new/`), each printing
`REPRO_<id>_DONE` per block and dumping every e() cell in %21x to `<config>.csv`. Logs:
`old/repro_inference.log`, `new/repro_inference.log`. Decisive lines are pasted under each
section; values are quoted from those logs.

---

### 1. BS-15: subgroup pairwise standard errors and both across-subgroup equality tests computed from a rank-deficient contrast covariance, uncautioned (`BS15-pw`, `BS15-jt`, `BS15-jtp`)

- **Symptom / impact:** with `subgroups()` and the clustering at or above the cohort
  variable `i` — the `cluster(i) subgroups()` call the help recommends — the saturated
  cell dummies leave each cell's residuals with no within-cluster variation, so the two
  subgroups' level standard errors agree to several digits and the analytic pairwise
  contrast standard error is the cancellation residue between them, not an estimated
  variance. v1.3 printed and posted that residue as inference: every `e(pwsgs)` row's se,
  t, p and 95% CI (columns 2-6) plus `e(joint_att_allsgssame)` and
  `e(joint_pre_att_allsgssame)`, with no caution. On the sweep fixture
  (`dta_binary`, `att cluster(i) subgroups(sg)`) the contrast se was `.0022844` — 1/28 of
  either level se and **1/44** of the same contrast's robust se (`.1002352`), implied
  correlation `.99933` — printed with `t 30.8`, `p 9.2e-46` and a CI excluding zero for a
  difference whose honest interval is `[-.127, .270]`; under `robust` the same equality
  tests give `.4831`/`.4393` (bug_sweep_2026-09-06.md, BS-15, "Reproduction"). Six wrong
  numbers per contrast row plus two joint p-values on the default analytic path: a user
  comparing subgroup effects published false-certainty differences. The bootstrap,
  `boottstat` and `permtstat` paths already blanked or cautioned there; the analytic site
  was the only silent one of the four.
- **Mechanism (v1.3):** `old/wooldid.ado` L2456-2466: `if wordcount("`sglist'") > 1 {` →
  per estimand `margins ... over(`subgroups') contrast(overjoint)` then
  `local joint_`item'_sdiff = r(p)[1,2]` (L2459-2460), then
  `margins ... over(`subgroups') pwcompare(effects)` and
  `matrix `pwsgs_`item'' = r(table_vs)[1..6, ...]'` (L2464-2466) — `r(table_vs)` sliced
  verbatim into `e(pwsgs)`; pre-treatment twin L2472-2481; continuous-treatment twins
  L2650-2670; `ereturn matrix pwsgs` at **L3188** (L3191/L3193 are the
  `joint_*_allsgssame` / `joint_pre_*_allsgssame` scalars). No screen on the covariance's rank anywhere, no
  caution; the SG display block prints the rows. (The E1 chain-test rank guard added in
  2026-07 fires only when `test` *drops* constraints; the residue keeps numerical full
  rank, so it passed — bug_sweep BS-15 "Verifier".)
- **Fix + proof:** author ruling 2026-09-07 (`rulings_2026-09-07.md` §A, BS-15): **"Blank
  with a caution"** — when the subgroup pairwise covariance is numerically rank-deficient,
  blank the pairwise SEs and both equality tests with a named caution, matching the
  bootstrap/boottstat/permtstat paths. Built as `cfac3f5` (2026-09-10): new subprogram
  `wooldid_sgdeg` (`new/wooldid.ado` L10614), called at L4319 on the analytic path only and
  only with >1 subgroup. It rebuilds each pairwise row's implied correlation by polarisation,
  `rho = (V_a + V_b - V_ab) / (2*sqrt(V_a*V_b))`, from the level-se columns of
  `e(wdidsgresults_<sg>)` and the contrast se in `e(pwsgs)`; a row with `1 - rho < 1e-2`
  has columns 2-6 set to missing and the equality scalar covering it left unposted, with one
  Caution naming the count and the largest rho (`e(wdid_sgdegen)`, `e(wdid_sgdegenr)`,
  `e(wdid_sgdegenact)` = `blank`); `safety(off)` reports instead (`act` = `report`). Inert
  under `bootstrap()`/`permute()` (which read the analytic se as their own input),
  `fastnose` and `suppresssgprimaryeffects`. The cut `1 - 1e-2` was put to the author as a
  judgment call the ruling had not fixed (`held_for_ruling_2026-09-09.md`, "BS-15 cut":
  rho `.99940`/`.99715` on the two dta_binary rows, `1 - 1e-3` caught only the first) and
  ruled **"Keep 1 - 1e-2"** on 2026-09-10. Tests `test/test_sgdegen.do` G1-G6 (G2: a
  healthy robust design is bit-identical, `.48310953`, `.43928221`, contrast se
  `.1002352`). Commit record: 20-leg fingerprint identical, battery moved on seven
  subgroup configurations at rho `.9990`-`1.0000`, withheld equality p's `1.1e-12`,
  `3.3e-88`, `9.2e-46`, `4.2e-25`; baselines re-blessed.
- **This run:** rows `BS15-pw-adm_svpsg` / `BS15-jt-adm_svpsg` on `adm_svpsg` (realistic
  admin panel, `cluster(state) subgroups(region)`): `e(pwsgs)` 25 cells — the five
  main-att rows `att:2vs1`, `3vs1`, `4vs1`, `4vs2`, `4vs3` lose columns 2-6 (old e.g.
  `att:2vs1` est `.0657`, se `.00629`, t `10.45`, p `7.6e-14`); `att:3vs2` is retained;
  the six pre rows are retained; `e(wdid_sgdegen)` = 5, largest rho `.9971`;
  `joint_att_allsgssame` `3.97e-20` → `.`; `joint_pre_att_allsgssame` `7.7e-68`
  unchanged. `BS15-jt/jtp-cont_sg_gap` (`1.11e-12` → `.`, `2.76e-12` → `.`;
  `e(wdid_sgdegen)` = 6, rho `.9990`), `BS15-jt/jtp-cont_sg_tidose` (`3.25e-88` → `.`,
  `5.03e-10` → `.`), `BS15-pw/jt/jtp-esptrip` and `BS15-pw/jt/jtp-plot_sg` (both on
  dta_binary: `9.16e-46` → `.`, `4.17e-25` → `.`; `pwsgs` 10 cells, old se `.0022844`,
  t `30.80` — the sweep's own fixture numbers).
- **Reproduction (`adm_svpsg`, exact suite options incl. `makeplots saveplots()`):**

```
OLD  old/repro_inference.log  (no caution of any kind in the block; both joint p's printed as real)
     console:  p-value for joint test that all SGs have the same att: 0.0000
               p-value for joint test that all SGs have the same pre-treatment att: 0.0000
     e(pwsgs) att rows     estimate        se          t           p        lb_95ci   ub_95ci
       2vs1.region          .06572      .006288     10.45     7.559e-14     .05307    .07837
       3vs1.region         -.06755      .00823     -8.208     1.248e-10    -.08411   -.051
       4vs1.region         -.003665     .004961    -.7387     .4637        -.01365    .006316
       3vs2.region         -.1333       .009792   -13.61      6.378e-18    -.153     -.1136
       4vs2.region         -.06939      .006594   -10.52      6.044e-14    -.08265   -.05612
       4vs3.region          .06389      .004473    14.28      1.019e-18     .05489    .07289
     BS15_joint_att_allsgssame=3.973e-20      BS15_joint_pre_att_allsgssame=7.717e-68
     BS15_sgdegen=. act=[] r=[]
     implied correlation per row, rho = (V_a + V_b - V_ab) / (2*sqrt(V_a*V_b)), from the level se in
     e(wdidsgresults_<sg>)[main:att, se] and the contrast se above (log lines re-joined at the wrap):
     BS15_row att:2vs1  se_contrast= 0.006288  se_a= 0.055532  se_b= 0.055878  rho= 0.99365  p= 7.56e-14
     BS15_row att:3vs1  se_contrast= 0.008230  se_a= 0.057745  se_b= 0.055878  rho= 0.99004  p= 1.25e-10
     BS15_row att:4vs1  se_contrast= 0.004961  se_a= 0.056689  se_b= 0.055878  rho= 0.99622  p= .464
     BS15_row att:3vs2  se_contrast= 0.009792  se_a= 0.057745  se_b= 0.055532  rho= 0.98581  p= 6.38e-18
     BS15_row att:4vs2  se_contrast= 0.006594  se_a= 0.056689  se_b= 0.055532  rho= 0.99331  p= 6.04e-14
     BS15_row att:4vs3  se_contrast= 0.004473  se_a= 0.056689  se_b= 0.057745  rho= 0.99711  p= 1.02e-18
NEW  new/repro_inference.log
     Caution: the analytic covariance of the subgroup effects is numerically rank-deficient in the
       contrast direction for 5 pairwise row(s) (largest implied correlation between two subgroups'
       estimates 0.9971). The saturated cell dummies leave each cell's residuals with no within-cluster
       variation to identify, so the contrast's standard error is cancellation residue rather than an
       estimated variance - the same degeneracy named in the jointtests caution. The standard error,
       statistic, p-value and interval of those rows, and the across-subgroup equality test(s) covering
       them, are reported as absent rather than computed from that covariance; the point estimates are
       unaffected. wooldid's randomization inference (permute()) or bootstrap-based clustered SEs ...
     console:  p-value for joint test that all SGs have the same att: .
               p-value for joint test that all SGs have the same pre-treatment att: 0.0000
     e(pwsgs) att rows     estimate        se          t           p        lb_95ci   ub_95ci
       2vs1.region          .06572      .           .          .            .         .
       3vs1.region         -.06755      .           .          .            .         .
       4vs1.region         -.003665     .           .          .            .         .
       3vs2.region         -.1333       .009792   -13.61      6.378e-18    -.153     -.1136
       4vs2.region         -.06939      .           .          .            .         .
       4vs3.region          .06389      .           .          .            .         .
     BS15_joint_att_allsgssame=.              BS15_joint_pre_att_allsgssame=7.717e-68
     BS15_sgdegen=5 act=[blank] r=[0.9971]
     BS15_row att:3vs2  se_contrast= 0.009792  se_a= 0.057745  se_b= 0.055532  rho= 0.98581  p= 6.38e-18   (retained)
     (the other five BS15_row lines print se_contrast= . rho= . p= .)
REPRO_BS15_DONE  (both builds. Point estimates, column 1, identical on both builds; the six pre_att rows
                  identical on both builds, e.g. 2vs1 .09069 / .02239 / 4.051 / .0001897. The new run also
                  prints two cautions the old block does not - "reference period for some units forced"
                  and "the cohort variable fips takes 288 distinct values" - which are other HEAD cautions on
                  this config, not part of BS-15.)
```

- **Classification: old-bug-fixed.** Two claims, proven in two places: the `adm_svpsg`
  reproduction above proves the *behavioural* change (v1.3 prints and posts the residue as
  inference with no caution; HEAD blanks columns 2-6 and the covering equality test and says
  why), while the proof that the residue is *wrong* — the same contrast's se 44x larger under
  `robust` than under `cluster(i)` — lives on the sweep's dta_binary fixture (bug_sweep BS-15
  "Reproduction"; the suite's `esptrip`/`plot_sg` rows carry those numbers) and was cited, not
  re-run. The discriminating fact is the sweep verifier's
  reclassification of BS-15 from "judgment" to **econometric** (bug_sweep_2026-09-06.md §
  BS-15): a contrast se that shrinks 44-fold when the VCE is coarsened from `robust` to
  `cluster(i)` is not a number any valid variance estimator produces, so what v1.3 printed
  was wrong output, not a convention. What IS a ruled convention is the remedy — blank
  rather than repair, and the cut `1 - 1e-2` — so new's "right" value is an absence plus a
  caution, not a corrected number; the point estimates never moved.
- **Open questions:** (i) the cut is a judgment call by the commit's own admission
  (`cfac3f5`, held_for_ruling "BS-15 cut"); on `adm_svpsg` the retained row `att:3vs2`
  sits at rho `.98581` (both builds compute the same value) against a cut at `.99` — its
  contrast se `.009792` is 1/5.7 of the level se `.057745`/`.055532` and its p `6.4e-18`
  prints with no per-row caution, while the flagged `att:3vs1` sits at `.99004`, 4e-5 over
  the cut. The sweep's discriminator (the same contrast under `robust`, 44x on dta_binary)
  was NOT run on `adm_svpsg` here, so the retained row is measured, not adjudicated: it is
  either a milder instance of the same residue or a genuine contrast, and this
  reproduction cannot tell which. (ii) The brief's example figure for this config ("the
  analytic contrast se old prints, e.g. 3.7e-9 with p = 0.0000") does not occur on
  `adm_svpsg`: old's att contrast SEs run `.004473`-`.009792` with p `1.0e-18`…`.464`, and
  neither the ledger reason (`BS15-pw-adm_svpsg`) nor the verbose re-judge nor any record
  grepped carries a `3.7e-9`; the "p = 0.0000" is old's console rendering of
  `joint_att_allsgssame = 3.97e-20`. If the parent had another config in mind, it is not
  among the BS15 ledger rows. (iii) On
  `cont_sg_gap` the six att pairwise rows' columns 2-6 (old p's `3.5e-15` … `1.3e-12`)
  are blanked by this same screen (`e(wdid_sgdegen)` = 6, `act` = `blank`), but the ledger
  carries no `BS15-pw-cont_sg_gap` row: the pattern-`pwsgs` row `SIA4-gp-p` (first match)
  absorbs them under a reason that describes the ame rows 7-12 (`0 → .`), not the att rows.
  Attribution imprecision only — no unexplained cell — but per the project's ledger rule a
  class should carry its own named row.

### 2. BS-10: `e(df_r)` under `impvalues` was the stage-1 residual df, not the Ibragimov-Müller reference (`BS10-dfr`)

- **Symptom / impact:** under `impvalues` every printed t, p and interval comes from a
  per-row IM t-test on `t(q-1)`, `q` the clusters contributing a non-missing per-cluster
  estimate to that row; but v1.3 posted the stage-1 regression's residual df — always
  larger — as `e(df_r)`. So `ereturn display`, `lincom`, `test`, `esttab` and `coefplot`
  reported p-values **too small** and intervals **too narrow** relative to the table
  wooldid had just printed, from the same run. On the suite's own `impv` config
  (`dta_binary`, `att cluster(i) esprelength(2) espostlength(2) impvalues`): `e(df_r)` =
  79 (= G-1) while every printed row reproduces only at `t(49)`; the `pretreatment:att`
  row prints p `.4068168` / `[-.1555138, .0640822]` and `lincom` on the same run gives
  `.405` / `[-.1544689, .0630373]`. Widening to a 6x6 window makes the row-dependence
  visible (df 49 on most rows, 24 on the far horizons; `post4` printed `2.475e-09` vs
  `2*ttail(79, 9.19) = 4.08e-14`) — bug_sweep_2026-09-06.md, BS-10.
- **Mechanism (v1.3):** `old/wooldid.ado` L1170 `local df_r = e(df_r)` captured after the
  stage-1 regression; the eight IM sites (L1566, 1597, 1749, 1823, 1852, 2127, 2206, 2245)
  each run `qui ttest `ibmucoefs'1 == 0` and build the row from `r(se)`, `r(t)`, `r(p)`
  and `invt(r(df_t), .025)` — i.e. that row's own `t(q-1)` — but nothing records
  `r(df_t)`; L2987 `cap ereturn post `eb' `ev', ... dof(`df_r')` (and the two rescue posts
  L2989-2990) then post the stage-1 value. No `ereturn scalar df_r` reset; no IM df stored
  anywhere.
- **Fix + proof:** author ruling 2026-09-07 (`rulings_2026-09-07.md` §A, BS-10): **"Post
  the IM degrees of freedom"** — the minimum across rows, so `ereturn display` and
  `lincom` agree with the printed table; **do not post missing** (the normal fallback is
  anti-conservative). This chose the sweep's reading (b) *fix* over (a) *document*. Built
  as `92e083b` (2026-09-10): at the six live IM ttest sites each surviving row appends
  `r(df_t)` to `imdflist` under the exact negation of the degeneracy scrub
  (`new/wooldid.ado` L1972, 2023, 2138, 2238, 2293, 2631); `wooldid_imdf` (L10345, called
  at L3948) writes the **minimum** over `e(df_r)`, keeps the stage-1 value as hidden
  `e(df_r_analytic)` (the bootstrap wire's convention), posts `e(imdf_min)`/`e(imdf_max)`
  and prints a footer when the df is not common across rows (L10895-10896); `e(gteffects)`
  cells stay on the analytic df. Tests `test/test_impvalues.do` I1/I1b/I2/I2b/I3/I5.
  Commit record: the 20-leg fingerprint moved on exactly one line (`e(df_r)` 79 → 49) and
  the battery on exactly one configuration.
- **This run:** `BS10-dfr-impv` on `impv` (`df_r` 79 → 49), `BS10-dfr-trig` on
  `trig_imdeg_part` (29 → **1**), `BS10-dfr-extbin` on `ext_bin_ext_impvalues` (11 → 8),
  `BS10-dfr-extcont` on `ext_cov10_cont_impvalues` (11 → 8). New-only additions
  `df_r_analytic`, `imdf_min`, `imdf_max` on each. No reported table number moved.
- **Reproduction (`impv`):**

```
OLD  old/repro_inference.log  (table row, then the same run's lincom)
  pretreatment att | -.0457158   .0546374  -.8367117   .4068168  -.1555138   .0640822
  BS10_df_r=79
  BS10_df_r_analytic=.
  BS10_table_pretreatment_p=   .4068167917 lb=  -.1555137957 ub=  .06408224385
  lincom pretreatment_att:  (1) | -.0457158  .0546374  -0.84  0.405  -.1544689  .0630373
  BS10_lincom_p=   .4052785219 lb=  -.1544688535 ub=   .0630373016 df=79
  BS10_2ttail_at_e_df_r=   .4052785219 at_t49=   .4068167917 at_t79=   .4052785219
NEW  new/repro_inference.log
  pretreatment att | -.04571578  .05463743  -0.84  0.407  -.15551379   .06408224
  BS10_df_r=49
  BS10_df_r_analytic=79
  BS10_imdf_min=49 imdf_max=49
  BS10_table_pretreatment_p=   .4068167951 lb=  -.1555137947 ub=  .06408224391
  lincom pretreatment_att:  (1) | -.0457158  .0546374  -0.84  0.407  -.1555138  .0640822
  BS10_lincom_p=   .4068167951 lb=  -.1555137947 ub=  .06408224391 df=49
  BS10_2ttail_at_e_df_r=   .4068167951 at_t49=   .4068167951 at_t79=   .405278525
REPRO_BS10_DONE  (both builds)
(the 9th-digit drift .4068167917 vs .4068167951 is the documented svmat storage-noise
 class A9-svmat-w on this config, worst 2.12e-07; the lincom-vs-table gap is the defect)
```

- **Classification: old-bug-fixed.** One run's `ereturn display`/`lincom` contradicted its
  own printed table, in the anti-conservative direction; the author ruled it a fix, not a
  documentation gap. The min-across-rows rule is the ruled convention on top of the fix.
- **Open questions:** (i) the sweep itself offered reading (a) — `e(df_r)` = the regression's
  residual df is Stata's usual contract and the `impvalues` documentation made no promise that
  `lincom` would match the IM table — as a defensible alternative to (b) fix (bug_sweep BS-10,
  "Verifier"); the 2026-09-07 ruling ("Post the IM degrees of freedom ... Do not post missing")
  is the discriminating fact for the old-bug-fixed call, and the min-across-rows rule is
  convention layered on it. (ii) On `trig_imdeg_part` (imbalanced panel, `espostlength(9)`) the rule
  yields `e(df_r)` = 1 (one far horizon is reached by two clusters), so a hand-run `lincom`
  on that run's healthy `att` row (q = 30) is now on `t(1)`. That is the convention working
  as ruled at its extreme — conservative for every row — and the footer plus
  `e(imdf_min)`/`e(imdf_max)` disclose it; noted here because the number is startling.

### 3. CS-60: per-subgroup joint pre-trend p-values were `e()` macros; now scalars — value-compared at `HEAD`, no blanket row needed

> **What this section is judged against (framing rewritten 2026-09-14).** At `HEAD` (`4280f3b`)
> there is **no ledger row for this class, and none is needed**: the comparer now compares the old
> macro against the new scalar **across the kind change**, so these cells are value-compared at
> `--rtol 1e-10` and come back MATCH instead of "a removal plus an informational addition". That is
> exactly the reclassification the headline block reconciles — `399 EXPLAINED` at `26ce421`
> `= 374 EXPLAINED + 25 MATCH` at `4280f3b`, the 25 being the configs on which MATCH became
> reachable once the comparer read these two sides as one value (`4280f3b`, verbatim: *"The comparer
> keys its kinds separately, so 361 values on 103 configs were being blessed as a removal plus an
> informational addition and never value-compared. It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)."* The commit does not
> say which 25, and neither does this section.) An earlier draft was framed on four blanket rows,
> `CS60-jhm-trig` / `-adm` / `-sg` / `-ext`, which `4280f3b` retired for that reason — a reader will
> not find them at `HEAD` and does not need to. What *is* live, and is what a reader should look up,
> is the narrower T-M2 subset: `TM2-sg-p` (`trig_tm2_sgjoint`) and `TM2-ext-p` (`ext_c2_*`), both now
> `kind = scalar` with the pattern narrowed to `joint_Choriz*_s*` (§4 is about those). **Not one
> measurement below changed; only the framing did.** See `LEDGER_DELTA.md` §A and §D.

- **Symptom / impact:** `e(joint_horiz)`, `e(joint_horiz_i)`, `e(joint_Choriz)`,
  `e(joint_Choriz_i)` were `ereturn scalar`, and the adversarial and `allsgssame` families
  were scalars at both the pooled and the subgroup level — but the four subgroup joint
  pre-trend twins `e(joint_horiz_s#)`, `e(joint_horiz_i_s#)`, `e(joint_Choriz_s#)`,
  `e(joint_Choriz_i_s#)` were `ereturn local`. A loop over subgroups got numbers for one
  family and strings for the other; `if e(joint_horiz_s1) < .05` is a **type mismatch
  r(109)** on v1.3; `esttab`/`estout` and any numeric read of e() saw a string. The values
  themselves were identical. The help codified the split ("presented by subgroup, when
  used, in the local macros") — consistency_sweep_2026-09-06.md, CS-60.
- **Mechanism (v1.3):** `old/wooldid.ado` L3146, 3150, 3156, 3161
  `ereturn local joint_horiz_s`sgcase' = ...` (and `_i_s`, `Choriz_s`, `Choriz_i_s`)
  inside the per-subgroup display loop, against L3098/3102/3106/3110
  `ereturn scalar joint_horiz` etc. for the pooled four.
- **Fix + proof:** author ruling 2026-09-07 (`rulings_2026-09-07.md`, CS-60): **"Post the
  subgroup joint p as a scalar at all three posting sites, so a loop over subgroups gets one
  type."** Built as `605d6d7` (2026-09-09, batch 3 "stored results"): `new/wooldid.ado`
  L4177-4180 `ereturn scalar joint_horiz_s`sgcase' = ...` (x4), plus the bootstrap and
  permutation engines' twins (`open_defects.md`: "scalar at all THREE posting sites (the
  audit found two; the permutation engine had the same pair as the bootstrap one)").
  `test/test_stored.do` asserts by regex that no `joint_*_s*` name is still a macro. Commit
  record: cert 40/40, battery re-blessed; every battery difference line was an added name or
  the CS-60 macro-to-scalar swap.
- **This run, as the `26ce421` comparer reported it:** that comparer keyed macro and scalar cells
  separately, so it reported each `joint_horiz[_i]_s#` macro REMOVED and the identically-named
  scalar ADDED with the same value, and four blanket rows blessed the pairs — on
  `trig_tm2_sgjoint` 3 lines (`.1036776642908662`, `8.53032953287e-16`, `1.9389369068e-45`), on
  `adm_base` 4, on `sg_attit`, `sg_degpre`, `sg_missif` 9, and 260 lines across the `ext_b2_*`,
  `ext_b3_*`, `ext_b4_*`, `ext_c2_*` grids. **Zero p-values moved.** The line counts are that
  comparer's bookkeeping; the values are the same values on both sides, which is why `HEAD`'s
  comparer — reading the pair as one value — returns MATCH on these cells and needs no row.

- **Reproduction (`trig_tm2_sgjoint`, shared with §4):**

```
OLD  CS60_joint_horiz_s1_macroform=[.1036776642908662]   CS60_joint_horiz_s1_exprform=.1036776642908662  (a string)
     CS60_scalars_named_joint_horiz_s: no      CS60_macros_named_joint_horiz_s:  YES
     CS60_numeric_compare_result=type mismatch
     CS60_numeric_compare_rc=109              <- `e(joint_horiz_s1) < .05' on v1.3
NEW  CS60_joint_horiz_s1_macroform=[.1036776642908662]   CS60_joint_horiz_s1_exprform=.10367766  (a scalar)
     CS60_scalars_named_joint_horiz_s: YES     CS60_macros_named_joint_horiz_s:  no
     CS60_numeric_compare_result=0
     CS60_numeric_compare_rc=0
REPRO_CS60_DONE  (both builds; s2 = 8.53032953287e-16 and s3 = 1.93893690680e-45 identical)
```

- **Classification: intended-change.** Values bit-identical on both builds; only the e()
  type moved, and the old convention was documented in the help. The r(109) demonstration
  shows why it was worth changing without making the old behaviour wrong output.

### 4. T-M2 / SI-10: the sg-continuous joint pre-trend p was read at `r(p)[1,2]` of a 1x1 matrix — silently missing (`TM2-sg-p`, `TM2-ext-p`)

- **Symptom / impact:** with `contreat()` + `subgroups()` and at least two pre-periods,
  the per-subgroup continuous joint pre-trend p `e(joint_Choriz_s#)` (and `_i_s#`) was
  always `.` — displayed as `p-value for joint test that sg=2 pre-treatment ame horizon
  estimates are 0: .` and posted as the macro `.` — even where the test is estimable and
  decisive; the pooled `e(joint_Choriz)` and the binary subgroup twins were real. A user
  reading the subgroup table saw "no evidence" where sg2's real p is `1.4e-13`. The
  battery's `cont_sg.csv` was the sole baseline pinning the `.`, and because the old build
  shares the defect, old-vs-new comparison alone could never flag it
  (`wooldid_solved_issues.md` L427-431).
- **Mechanism (v1.3):** `old/wooldid.ado` L2617 runs the sg-continuous
  `margins ... `contreatelasticitytype'(`contreat') over(`_t') ... contrast(overjoint)`
  whose `r(p)` for a continuous dy/dx is **1x1** (the joint column only), unlike the binary
  twin's 1x2 (vacuous base column + joint); L2620 then reads
  `local joint_Choriz_s`sgcase' = r(p)[1,2]` and L2631 the same for `_i_s`. Stata's lenient
  out-of-range matrix subscript in an expression returns missing with rc 0 (probed:
  solved_issues L429). The pooled continuous sites L2149/L2154 read `[1,1]` correctly, and
  the binary sg sites L2430/L2440 read `[1,2]` correctly for their 1x2 `r(p)`.
- **Fix + proof:** found as ticket T-M2 by the 2026-08-23 mini-assessment; the author
  challenged it and the challenge was empirically adjudicated (instrumented build on the
  battery `cont_sg` configuration: the continuous overjoint `r(p)` is 1x1 with real p's at
  `[1,1]` — sg1 `6.209e-06`, sg2 `4.109e-07` — that the `[1,2]` read discards); FINAL
  RULING 2026-08-23 "DEFER — address during the full assessment run" with a mandated
  evidence protocol (solved_issues L429-431). Landed 2026-08-25 as `413eabc` (author-
  ordered): evidence banked PRE-fix at `review_process/tm2_evidence/` (3+ datasets incl.
  rich non-degenerate data, instrumented + stock builds); `[1,2]` → `[1,1]` at the two sites
  (`new/wooldid.ado` L3235, L3246); the binary twins untouched; suite config
  `trig_tm2_sgjoint` (planted sg-2/sg-3 dose pre-trends) and ledger rows `TM2-sg-p` /
  `TM2-ext-p` pinned with the expected direction recorded — old `.`, new real; a new-side
  regression to `.` leaves no diff, strands the row STALE and fails the run. Validated
  2026-08-25 on the SI-A1 full run (old `joint_Choriz_s1 = .` vs new `.8140497...`)
  (solved_issues L485-486).
- **This run:** `TM2-sg-p` on `trig_tm2_sgjoint`: `joint_Choriz_s1..3` macro `.` on old →
  scalars `.8140` (sg1, no planted trend, correctly unrejected), `1.41e-13` (sg2),
  `1.44e-36` (sg3); the att twins `joint_horiz_s1..3` (`.1037`, `8.53e-16`, `1.94e-45`)
  and the pooled `joint_Choriz` (`6.17e-22`) are bit-identical across builds. `TM2-ext-p`
  fired on **24 of the 32** `ext_c2_*` configs (`h10`…`h55` x `f1/f2/f4/p`; `h22_f2` and
  `h55_f4` also post the `_i_s#` twins), e.g. `ext_c2_h10_f1` new `s1 1.62e-05`,
  `s2 .603`, `s3 .00112`; `ext_c2_h22_f2` `s1 .0560`, `s2 4.04e-05`, `s3 .230`. The eight
  `h00_*`/`h03_*` configs have `esprelength(0)`/`(1)`, so no sg pre-trend joint test exists
  there and nothing fires — consistent with the row's scope.
- **Reproduction (`trig_tm2_sgjoint`, data generated in-session by the suite's prep):**

```
OLD  console:  p-value for joint test that sg=1 pre-treatment ame horizon estimates are 0: .
               p-value for joint test that sg=2 pre-treatment ame horizon estimates are 0: .
               p-value for joint test that sg=3 pre-treatment ame horizon estimates are 0: .
     TM2_joint_Choriz_pooled=6.175e-22
     TM2_joint_Choriz_s1_macroform=[.]                   TM2_joint_Choriz_s1_exprform=.
     TM2_joint_Choriz_s2_macroform=[.]                   TM2_joint_Choriz_s2_exprform=.
     TM2_joint_Choriz_s3_macroform=[.]                   TM2_joint_Choriz_s3_exprform=.
NEW  console:  ... sg=1 pre-treatment ame horizon estimates are 0: 0.8140
               ... sg=2 pre-treatment ame horizon estimates are 0: 0.0000
               ... sg=3 pre-treatment ame horizon estimates are 0: 0.0000
     TM2_joint_Choriz_pooled=6.175e-22
     TM2_joint_Choriz_s1_macroform=[.814049717037408]    TM2_joint_Choriz_s1_exprform=.81404972
     TM2_joint_Choriz_s2_macroform=[1.41265752805e-13]   TM2_joint_Choriz_s2_exprform=1.413e-13
     TM2_joint_Choriz_s3_macroform=[1.44227544777e-36]   TM2_joint_Choriz_s3_exprform=1.442e-36
REPRO_TM2_DONE  (both builds; the att twins joint_horiz_s1..3 print .1036776642908662,
                 8.53032953287e-16, 1.93893690680e-45 on BOTH builds)
```

- **Classification: old-bug-fixed.** A real, estimable p-value was replaced by missing
  through an out-of-range subscript; the author's own challenge was disproven on the
  evidence and the fix was author-ordered.

### 5. DDG: pairwise-subgroup cells missing on v1.3, valued on HEAD — on a spec-degenerate design (`DDG-pw-1..3`, `DDG-jt-1..3`)

- **Symptom / impact (as recorded):** the 2026-08-23 report §1.13 records that on the
  RCS + pweight / twoway-cluster + `subgroups()` + `fe(group_time)` + `safety(off)`
  designs (`ext_cov10_rcs_*`), v1.3's double-margins pairwise machinery returned MISSING
  `e(pwsgs)` cells (the `%21x` encoding of Stata missing, ticket 1 of REPORT_2026-08-18)
  and missing `*_allsgssame` scalars, plus one garbage negative se (`-8.9e-16`), where
  HEAD's E1 posted-pwcompare + rank guard (`cd8ce95`/`2a1c9ca`/`7ea1533`/`184a6ff`/
  `ed780ac`, 2026-07-26; solved_issues L310) computes values.
- **Mechanism (v1.3):** the same L2459-2481 block as §1: `contrast(overjoint)` for the
  homogeneity p and a second `pwcompare(effects)` pass for the table, each a full
  delta-method margins pass with no post/restore; where the pass returned an empty or
  rank-deficient `r(table_vs)`, the sliced cells are missing. HEAD runs one posted
  `pwcompare` pass, reads the homogeneity p from a chain `test` on the posted `e(b)`
  (`wooldid_sgchain`), and falls back to the contrast pass when `r(df) != k-1`
  (`new/wooldid.ado` L2990-3012, pre leg L3080-3101).
- **What this run actually shows:** `DDG-pw-1`/`DDG-jt-1` on
  `ext_cov10_rcs_pweight_ddd_c1`: `joint_att_it_allsgssame` `.` → `2.27e-199` **and**
  `joint_pre_att_allsgssame` `.99999787` → `.` (real → missing); `pwsgs` 25 cells, e.g.
  main `att:2vs1` old `[1, ., ., ., ., .]` → new est `1.0`, se `2.33e-10` (= 2^-32), t `4.29e9`
  (= 2^32), p `2.1e-207`, CI `[1, 1]` — the row v1.3 had carried under `att_i`, which on HEAD
  goes missing in its turn. `DDG-pw-2`/`jt-2` on `..._ddd_c2`: old degenerate → new
  differently degenerate (`joint_att_allsgssame` `1.85e-198` → `6.38e-202`; pwsgs t
  `1.75e9` → `2.48e9`). `DDG-pw-3`/`jt-3` on `..._twoway_cluster_c1`: `pwsgs` 10 cells
  `.` → the same `2.33e-10 / 4.29e9 / 2.1e-207` row, and again `joint_pre_att_allsgssame`
  `.99999787` → `.`. HEAD's own BS-15 screen (§1) flags **8 of 8** pairwise rows on the
  `_c1` config with "largest implied correlation `2.0000`" — a correlation above 1, i.e.
  the polarisation identity itself is broken by the variances in hand — but under
  `safety(off)` it only reports (`e(wdid_sgdegenact)` = `report`), so nothing is blanked.
  A subgroup contrast of exactly `1.0` with se `2e-10` is `fe(group_time)` absorbing the
  cells so the treatment effect becomes a deterministic identity (the SI-7 knife-edge the
  `survey_rcs` data-class note attributes to the `fe(i.group##i.t)` spec, "spec-induced").
- **Reproduction (`ext_cov10_rcs_pweight_ddd_c1`):**

```
OLD  (no caution printed anywhere on the run)
     e(pwsgs)     estimate         se           t            p       lb        ub
       att          1            .            .            .
       pre_att     -8.882e-16    3.293e-10   -2.697e-06    1      -6.8e-10   6.8e-10
       att_it       1            .            .            .
       pre_att_it   0            .            .            .
       att_i        1            2.328e-10    4294967296   2.080e-207   1        1
       pre_att_i   -8.882e-16    1.646e-10   -5.395e-06    1
       att_itime    1            .            .            .
       pre_att_it~e -8.882e-16   1.646e-10   -5.395e-06    1
     DDG_joint_att_allsgssame=.   joint_pre_att_allsgssame=.99999787
     DDG_joint_att_it_allsgssame=.   joint_pre_att_it_allsgssame=.      (att_i, att_itime: . / .)
     DDG_sgdegen=. act=[] r=[] fcguard=[]
NEW  Caution: wooldid found subgroup contrasts with zero variance at a subgroup margins call.
       Wooldid treats that as a numerically singular direction of the VCE, where an exact zero
       margins returns can be cancellation residue rather than an estimate, so which subgroup
       horizons and joint tests are reported can depend on the summation order - which the
       default cell collapse (see fastcollapse in the help file) changes. ... with safety(off)
       it proceeds on the collapsed rows as requested. ...
     Caution: 4 reported standard errors appear to be identical under a cluster-type VCE. This
       usually indicates a failure of analytic clustered standard error estimation; ...
     Caution: the analytic covariance of the subgroup effects is numerically rank-deficient in
       the contrast direction for 8 pairwise row(s) (largest implied correlation between two
       subgroups' estimates 2.0000), so those contrast standard errors are cancellation residue
       rather than estimated variances. safety(off) was specified, so they and the
       across-subgroup equality tests are reported rather than withheld; ...
     e(pwsgs)     estimate         se           t            p       lb        ub
       att          1            2.328e-10    4294967296   2.080e-207   1        1
       pre_att     -8.882e-16    2.328e-10   -3.815e-06    1      -4.8e-10   4.8e-10
       att_it       1            5.206e-10    1920767767   2.271e-199   1        1
       pre_att_it   0            .            .            .
       att_i        1            .            .            .
       pre_att_i   -8.882e-16    1.646e-10   -5.395e-06    1
       att_itime    1            2.328e-10    4294967296   2.080e-207   1        1
       pre_att_it~e -8.882e-16   1.646e-10   -5.395e-06    1
     DDG_joint_att_allsgssame=.   joint_pre_att_allsgssame=.
     DDG_joint_att_it_allsgssame=2.27e-199   joint_pre_att_it_allsgssame=.   (att_i, att_itime: . / .)
     DDG_sgdegen=8 act=[report] r=[2.0000] fcguard=[caution]
REPRO_DDG_DONE  (both builds)
(the "recovered" row is se = 2^-32 exactly, t = 2^32 exactly: a floating-point residue; on v1.3
 the identical row sat under att_i, on HEAD under att and att_itime while att_i went missing)
```

- **Classification: unclear.** (**The disagreement with `legacy` L.8 item 5 is settled, in this
  section's favour, 2026-09-14.** That item used to hand the DDG family off as "re-verified by
  another agent" and therefore settled; it now adopts *unclear* and defers here. The settlement was
  made from the primary sources — the two builds, the 2026-09-13 `run_old`/`run_new` dumps on the
  same pinned pair (old `301b706`; new `6a05127`, whose ado is byte-identical to this section's pin
  `08ec0b4`, md5 `382c1e0e…`) and the ledger row — and every measured claim below reproduced
  exactly: `joint_pre_att_allsgssame` `+1.ffffb8906f474X-001` (= `.999997871047198`) →
  `+1.0000000000000X+3ff` (= `.`) on **two** configs; `joint_att_it_allsgssame` `.` →
  `+1.1630ccf698a1aX-294` (= 2.27e-199); and the `2^-32 / 2^32 / 2.1e-207` residue row
  **byte-identical** across builds, sitting under `att_i` on old and under `att` *and* `att_itime`
  on new while new's `att_i` goes missing. Census on `..._pweight_ddd_c1`'s 48 `pwsgs` cells: 15
  missing → value, **5 value → missing**, 20 bit-identical. See `CRITIQUE_CLOSEOUT.md` item 1.)
  The old → new movement is real and traces to the E1 chain
  machinery replacing two margins passes, but §1.13's framing "old returned missing where
  HEAD returns valid values" does not describe these cells: the "recovered" values are
  numerically degenerate (t-statistics of order 1e9 on contrasts of exactly 1.0, p `1e-207`)
  that v1.3 had left missing, sibling configs move from one degenerate value to another,
  and one equality p moves the other way (real → missing) on two configs. Neither build's
  numbers on this spec are meaningful inference; by the pathological-data effort cap this
  is "observable failure?" evidence, not feature-quality evidence, and was not diagnosed
  further here.
- **Open questions:** (i) the `DDG-jt-1`/`DDG-jt-3` `joint_*` wildcard absorbs a
  real → missing move (`joint_pre_att_allsgssame` `.99999787` → `.`) under a reason that
  says "recovered"; its cause on HEAD is not BS-15 (report mode there) — HEAD's console on this block prints, before the table, the cell-collapse degeneracy caution (`wooldid found subgroup contrasts with zero variance at a subgroup margins call ... which subgroup horizons and joint tests are reported can depend on the summation order - which the default cell collapse changes ... with safety(off) it proceeds on the collapsed rows as requested`; `e(wdid_fcguard)` = `caution`), i.e. the F-collapse guard `wooldid_fcdeg` let the subgroup tests be evaluated on collapsed rows under `safety(off)`, and the identical `2^-32 / 2^32 / 2.1e-207` residue row that v1.3 carried under `att_i` appears under `att` and `att_itime` on HEAD while `att_i` goes missing — residue relabelled across rows, not recovered. Read from the console only; not instrumented.
  **[Console reading CONFIRMED against the dumps 2026-09-14]** — the relabelling this question
  reads off the console is now measured cell by cell in the `%21x` dumps: old's `att_i` row is
  byte-identical to new's `att` and `att_itime` rows, and five of new's `att_i` cells are missing.
  **[PARTIALLY CLOSED by `4280f3b`, 2026-09-13]** — the `joint_*` half of the complaint is fixed:
  `DDG-jt-1`'s reason now records both directions explicitly (*"this one row blesses both
  directions; the withheld move is recorded here not re-diagnosed"*) and `DDG-jt-3` says only the
  withheld direction fires. Retained, not deleted: the question as written was that the wildcard
  absorbed the move under a reason that said "recovered".
  (ii) The class deserves a ledger reason that says "spec-degenerate on both builds; cells
  move" rather than "missing → valid". **Still open on the `-pw` rows.** `4280f3b` gave the three
  `-jt` reasons the two-direction wording but left the `-pw` framing as "recovered":
  `DDG-pw-1`'s reason is **byte-unchanged** from `26ce421` and still reads *"v1.3 returns
  MISSING pairwise-SG cells … where HEAD returns valid values"*, while `DDG-pw-2` / `-pw-3` were
  merely expanded from "same" to *"as DDG-pw-1 … missing-cells-recovered class"*. On that same
  config five
  `pwsgs` cells move value → missing and the suite's own degeneracy scan of the **new** dump
  records `PWSGS_DEAD(0v+2i/8)` against `degeneracy_expected.csv`'s RULED row citing `DDG-pw-1`
  (line 173). Corroborating, not contradicting, this section's *unclear*.

### 6. E1b / A18: the whole-sample joint pre-trend tests ignored the requested unconditional VCE (`A18-jt`)

- **Symptom / impact:** the two binary whole-sample pre-trend joint tests (`joint_horiz`,
  `joint_horiz_i`) hardcoded `df()` and silently used the conditional VCE even when the
  user specified `unconditionalse` — the reported joint p was computed under a VCE the
  user did not ask for (2026-08-23 report §1.12).
- **Mechanism (v1.3):** `old/wooldid.ado` L1521 `local marginsdof df(`df_r')` and L1524
  `if "`unconditionalse'" == "unconditionalse" local marginsdof `marginsdof'
  vce(unconditional)` — every estimand margins call uses `marginsdof`, but the two joint
  tests at L1773/L1777 end in a literal `df(`df_r')` instead, so `vce(unconditional)`
  never reaches them.
- **Fix + proof:** `7ea1533` (2026-07-26, separable fold-in of the E1 chain): a dedicated
  `jointmarginsdof` local (`new/wooldid.ado` L1880-1881, used at L2170/L2174) threads
  `vce(unconditional)`; deliberately not `marginsdof`, whose `impvalues` value would kill
  the test. Acceptance `test/test_dffix.do`, tolerance-based because the unconditionalse
  path is not bit-reproducible run-to-run (solved_issues L310).
- **This run:** `A18-jt` on `trig_uncond_jt` (`[aw=w] att robust unconditionalse
  jointtests`): `joint_horiz` `.0993522568809689` → `.1269868572066457`.
- **Reproduction (`trig_uncond_jt`):**

```
OLD  console: p-value for joint test that all pre-treatment att horizon estimates are 0: 0.0994
     A18_joint_horiz=  .09935225688
NEW  A18_joint_horiz=   .1269868572
REPRO_A18_DONE  (both builds; unconditionalse path, not bit-reproducible run-to-run -
                 both values reproduce the 2026-09-11 dumps to the printed digits)
```

- **Classification: old-bug-fixed.**

### 7. IM degeneracy guard and the e(V)-preserve post (`A16-imd-*`, `A17-*`, `A17-part-*`)

- **Symptom / impact:** (A16) with `impvalues` on noiseless / near-deterministic outcomes,
  v1.3 reported a collapsed se of exactly 0 with a CI of zero width and posted a
  zero-diagonal `e(V)` — a storage artifact presented as certainty (2026-08-23 report
  §1.9). (A17) whenever any estimand's se came back missing — a far event-study horizon
  reached by a single cluster under `impvalues`, missing per-subgroup VCV diagonals on the
  continuous-subgroup configs, or plain continuous runs with missing `*_ame` SEs — v1.3
  silently posted **b-only** (no `e(V)`, `e(properties)` = `b`) with no caution, so
  `lincom`/`test`/`esttab` downstream lost the VCV for every healthy estimand too (§1.10).
- **Mechanism (v1.3):** the eight IM ttest sites (L1566 …) take `r(se)` with no
  degeneracy guard, so a numerically identical per-cluster vector gives se 0 (A16) and a
  single contributing cluster gives se `.` (A17); the post ladder L2987-2990
  (`cap ereturn post `eb' `ev'` → on failure `cap ereturn post `eb'` → `ereturn post`)
  accepts a zero diagonal but not a missing one, and falls to b-only in silence.
- **Fix + proof:** `16263da` (2026-07-22): at the six live IM sites, `se < 1e-10*|b|`
  nulls se/t/p/CI to missing with one Caution; because `ereturn post` rejects a missing
  diagonal, `e(V)` is dropped for a fully-degenerate run (`properties` `b V` → `b`);
  verified 3/3 with 22/22 realistic configs bit-identical, zero false fires (solved_issues
  L236). `47defa0` + doc `c53e25b` (2026-07-26): partially-degenerate rows get a **negative
  sentinel** (`-1e100`) V diagonal that renders as missing on every consumer surface,
  healthy rows keep their true VCV, and a Caution reports the affected count; zero-diagonal
  and identified-rows-only alternatives proven harmful; 3+3 unanimous, 32/32 clean configs
  bit-identical (solved_issues L262).
- **This run:** `A16-imd-V/Vn/w/p/b` on `trig_imdegen` (pathological by class): `e(V)`
  REMOVED, `properties` `[b V]` → `[b]`, `wdidmainresults` 13 cells (se `0.0` → `.`, CI
  bounds → `.`), b at svmat storage noise `3.01e-08` (bound 1e-6); `A3-B4b` rides on the
  same config (§8). `A17-part-p/w/b` on `trig_imdeg_part` (stress class, `espostlength(9)`):
  `properties` `[b]` → `[b V]`, `wdidmainresults` 73 cells, b worst `2.41e-06` within the
  1e-4 bound; `A17-sgv-1/2` on `cont_sg_tidose` / `cont_sg_gap`, `A17-uz` on
  `ext_cont_ext_untreatedzero`, `A17-cd` on `ext_path10_cont_constant_dose`: all
  `[b]` → `[b V]`.
- **Reproduction (`trig_imdegen`, `trig_imdeg_part`):**

```
trig_imdegen (A16)
OLD  A16_properties=[b V]
     A16_main_att_est=   1.826923132 se=0 t=. p=. lb=1.8269231 ub=1.8269231
     A16_eV_confirm_rc=0   A16_eV_11=0        (zero-diagonal e(V) posted)
     A16_ll=         1.875                    (B4: an IM bound, see section 8)
NEW  Caution: for one or more Ibragimov-Muller estimates the per-cluster estimates were
       numerically identical, so the cross-cluster (IM) standard error is not identified;
       wooldid has set the standard error, test statistic, p-value, and 95% confidence
       interval to missing for the affected estimate(s) while retaining the point estimate. ...
     Caution: standard errors were not identified for 4 of 4 reported estimates, so no
       variance-covariance matrix e(V) was posted. Point estimates remain in e(b); ...
     A16_properties=[b]
     A16_main_att_est=   1.826923077 se=. t=. p=. lb=. ub=.
     A16_eV_confirm_rc=111                    (no e(V))
     A16_ll=     4106.6423
REPRO_A16_DONE  (both builds)
trig_imdeg_part (A17)
OLD  (no caution)
     A17_properties=[b]   A17_eV_confirm_rc=111   A17_df_r=29
     A17_table post8   est=   1.13442 se=         . p=         .
NEW  Caution: standard errors were not identified for 1 of 13 reported estimates. e(V)
       contains placeholder entries for the affected estimates that display as missing
       standard errors; all other estimates retain their variances in e(V). ...
     A17_properties=[b V]   A17_eV_confirm_rc=0
     A17_Vdiag main_att      = .12760265   ...   A17_Vdiag ES_att_post7 = .56553546
     A17_Vdiag ES_att_post8  = -1.00e+100         (the 47defa0 sentinel)
     A17_df_r=1 df_r_analytic=29 imdf_min=1 imdf_max=15      (BS-10 footer printed:
       "Ibragimov-Muller degrees of freedom are not common across the rows above: they
        run from 1 to 15 ... e(df_r) carries the minimum, 1 ...")
REPRO_A17_DONE  (both builds)
```

- **Classification: old-bug-fixed** (both).
- **Open questions:** A17 is simultaneously a deliberate contract change — partially-degenerate
  runs now post `e(V)` with `-1e100` sentinel diagonals plus a Caution instead of nothing — and
  the 08-23 report lists it under intended changes as well (§3.6), by design; a reader who
  treats "e(V) posted with a placeholder" as the convention rather than "b-only in silence"
  as the bug will file it under intended-change. The silent b-only post with no caution is
  what makes it a defect here.

### 8. B4: `impvalues` clobbered the returned `e(ll)` with an IM confidence bound (`A3-B4*`, `B4-ext-*`)

- **Symptom / impact:** with `impvalues`, `e(ll)` silently returned the **last
  Ibragimov-Müller CI lower bound** instead of the log-likelihood — any model comparison or
  information criterion built on `e(ll)` was garbage (2026-08-23 report §1.2).
- **Mechanism (v1.3):** `old/wooldid.ado` L1287 `local ll = e(ll)` stores the stage-1
  log-likelihood; the eight IM blocks then reuse the same local for the CI lower bound,
  `local ll = `b' + invt(r(df_t),.025)*`se'` (L1568, 1599, 1751, 1825, 1854, 2129, 2208,
  2250); L3070 `foreach scalar in N Nclus ll ll_0 ... { ereturn scalar `scalar' = ... }`
  posts whichever bound was computed last. On `impv` that is the `post1` row's lower bound:
  old's table prints `post1 lb_95ci = .6118055` and old's `e(ll)` = `.611805501489`.
- **Fix + proof:** `4de841d` (2026-07-20), IM locals renamed; battery showed exactly one
  changed value across 18 configs — `e(ll)` under `impvalues`, now the true log-likelihood
  (solved_issues L35-42).
- **This run:** `A3-B4` on `impv` (`.6118` → `-400.91`), `A3-B4b` on `trig_imdegen`
  (`1.875` → `4106.64`; positive ll is correct on that near-noiseless panel), `A3-B4c` on
  `adm_impv` (`.4574` → `-38745.48`), `B4-ext-1` on `ext_bin_ext_impvalues` (`1.0356` →
  `487.18`), `B4-ext-2` on `ext_cov10_cont_impvalues` (`-.0244` → `1684.31`).
- **Reproduction (`impv`, shared with §2):**

```
OLD  B4_ll=        .611805501489    <- old's own table above it: post1 lb_95ci = .6118055
NEW  B4_ll=       -400.910123585
REPRO_B4_DONE  (both builds)
```

- **Classification: old-bug-fixed.**

---

### Cross-check of the reproduction dumps against the suite dumps

`test/compare.py` takes two DIRECTORIES of `%21x` dumps (an earlier file-path call in this
session compared empty sets and printed a vacuous `TOTAL_DIFF_GROUPS=0`; it was redone as
below). `suite_old/` and `suite_new/` hold copies of the seven 2026-09-11 suite dumps
`test/suite/out/run_old|run_new/<config>.csv`; `old/` and `new/` hold this reproduction's
dumps. The harness's `dumpres` omits the `matnames` and `varlist` kinds and the suite's
`esamp`/`files`/`prov`/`end` bookkeeping, so those show as REMOVED; every `scalar`, `macro`
and `mat` cell is compared.

```
python3 test/compare.py suite_old old --rtol 1e-12  -> xcheck_old.txt
  TOTAL_DIFF_GROUPS=39 = 32 "matnames ...: REMOVED" + 7 "varlist: removed=[...] added=[]"
  (adm_svpsg 9, ext_cov10_rcs_pweight_ddd_c1 7, impv 4, trig_imdeg_part 3, trig_imdegen 4,
   trig_tm2_sgjoint 8, trig_uncond_jt 4) - no scalar, macro or mat line
python3 test/compare.py suite_new new --rtol 1e-12  -> xcheck_new.txt
  TOTAL_DIFF_GROUPS=95 = 88 "matnames ...: REMOVED" + 7 "varlist: removed=[...] added=[]"
  (adm_svpsg 24, ext_cov10_rcs_pweight_ddd_c1 22, impv 5, trig_imdeg_part 5, trig_imdegen 4,
   trig_tm2_sgjoint 28, trig_uncond_jt 7) - no scalar, macro or mat line
kind counts, adm_svpsg: suite_new 30 scalar / 23 macro / 742 mat (+23 matnames, 1 varlist, 13 esamp,
  files, prov, rc, end); new/ 30 scalar / 23 macro / 742 mat / rc. suite_old 14 / 5 / 632 (+8 matnames ...);
  old/ 14 / 5 / 632 / rc.
```

So on all seven configs, on both builds, every e() scalar, macro and matrix cell this
reproduction produced is identical (rtol 1e-12) to the 2026-09-11 suite dump the ledger rows
and the verbose re-judge were written against — including the `unconditionalse` config
`trig_uncond_jt`, for which the brief's tolerance exception (not bit-reproducible run-to-run)
did not in fact bite at 1e-12 here. The magnitudes quoted under "This run" above are therefore
the same numbers the reproduction prints.

### Relevant but untested here

- `BS15-*` on `esptrip`, `plot_sg`, `cont_sg_tidose`; `BS10-dfr-extbin/extcont`;
  the `adm_base` / `sg_*` / `ext_*` CS-60 cells (rows `CS60-jhm-adm/sg/ext` at `26ce421`, retired at
  `4280f3b`); `TM2-ext-p` on the 24 `ext_c2_*` configs; `DDG-*-2/3`;
  `A17-sgv/uz/cd`; `A3-B4c`, `B4-ext-*` — cited from the 2026-09-11 dumps and the verbose
  re-judge, not re-run (one representative config per family was re-run, per the brief).
- The BS-15 cut's behaviour on designs between rho `.99` and `.997` (the retained
  `adm_svpsg` row) was measured, not adjudicated; the sweep's robust-vs-cluster discriminator
  was not run on `adm_svpsg`, and the 44x figure is cited from the dta_binary fixture.
- The DDG pre-att equality p's real → missing cause on HEAD (§5) was read from HEAD's
  console output only, not instrumented. **Updated 2026-09-14:** the *movement* is no longer
  console-only — the real → missing pair and the byte-identical residue-row relabelling are now
  measured cell by cell in the `%21x` dumps of the 2026-09-13 run on the same pinned build pair
  (`CRITIQUE_CLOSEOUT.md` item 1). What remains console-only is the attribution of the cause to
  the F-collapse guard `wooldid_fcdeg` under `safety(off)`; that was not instrumented.


---


# Part 6 - intended changes, storage-noise tiers and the rc-by-design class


> ## Headline of record and evidence provenance (added 2026-09-13; identical in all six sections)
>
> **Superseding note added 2026-09-14.** Two things below have since moved and are corrected
> here rather than rewritten in place, because the block is accurate about the run it describes.
> (i) `4280f3b` is no longer `HEAD`; the report's own header names the run of record by its
> `prov[new]` sha, and **that** is the build the report's claims attach to. (ii) The closing run
> of this assessment is the **2026-09-14** one over a **559-config** universe against the same
> v1.3 old side, reported in Part 7.5: **559 of 559 compared, zero DIFFERS, unexplained 0,
> stale_ledger 0**. The 2026-09-11 run described below is where this section's EVIDENCE was
> gathered and its magnitudes measured; the 2026-09-14 run is what verifies that nothing still
> differs. Neither supersedes the other - they answer different questions.
>
> **Run of record.** The 2026-09-11 dual-build run: **old = `301b706`** (published v1.3),
> **new = `14c30d3`**, whose `wooldid.ado` is byte-identical to this section's pin `08ec0b4`
> and to current `HEAD` apart from a nine-line comment block and a shortened `dis as text`
> permutation caution (10 insertions, 1 deletion; no estimation code).
>
> **The report judges those dumps against the ledger and comparer at `4280f3b` = `HEAD`.** That
> suite defines a **538-config universe**; the 2026-09-11 dumps cover **507** of them:
> **507 compared — EXPLAINED = 374, MATCH = 25, RC_EXPLAINED = 108 — unexplained 0,
> stale_ledger 0, unscoped 0, prov_fail 0, COVERAGE_OK.** The other 31 configs are new at
> `4280f3b`, have no 2026-09-11 dumps, and are outside every claim in this report.
>
> An earlier headline — **507/507 compared, EXPLAINED = 399, RC_EXPLAINED = 108**, at ledger
> `26ce421` — is quoted by three of the six drafts (`legacy` and `intended` give the full
> 399/108 line; `inference` gives the 507/507). It is the **same dumps and the same 507 configs**
> judged by the pre-`4280f3b` comparer. The two reconcile exactly: `399 = 374 + 25`, the 25 being
> the configs on which MATCH became reachable once the comparer began comparing across batch 3's
> macro→scalar kind change (`4280f3b`, verbatim: *"It now compares across the kind change; six
> blanket ledger rows retire and MATCH becomes reachable again (25 configs)"*; the same paragraph
> puts the cells affected at *"361 values on 103 configs"*). **Both numbers are measured and
> neither is withdrawn.** See `LEDGER_DELTA.md` §0.
>
> **Evidence provenance.** Every `test/suite/out/run_old|run_new/<config>.csv` cited below is a
> **regenerable** artifact, not a stored one: `test/suite/out/` is gitignored (`test/.gitignore`
> line 5) and the 2026-09-11 dumps were overwritten in place by the 2026-09-13 run. Regenerate
> them with `./run_full.sh oldnew` from `test/suite`, which stages old from
> `git show main:wooldidinstall/wooldid.ado` and new from `git show HEAD:…`. `compare_verbose_0911.txt`
> — cited by four sections, twice as `scratchpad/compare_verbose_0911.txt`, which is **not**
> repo-relative and never was — is likewise a recipe: rebuild it with `compare_oldnew.py … --verbose`,
> passing `--ledger`/`--configs` from whichever commit's judgement you want. The command, the build
> table and the one harness caveat (`4280f3b` changed `mkrun.py`'s `dumpres`, so a regenerated dump
> may differ in its `files|` line or gain a `matnames|…|UNREADABLE` line — **value cells are
> unaffected**) are in **`review_process/assessment_2026-09-13/repro/README.md`**.
>
> **Reproductions.** The do-files this section cites "in this directory" are now parameterised —
> repo root and scratch workdir as positional arguments — and live at
> **`review_process/assessment_2026-09-13/repro/<section>/`**. Their exact invocations are in each
> file's header and in that README. Bulk artifacts (dumps, logs, copied `.ado` builds) are not
> copied: they are regenerable, and the README says how.
>
> **Per-edit log:** `EDITS.md`. **Ledger delta and closed open questions:** `LEDGER_DELTA.md`.


Group: intended-change / storage-noise / rc-by-design. Sections I.1–I.8 below follow the
template of `wooldid_old_bug_report.md` (§0, §1.1–1.5, §1.14, §1.16, §3): symptom / impact,
mechanism at v1.3 (line numbers located in the old build itself), fix + proof, this run's
ledger firings with magnitudes, a runnable reproduction with its decisive lines, and a
classification.

**Builds.** old = published v1.3 (`git main`, `301b706`); new = HEAD `08ec0b4`. Both copies
of `wooldid.ado` live in this directory (`old/`, `new/`) and were run in SEPARATE Stata
processes (`adopath ++` its own directory each).

**Sources of record.** The 2026-09-11 dual-build dumps `test/suite/out/run_old`
(prov `old|301b706`) and `run_new` (prov `new|14c30d3`); the verbose re-judge
the verbose re-judge **at ledger `26ce421`** (regenerable — see the headline block above;
the filename `scratchpad/compare_verbose_0911.txt` is a session path, not repo-relative)
(507/507 compared: **399 EXPLAINED, 108 RC_EXPLAINED
= 66 rc-by-design + 42 both-error parity, 0 unexplained, 0 stale, 0 unscoped, 0 provenance
failures**); the divergence ledger `test/suite/ledger.csv` (row ids cited throughout);
`wooldid_solved_issues.md`; `review_process/rulings_2026-09-07.md` and
`review_process/held_for_ruling_2026-09-09.md` (PS-62). The ado diff `14c30d3..08ec0b4` is one
comment block and one permute caution string (non-numeric), and my own `08ec0b4` dumps equal
`run_new` at rtol 0 on all seven reproduced configs (the `basic` e(V) verified cell-by-cell in
%21x); my old dumps equal `run_old` on six, and on `trig_cfx` differ only in the margins
tempvar row stripes (§I.3's volatility itself). The 0911 verbose therefore stands for HEAD.

**Reproductions.** `repro_intended.do` and `repro_svpsg.do` are at
`review_process/assessment_2026-09-13/repro/intended/`, parameterised on repo root and scratch
workdir. `repro_intended.do` (basic, cov_icbad, cov_icbad2, pois, allest_aw, impv,
trig_cfx, four rc probes; run as `cd old; stata-mp -e do ../repro_intended.do old`, then the
same with `new`) and `repro_svpsg.do` (adm_svpsg). Logs: `old/repro_intended.log`,
`new/repro_intended.log`, `old/repro_svpsg.log`, `new/repro_svpsg.log`; %21x dumps in
`old/dumps/`, `new/dumps/`. Every block prints a `REPRO_<id>_DONE` sentinel; all sentinels
fired on both builds. The decisive lines below are pasted verbatim from those logs.

---

### I.1 `e(V)` now carries the within-family event-study covariance blocks (PS-62, ledger `PS62-esblock`) — INTENDED CHANGE

- **Symptom / impact (old behaviour, not a defect in the sense of a wrong number):** v1.3
  posted `e(V)` as `diag(se^2)` — every off-diagonal exactly 0. `test` and `lincom` across two
  event-study horizons therefore ANSWERED, using an independence covariance wooldid never
  claimed; the 2026-09-06 practice sweep flagged this (PS-62: "`e(V)` beside them is
  deliberately false"). Old also has no `e(V_es_*)` at all (those are F2 additions,
  `d4b2118`), so nothing in old's `e()` carried the horizon covariance that its own margins
  call had computed.
- **Mechanism (v1.3):** old line 2976
  `` matrix `ev' = diag(hadamard(`ebformer'[...,2],`ebformer'[...,2])) `` squares the
  standard-error column of the assembled results table, and line 2987
  `` cap ereturn post `eb' `ev', ... `` posts it. New: `wooldid_vsplice` (new lines
  10557–10608, Mata `wooldid_vspl` at 16416) copies each family's `e(V_es_<family>)` into the
  matching diagonal block of `e(V)` through `ereturn repost V=`, behind a per-family gate (every
  positive-finite-variance horizon must match a row by name AND agree with the V_es diagonal to
  1e-8 relative); `e(wdid_vcvstate)` = `analytic-blocks`, `e(wdid_vblocks)` names the families
  spliced. Cross-family and family-vs-pooled off-diagonals stay zero by ruling (three modes —
  impvalues, unconditionalse, fastnose — have no common delta-method object, and a pooled
  effect is a linear combination of its own horizons, so the general matrix would be
  rank-deficient by construction).
- **Fix + proof:** `a366c46` (2026-09-10), implementing the author's 2026-09-07 rulings after
  the six-lens feasibility review (`review_process/rulings_2026-09-07.md` §F: scope = ES blocks
  only; labelling = help text plus stored `e()` macro, deliberately no footer). The one
  deviation from the literal ruling — the −1e100 sentinel diagonal is KEPT and only its
  off-diagonals held at zero, because `ereturn repost` was measured to accept an indefinite
  matrix verbatim and zeroing would turn a missing SE into a fabricated 0 — is written up in
  `review_process/held_for_ruling_2026-09-09.md` and was confirmed "keep as built". PS-62 is a
  ruling record, not a `wooldid_solved_issues.md` entry.
- **This run:** `PS62-esblock` fired on **335 configs, 18,554 cells**, every one
  "event-study within-family off-diagonal cell(s) 0 -> nonzero" in `matrix V` (basic: 30 cells
  = the 6×5 off-diagonals of one 6-horizon family; adm_svpsg: 10). No other object is touched by
  the row. The ledger row text says "18554 cells over 360 configs on the 2026-09-10 dual-build
  run": the cell count matches this run exactly, the config count does not (open question
  below).
- **Reproduction (`repro_intended.do`, block `basic`: `y i t ttre, att cluster(i)
  esprelength(3) espostlength(3) summarystats jointtests safety(off)` on `test/data/dta_binary.dta`):**
  ```
  REPRO_PS62 build=old e(V) is 8x8; names: main_att pretreatment_att ES_att_pre3 ES_att_pre2 ES_att_pre1 ES_att_post0 ES_att_post1 ES_att_post2
  REPRO_PS62 build=old e(wdid_vcvstate)=[] e(wdid_vblocks)=[]
  REPRO_PS62 build=old diag[1,1] = +1.52a96a2b6eba4X-008  ( .0051675686076419)
  REPRO_PS62 build=old diag[3,3] = +1.46d6d15b60761X-008  ( .0049871693730682)
  REPRO_PS62 build=old diag[8,8] = +1.6f38b9b4633b0X-007  ( .0112067133675512)
  REPRO_PS62 build=old ES block rows/cols 3-8: sum|offdiag| =                     0, nonzero offdiag cells = 0
  REPRO_PS62 build=old outside the ES block (main_att / pretreatment_att rows+cols): sum|offdiag| =                     0, nonzero = 0
  REPRO_PS62 build=new e(wdid_vcvstate)=[analytic-blocks] e(wdid_vblocks)=[att]
  REPRO_PS62 build=new diag[1,1] = +1.52a96a2b6eba4X-008  ( .0051675686076419)
  REPRO_PS62 build=new diag[3,3] = +1.46d6d15b60761X-008  ( .0049871693730682)
  REPRO_PS62 build=new diag[8,8] = +1.6f38b9b4633b0X-007  ( .0112067133675512)
  REPRO_PS62 build=new ES block rows/cols 3-8: sum|offdiag| =     .0652772466737125, nonzero offdiag cells = 30
  REPRO_PS62 build=new outside the ES block (main_att / pretreatment_att rows+cols): sum|offdiag| =                     0, nonzero = 0
  ```
  All eight diagonal cells are %21x-identical across builds (the full eight lines are in both
  logs). Offline check on the %21x dumps: old has 0 nonzero off-diagonals; new has 30, all
  inside rows/cols 3–8; new's `V[3..8,3..8]` equals `e(V_es_att)` — upper triangle bitwise,
  diagonal kept from `e(V)` (≤ 1.7e-16 relative from V_es, the "same quantity computed twice"
  of the commit message), lower triangle mirrored so the posted matrix is exactly symmetric
  where margins' raw `r(V)` is not (worst 8.3e-16); e.g. new `V[ES_att_pre3,ES_att_pre2]` =
  `+1.56756478d1fdeX-00a` (.00130637574894) where old holds `+0.0000000000000X-3ff`. My old and
  new `V` equal the suite's `run_old`/`run_new` `basic` dumps bitwise.
- **Classification: intended-change.** Author-ruled design change; the diagonal — every
  reported standard error — is unchanged to the bit on every config; no new estimation.

### I.2 `e(cmdline)` begins with the command name (CORE-10) — INTENDED CHANGE; at `HEAD` the comparer normalises the prefix and no row is needed

> **What this item is judged against (framing rewritten 2026-09-14).** At `HEAD` (`4280f3b`) there
> is **no ledger row for this change, and none is needed**: the comparer strips the uniform
> `"wooldid "` prefix before comparing, so the difference is not reported at all and the *rest* of
> `e(cmdline)` is compared instead of being blanket-blessed by an unbounded all-config row. The code
> is `test/suite/compare_oldnew.py` **lines 157-158** — `if p[1] == "cmdline" and
> v.startswith("wooldid "):` / `v = v[len("wooldid "):]` — under a comment that states the same
> 399 this item measures (*"new
> e(cmdline) = 'wooldid ' + the v1.3 text (equal on all 399 rc-0/0 configs, 2026-09-12)"*). An
> earlier draft was framed on the blanket row `CORE10-cmdline`, retired in the same commit, together
> with `B21b-ic8-p`, which the unbounded CORE-10 row had been hiding (`B21b-ic8-s` and `-m` survive).
> **The change to `e(cmdline)` and every measurement below are unaffected; only the framing moved.**
> `legacy` L.1 and `si20` S.3 also cite the retired id; both already annotate it inline as *"retired
> at `4280f3b`"*, and neither was rewritten in this pass. See `LEDGER_DELTA.md` §A and §G(3).

- **Symptom / impact:** old's `e(cmdline)` omitted the leading `wooldid `, contrary to Stata
  convention (every official estimator's `e(cmdline)` is re-runnable as typed).
- **Mechanism (v1.3):** old line 2995 `` ereturn local cmdline `"`0'"' ``; new line 4035
  `` ereturn local cmdline `"wooldid `0'"' ``.
- **Fix + proof:** `aec904e` (2026-08-19, "isolated re-bless", verified row-by-row —
  `wooldid_solved_issues.md` batch record at line 402 ff.). Uniform 8-byte prefix, no other
  content change; this is why MATCH has been 0 since the 08-23 run.
- **This run:** **all 399 rc-0 configs carry the `"wooldid "` prefix on `e(cmdline)` and nothing
  else in that macro** — that is the standing measurement, and it is the same number the comparer's
  own CORE-10 comment records. At `26ce421` the bookkeeping split it: `CORE10-cmdline` fired on
  **398** configs, and the 399th, `cov_icbad2`, carried the identical prefix but had its cmdline
  diff claimed by the earlier-ordered config-scoped wildcard `B21b-ic8-p` (first match wins). That
  398/1 split was an artifact of which row claimed the cell, never of the builds; at `4280f3b`
  neither row exists and the cell is not reported at all.
- **Reproduction (`repro_intended.do`, block `basic`):**
  ```
  REPRO_CORE10 build=old rc=0 e(cmdline)=[y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats jointtests safety(off)]
  REPRO_CORE10 build=new rc=0 e(cmdline)=[wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats jointtests safety(off)]
  ```
- **Classification: intended-change.**

### I.3 Stable cfx row keys (NFFAST-02, ledger `A6-cfxstripe`) — INTENDED CHANGE

- **Symptom / impact:** old's `e(cfx*)` matrices carried row names taken from margins'
  temporary factor variables — `1.__00000G:1bn.__00000H`-style stripes that change from run
  to run and mean nothing to a consumer (a script keying on row names could not match the
  same lattice bin across two runs).
- **Mechanism (v1.3):** the cfx matrices are transposed slices of margins' `r(table)` — old
  lines 1913, 1937 (main) and 1972, 1993 (per-subgroup), e.g.
  `` matrix `cfx`item'_`nametext'' = r(table)[1..`cfxmatendpoint', ...]' `` — so their row
  names inherit margins' column labels, which are built from the tempvar names of the lattice
  and target dummies. New renames them to `lattice_#` with null roweq (new line 10221,
  `` local newnames `newnames' lattice_`=regexs(1)' ``).
- **Fix + proof:** `cb1996c` (2026-08-19, the final issue-file batch incl. NFFAST-02;
  `wooldid_solved_issues.md` line 402 ff.). Matnames-only: the numeric cfx cells are covered
  by their own rows (A5-cfx-n for the B7 width restoration, SI20-P-* and X-est-cfx for values).
- **This run:** `A6-cfxstripe` fired on **7 configs / 30 matnames lines** (`trig_cfx`,
  `plot_cont`, `ext_char_cfx_nonatt`, `ext_cont_ext_cfx`, `ext_cont_ext_cfx_subgroups`,
  `ext_cov10_cfx_inference`, `si20_lat_cfx_pw`; the seventh is new since 08-23). The
  volatility is directly observed: on the SAME config and data, the suite's `run_old` dump and
  my old run carry different stripes (`__000009` vs `__00000H`; the verbose prints them in the
  comparer's normalized `__TMP1` spelling) — so my old `trig_cfx` dump is the only one of the
  seven that does not equal `run_old` at rtol 0, and only in those names.
- **Reproduction (`repro_intended.do`, block `trig_cfx`: `att att_it contreat(intensity)
  contreatcontrols(basic) lattice(all 3) makecfxplots cfxplottypes(post 2) safety(off)
  cluster(i) esprelength(2) espostlength(2)` on `dta_cont.dta`):**
  ```
  REPRO_A6 build=old rc=0 e(cfxatt_cfxpost) 3x6 rownames=[1.__00000G:1bn.__00000H 1.__00000G:2.__00000H 1.__00000G:3.__00000H] colnames=[b se t pvalue ll ul]
  REPRO_A6 build=old e(cfxatt_it_cfxpost) 3x1 rownames=[1.__00000G:1bn.__00000H 1.__00000G:2.__00000H 1.__00000G:3.__00000H] colnames=[b]
  REPRO_A6 build=new rc=0 e(cfxatt_cfxpost) 3x6 rownames=[lattice_1 lattice_2 lattice_3] colnames=[b se t pvalue ll ul]
  REPRO_A6 build=new e(cfxatt_it_cfxpost) 3x6 rownames=[lattice_1 lattice_2 lattice_3] colnames=[b se t pvalue ll ul]
  ```
  (The 3x1 → 3x6 width change on the second matrix is the B7 fix, §1.3 of the old report, not
  this row.)
- **Classification: intended-change** (a naming contract, author-approved; the old names were
  unusable rather than numerically wrong).

### I.4 `_ppmlhdfe_d` left in the user's dataset on every clean hdfe-poisson run (B12 + S6, ledger `E-B12`, `E-B12e`) — OLD BUG, FIXED

- **Symptom / impact:** every clean `poisson` run on v1.3 left ppmlhdfe's internal
  `_ppmlhdfe_d` variable in the user's data (and the same inverted condition actively
  deleted a user's own pre-existing variable of that name on the plain-poisson path).
- **Mechanism (v1.3):** old line 3440
  `` if "`poisson'" == "poisson" & ("`unconditionalse'" != "" | "`customvce'" != "" ) qui cap drop _ppmlhdfe_d ``
  — but ppmlhdfe's `d` option runs only when `unconditionalse`/`customvce` are EMPTY, so the
  drop was a no-op on the path that creates the variable and fired on the path that never does.
- **Fix + proof:** `9d9d26f` (2026-07-20; condition mirrors the creation gate; verified 3+3
  incl. the rc-3200 fallback — `wooldid_solved_issues.md` B12) plus S6 `edebab1` (2026-07-27;
  cleanup obligation registered at creation, epilogue drop frame-qualified; two verifiers,
  decoy variable survives on the fixed build and dies on the parent — tranche-13 record).
  Under `verbose` the variable is deliberately retained (author directive).
- **This run:** `E-B12` (varlist) and `E-B12e` (e(sample) digest) each fired on **8 configs**:
  `pois`, `trig_tt_pois`, `adm_pois_pw`, `ext_contract10_poisson_converged`,
  `ext_contract10_ppml_d_collision`, `ext_path10_ppml_separation`, `si20_ic_pois`,
  `si20_lat_pois` (the last two new since 08-23). E.g. `pois`: "varlist var REMOVED in new:
  _ppmlhdfe_d" and "esamp _ppmlhdfe_d: REMOVED in new (was -63.779850742850115)".
- **Reproduction (`repro_intended.do`, block `pois`: `ycount i t ttre, att cluster(i) poisson`):**
  ```
  REPRO_E_B12 build=old rc=0 confirm-variable _ppmlhdfe_d after the run: rc=0   (0 = LEFT in the user's data; 111 = absent)
  REPRO_E_B12 build=new rc=0 confirm-variable _ppmlhdfe_d after the run: rc=111   (0 = LEFT in the user's data; 111 = absent)
  ```
- **Classification: old-bug-fixed** (data-hygiene defect; estimation numerics untouched —
  `pois` shows no numeric row).

### I.5 saveplots x subgroups filename scheme (ledger `C-esfn`) — INTENDED CHANGE

- **Symptom / impact:** the per-subgroup event-study data file was named
  `<stub>_ES_<var><value>.dta`, concatenating the subgroup variable and its value with no
  separator — ambiguous for multi-digit values (`region1` + value 2 vs `region` + value 12).
- **Mechanism (v1.3):** old line 3404
  `` if "`saveplots'" != ""  qui save "`saveplots'_ES_`subgroups'`sgcase'.dta", `replace' ``;
  new line 12950 `` qui save "`saveplots'_ES_`subgroups'_`sgcase'.dta", `replace' ``.
- **Fix + proof:** `06e1a5b` (2026-07-26, "saveplots: separate the subgroup variable and value
  in the ES data filename"). The `.gph` names (`_ES_att_sg#.gph`) are unchanged.
- **This run:** `C-esfn` on `adm_svpsg` (the `files|` dump line): old
  `svpadm_ES_region1.dta … region4.dta` → new `svpadm_ES_region_1.dta … region_4.dta`; the
  other six files identical.
- **Reproduction (`repro_svpsg.do`: `y fips year ttre, att cluster(state) subgroups(region)
  esprelength(1) espostlength(1) makeplots saveplots(svpadm)` on the suite's admin panel with
  its `adm_sub` unit-i subsample; the first attempt's `di` died on the embedded quotes of the
  file list, so the file is printed one name per line):**
  ```
  REPRO_C build=old rc=0 saveplots wrote the following files (one per line):
  REPRO_C_FILE build=old svpadm_ES.dta
  REPRO_C_FILE build=old svpadm_ES_att.gph
  REPRO_C_FILE build=old svpadm_ES_att_sg1.gph ... svpadm_ES_att_sg4.gph
  REPRO_C_FILE build=old svpadm_ES_region1.dta
  REPRO_C_FILE build=old svpadm_ES_region2.dta
  REPRO_C_FILE build=old svpadm_ES_region3.dta
  REPRO_C_FILE build=old svpadm_ES_region4.dta
  REPRO_C build=new rc=0 saveplots wrote the following files (one per line):
  REPRO_C_FILE build=new svpadm_ES.dta
  REPRO_C_FILE build=new svpadm_ES_att.gph
  REPRO_C_FILE build=new svpadm_ES_att_sg1.gph ... svpadm_ES_att_sg4.gph
  REPRO_C_FILE build=new svpadm_ES_region_1.dta
  REPRO_C_FILE build=new svpadm_ES_region_2.dta
  REPRO_C_FILE build=new svpadm_ES_region_3.dta
  REPRO_C_FILE build=new svpadm_ES_region_4.dta
  ```
- **Classification: intended-change.**

### I.6 Interactive-control copies stored as float: silent garbage at rc 0 for controls with levels above 2^24 (SI-15 / B21(b), ledger `B21b-ic8-s/m/p`, config `cov_icbad2`) — OLD BUG, FIXED

- **Symptom / impact:** with `interactivecontrols(x)` and a double-typed `x` whose LEVEL is
  large relative to its within-unit variation, v1.3 returned an absurd ATT with an inflated
  variance at rc 0, e(V) posted. For administrative-style controls (identifiers, dollar
  amounts, dense date codes at the 1e8 scale) a user would have published a sign-flipped
  estimate with a 31x-too-large standard error and no message.
- **Mechanism (v1.3):** old line 923 `` qui gen `ic`iter'' = `intercon' `` copies the control
  as FLOAT. Float is exact only to 2^24 = 16,777,216; at ~1.02e8 (between 2^26 and 2^27) the
  float spacing is 8, so each copied value is off by up to 4 (measured below). The
  interactive-control machinery is engineered so that the cell-slope columns are EXACTLY
  collinear with the absorbed terms and drop cleanly; with the quantized copy the collinearity
  is only approximate, the drop does not happen, and reghdfe fits an ill-conditioned model.
  New line 1132 `` qui gen double `ic`iter'' = `intercon' ``.
- **Fix + proof:** `3cb0c41` (B21, 2026-07-21; seven `gen` → `gen double` conversions:
  ringer copy + powers, IC copies, four w_att* weights). The B21 assessment found exactly this
  class at the IC site — "rc 0 with absurd coefficients ... because the float copy's ~1e-5
  quantization noise defeated the engineered exact-collinearity drop"; two verifiers reproduced
  it independently, one proving by single-site variant builds that the ic edit alone is the
  entire cure and the trigger boundary needs level-to-within-spread ratios in the hundreds
  (`wooldid_solved_issues.md` B21). The suite pin `cov_icbad2` (SI-15, `8f7e830`, 2026-08-25)
  rescales `icbad` by 100 to ~1.02e8 precisely to keep this class visible.
- **This run (`cov_icbad2`, ACCEPTED BLIND SPOT rows, unbounded):** `B21b-ic8-m` — `b[1,1]`
  −1.792004230359348 → 1.300034583785127 (reldiff **1.73e+00**), `V[1,1]` 3.509305756440952 →
  0.0036406763974571 (9.99e-01), `wdidmainresults` 6 cells (se 1.8733141104579745 →
  0.06033801784494699; t −0.9565957040281149 → 21.545861634465318); `B21b-ic8-s` — seven fit
  scalars (ll −975.9488635332602 → −979.6667135423407, r2 .97164487 → .97152750, rss
  311.70352748 → 312.99382067, …); `B21b-ic8-p` — fired on `cmdline` only.
- **Reproduction (`repro_intended.do`, blocks `cov_icbad` and `cov_icbad2`: `y i t ttre, att
  cluster(i) interactivecontrols(icbad)` on `dta_s_covars.dta`, where `icbad = 1e6 + 137*i + t`;
  cov_icbad2 = the same after `replace icbad = icbad*100`). A three-value triangle proves the
  old value WRONG, not merely different: rescaling one control by a constant cannot move the
  ATT; at ~1e6 the float copy is exact and both builds agree to the bit; at ~1e8 new preserves
  that value to 1e-15 while old flips its sign:**
  ```
  REPRO_B21b_ic8 build=old cov_icbad  level~1e6: max|icbad-float(icbad)| = 0
  REPRO_B21b_ic8 build=old cov_icbad  rc=0 att b[1,1]=+1.4ccf110506ad8X+000 =  1.300034583785125  se=.0603380178449468  properties=[b V]
  REPRO_B21b_ic8 build=old cov_icbad2 level~1e8: max|icbad-float(icbad)| = 4
  REPRO_B21b_ic8 build=old cov_icbad2 rc=0 att b[1,1]=-1.cac0ca0bafe28X+000 = -1.792004230359348  se=1.873314110457974  properties=[b V]
  REPRO_B21b_ic8 build=new cov_icbad  rc=0 att b[1,1]=+1.4ccf110506ad8X+000 =  1.300034583785125  se=.0603380178449468  properties=[b V]
  REPRO_B21b_ic8 build=new cov_icbad2 level~1e8: max|icbad-float(icbad)| = 4
  REPRO_B21b_ic8 build=new cov_icbad2 rc=0 att b[1,1]=+1.4ccf110506ae0X+000 =  1.300034583785127  se=.060338017844947  properties=[b V]
  ```
- **Classification: old-bug-fixed** (silent wrong answer at rc 0 on realistic-scale input).
- **Record disagreements (open questions):** (a) the ledger row `B21b-ic8-p` hedges "old
  garbage path may not post e(V) (properties)"; on this run and in my reproduction old DID post
  `[b V]`, and the row fired on `cmdline` only — the hedge is stale text, harmless because the
  row is a config-scoped wildcard. (b) the B21 solved-issues record describes the garbage mode
  as "rc 0 with absurd coefficients, no e(V), missing SEs throughout wdidmainresults";
  `cov_icbad2` exhibits a different symptom of the same class — garbage b WITH a posted garbage
  V (se 1.87 in place of 0.060). Same defect and same single-line cure; the two records describe
  two failure shapes, and neither is wrong about its own case.

### I.7 Float-precision storage defects: `gen` float truncation (B21) and `svmat` float truncation (`6267522`) — the bounded storage-noise tiers, magnitudes UPDATED for this run — OLD BUGS, FIXED (small magnitudes)

- **Symptom / impact:** old stored the engineered-collinearity (`ringer`) copies, the
  interactive-control copies, the four estimand weights and every `svmat`-materialized
  intermediate (Ibragimov–Müller per-cluster coefficients, adversarial fit inputs, plot frames)
  as FLOAT. On ordinary data the consequence is last-digit noise in the affected estimands
  (att_i/att_itime cell weights, IM b/se/t/p, adversarial p and line), tail-amplified in p-values
  and joint-test scalars; the one silent-garbage mode of the same class is §I.6.
- **Mechanism (v1.3):** `gen` (float) at old lines 877 and 884 (`ringer` and its powers), 923
  (interactive controls, §I.6), 1080/1087/1094/1100 (`w_att_i`, `w_att_itime`, `w_att_it`,
  `w_att`); `svmat` (float by default; ignores `set type double`) at the 16 old sites 1460,
  1491, 1563, 1594, 1652, 1746, 1820, 1849, 2036, 2124, 2203, 2244, 2343, 2535, 3232, 3332. New:
  `gen double` at 1086/1093/1132/1549/1556/1563/1569; every `svmat` in the new build is
  `svmat double` (11 occurrences; the 12th grep hit is a comment) and the four adversarial
  `predict` calls are `predict double`.
- **Fix + proof:** `3cb0c41` (+ docs `7672719`; battery 16 configs byte-identical, allest_aw
  and sg_attitime shifted at storage-noise magnitude ≤ ~1.2e-8 relative on estimates,
  tail-amplified p to ~5e-8 and one joint-test scalar at 1.8e-29 shifting 4.4e-6 relative; 3/3
  verifiers) and `6267522` (occurrence-count-asserted edits; 16/18 battery configs bit-identical,
  diffs confined to impv and esf_adv; typical ~4e-8 relative on estimates, tail-amplified on p;
  the IM variance shift unbounded only in degenerate corners; 3/3). New values are the correct
  ones; movement is location-confined.
- **This run — every bounded tier fired inside its bound** (the "worst fired reldiff" lines of
  the verbose tail; bounds are ~10x the calibrated ceilings). Changes from the 08-23 report
  are marked:
  | tier | configs | worst reldiff (bound) | note |
  |---|---|---|---|
  | `A8-B21-b/V/w` | `allest_aw` | 1.22e-08 / 4.73e-09 / 5.18e-07 (1e-6) | main_att bit-identical; only att_i/att_itime cells move |
  | `A8-sgit-b/V/w/s` | `sg_attitime` | 9.85e-09 / 3.51e-09 / 5.14e-07 / 5.78e-07 (1e-6) | |
  | `A8-sgit-p`, `A8-sgit-j` | `sg_attitime` | 3.30e-06 / 3.30e-06 (1e-4) | pwsgs p cells and `joint_att_itime_allsgssame` (1.79021054e-29 → 1.79021646e-29) tail-amplified |
  | `A9-svmat-b/V/w` | `impv` | 4.14e-08 / 2.39e-08 / 2.12e-07 (1e-6) | unchanged from 08-23 |
  | `A9-adm-b/V/w` | `adm_impv` | 7.35e-09 / 7.21e-07 (1e-5) / 1.66e-05 (1e-4) | IM t/p tail-amplified; unchanged |
  | `B21-cwm-b/w` | `cov_cwm` | 4.43e-08 / 8.71e-08 (1e-6) | |
  | `B21-adm-b/V/w` | `adm_aw` | 7.47e-10 / 2.06e-09 / 3.13e-08 (1e-6) | |
  | `XIC-f` | `ext_*interactive*` | 4.65e-08 (1e-3) | **CHANGED**: the 08-23 headline "IC-copy 1.33e-04" is no longer attributable to B21 — the ledger retired `XIC-b/V/w` into the `SI20-I-*` rows (ledger comment at line 128; SI-20 `1df8613` is a model-construction change outside this group). Only the fit-scalar row remains. |
  | `X-est-b/w` | 105 `ext_*` configs | 9.46e-05 / 9.47e-05 (2e-4) on `ext_c1_h22_p` | elasticity-transform amplification; 106 → 105 configs |
  | `X-est-V` | `ext_*` | 2.98e-07 (1e-5) | |
  | `X-est-sg` | `ext_*` | 4.19e-05 (4e-4) on `ext_c2_h22_f2` | |
  | `X-est-pw` | `ext_*` | 1.12e-04 (1e-3) on `ext_c2_h00_f2`/`h03_f2` | degenerate-SE pwsgs cells, documented run-varying |
  | `X-est-hist/cfx` | `ext_*` | 6.73e-09 (1e-7) / 3.86e-10 (1e-7) | |
  | `X-joint-s` | `ext_*` | 4.09e-06 (5e-5) | |
  | `X-fit-F/rss/rmse/ll0/ll` | `ext_*` | 2.12e-07 / 5.94e-08 / 2.97e-08 / 5.16e-08 / 1.59e-08 (2e-6) | |
  | `X-fit-r2/chi2` | `ext_*` | 3.24e-09 / 1.91e-09 (1e-7) | |
  | `ULP-ppw-b/V/w` | `adm_pois_pw` | 1.16e-10 / 5.56e-10 / 3.18e-09 (1e-6) | ppmlhdfe+pweight IRLS, run-varying on both builds |
  | `ULP-upw-b/V/w` | `adm_uncond_pw` | 6.01e-10 / 1.04e-08 (1e-6) / 2.92e-07 (1e-5) | unconditionalse 1–2 ulp class |
  | `A12-ttp-s` | `trig_tt_pois` | 2.06e-09 (1e-4) | tcen `f7019c9` poisson; documented run-varying (1.24e-5 on another run) |
  | (`A16-imd-b`, `A17-part-b`) | `trig_imdegen`, `trig_imdeg_part` | 3.01e-08 (1e-6), 2.41e-06 (1e-4) | svmat-class members owned by the IM-degeneracy sections, listed for completeness |
- **Reproduction (`repro_intended.do`, blocks `allest_aw` and `impv`):**
  ```
  REPRO_A8_B21 build=old b[1,1] = +1.73a2f01a5badcX-001  ( .7258524925399246)
  REPRO_A8_B21 build=new b[1,1] = +1.73a2f01a5badcX-001  ( .7258524925399246)      <- main_att bit-identical
  REPRO_A8_B21 build=old b[1,3] = +1.4a37990572df4X-001  ( .6449554270888584)      <- main_att_i
  REPRO_A8_B21 build=new b[1,3] = +1.4a37990f48c8cX-001  ( .6449554282338723)
  REPRO_A8_B21 build=old b[1,4] = +1.3dedc7e387394X-001  (  .620954748653832)      <- main_att_itime
  REPRO_A8_B21 build=new b[1,4] = +1.3dedc7a24f250X-001  ( .6209547410613343)
  REPRO_A9_svmat build=old rc=0 b[1,1] = +1.7d797a0f9999aX-001 ( .7450674194376916)   V[1,1] = +1.11bc330cea630X-008 ( .0041768669733897)   e(ll)=.6118055
  REPRO_A9_svmat build=new rc=0 b[1,1] = +1.7d797a08b0d78X-001 ( .7450674186333535)   V[1,1] = +1.11bc3321859f5X-008 ( .0041768669921311)   e(ll)=-400.91012
  ```
  (The `e(ll)` change on `impv` is B4, §1.2 of the old report.)
- **Classification: old-bug-fixed** — genuine float-storage defects with small, bounded,
  location-confined magnitudes on well-behaved data; the unbounded member of the class is §I.6.

### I.8 The rc-by-design class: 66 designed return-code changes, grouped — MIXED (see per-group classification)

**Reconciliation.** `configs.csv`/`configs_ext.csv` carry 67 rows with `rc_old != rc_new`; one
(`ext_contract10_nonpositive_time`, rc_old 0, rc_new BLANK) is an unset pin and was judged
EXPLAINED (rc 0 on both builds, cmdline only), leaving the **66** the verbose counts. Against
the 08-23 report's 32 designed pins: `ext_contract10_rcs_factor_fe` (0 → 198, the §2.2
naming-fragility defect) LEFT the class — SI-7 was fixed 2026-08-27 and the config is now
0/0 EXPLAINED; 31 carried forward (their attributions in §1.16/§3.4 of the 08-23 report
stand, the `b42b71f` trio adversarially verified there) and **35 arrived**: 33 new-only
options plus `cov_ccc_absorb` and `ext_contract10_poisexp_without_poisson`. Old-side r198 on a
new-only option is v1.3's unknown-option error, correct old behaviour.

**Reproduction (`repro_intended.do`, block `rc-by-design probes`, one probe per bucket):**
```
REPRO_RC build=old scr_fastnose (new-only option) rc=198          REPRO_RC build=new scr_fastnose (new-only option) rc=0
REPRO_RC build=old cont_uz_ss (crash cure B5) rc=198               REPRO_RC build=new cont_uz_ss (crash cure B5) rc=0
REPRO_RC build=old scr_lat0 (informative refusal d1f7542) rc=0     REPRO_RC build=new scr_lat0 (informative refusal d1f7542) rc=198
REPRO_RC build=old scr_impsg (both-error parity) rc=198            REPRO_RC build=new scr_impsg (both-error parity) rc=198
```

#### I.8a New-only options: old r198 unknown-option, new rc 0 (40 configs) — INTENDED CHANGE

| configs | option(s) | commit(s) |
|---|---|---|
| `scr_fastnose` | `fastnose` | F9 `fc0cadb` (2026-07-26) |
| `scr_antic` | `anticipation()` | F4 `1f521f6` (2026-07-26) |
| `scr_ptp`, `patho_ptp` | `pretrendpower()` | F7 `f9fd239` (2026-07-26) |
| `gt_basic`, `patho_gt` | `gteffects` | F3 `b65b91b` (2026-07-27) |
| `perm_basic`, `patho_perm` | `permute()` | F6 `125bbd0` + engine `86d37cb` (2026-07-27) |
| `scr_permexact_t` | `permexact permtstat permtci` | F6 `86d37cb`, `b544c55` |
| `scr_permstrata` | `permstrata()` | F6 `86d37cb` |
| `perm_absp` | `permabsp` | `310fa61` (2026-09-05) |
| `perm_ase`, `perm_save` | `permanalyticse`, `permsave()` | permute completeness chunk 1, `1561386` (2026-09-05) |
| `perm_side_left`, `perm_rawp` | `permside()`, `permrawp` | `dc3f345` (2026-09-06, F064/F052) |
| `perm_stepdown` | `permstepdown` (opt-in) | PS-45 step-down half `7a158ea` (2026-09-10) |
| `perm_cfxnose` | `permute()` + `cfxplotnose` | D-1 `82934b3` (2026-09-11) on the F6 base |
| `boot_std_att`, `boot_std_reid`, `boot_obs`, `patho_boot` | `bootstrap(standard)`, `reps()`, `rseed()`, `dots()` | bootstrap option surface `8cfd915` (2026-08-30) |
| `boot_bayes` | `bootstrap(bayes)`, `bootalpha()`, `level()` | `1ad2f65` (Bayesian engine, 2026-08-30) + `c9b47d8` (F001, 2026-09-06) |
| `boot_jt`, `boot_sgjt`, `boot_aw_ties` | bootstrapped `jointtests`, subgroup contrasts, `ties` | `8d8316f` (2026-08-30) |
| `boot_tdf` | `bootdf()` | `4906bd3` (2026-08-30) |
| `boot_disp`, `boot_save`, `boot_tstat` | `bootdisplay()`, `bootsave()`, `boottstat` | Layer 2 `1561386` (2026-09-05) |
| `boot_bca` | `bootdisplay(bca)` | F067 `366ec87` (2026-09-06) |
| `boot_ts_gt` | `boottstat` + `gteffects` | F033 `2f24ce9` (2026-09-10) + cells `82934b3` |
| `band_analytic`, `band_boot` | `simultaneousbands` | PS-45 bands half `8ea5a25` (2026-09-10) |
| `band_multi` | `simultaneousbands` + two families | `8ea5a25` + D-9 `82934b3` (2026-09-11) |
| `fcol_base`, `fcol_fn`, `fcol_pois`, `fcol_sg` | `fastcollapse` | `1eea605` (2026-09-01), subgroups `4b224f0` (2026-09-02) |
| `fcol_slow` | `fullprecisionslow` | `111ecbb` (2026-09-01) |
| `level_analytic` | `level()` on the analytic path | F054 `a969126` (2026-09-06) |

Seven of these (`gt_basic`, `perm_basic`, `scr_antic`, `scr_fastnose`, `scr_ptp`,
`scr_permexact_t`, `scr_permstrata`) are the 08-23 §3.4 pins; the other 33 are new since.
Classification: **intended-change** (old's r198 is correct old behaviour, not a defect).

#### I.8b Crash cures: old aborted, new completes correctly (7 configs) — OLD BUGS, FIXED

All seven carried from the 08-23 report's §1.16 bucket A; attributions unchanged:

| config | rc old → new | fixed defect |
|---|---|---|
| `cont_uz_ss` | 198 → 0 | B5 `4979abf`: contreat + untreatedZero + summarystats crashed r198 ereturning an undefined local after estimation completed |
| `adv_gapa` | 198 → 0 | B16 `c0e17f8`: adversarial reference-gap hunt aborted on gap-year data |
| `adv_r2001` | 2001 → 0 | **B16 `c0e17f8`** (the cut fix — the wrong cut had selected a single row): **at `c0e17f8` the config is already rc 0 with a posted `adversarial_pval_att` of 9.47227267238e-21**, one commit before `05a243b`, whose guards cover only the residual `min_t == 1` and gap-hollowed cells. *Corrected 2026-09-13:* this row read "B16 `05a243b`: degenerate 1-row pre-trend fit aborted; now clean skip + Caution", copied from the then-stale `configs.csv` note. Both halves are false — the new build **posts p 1.20615084072e-05 and a line for the config, with no skip and no Caution** (`adversarial` S1 and OQ8; `git show c0e17f8:test/baseline/adv_r2001.csv` gives `rc\|0` and `scalar\|adversarial_pval_att\|+1.65da3335557fcX-043`, and `git merge-base --is-ancestor c0e17f8 05a243b` succeeds). `configs.csv` was corrected at `4280f3b`. **`adversarial`'s account is the one of record.** |
| `trig_sgadv_deeppre` | 2000 → 0 | `4483758`: SG adversarial skiplocation crash |
| `trig_bigcodes` | 461 → 0 | ttre double storage `5a97f70`: YYYYMMDD onset grids aborted r461 (the silent cohort-merging variant proven in the `8c5bfd7` work) |
| `scr_badif` | 411 → 0 | B23 `8d7bd01`: input asserts ran over the whole dataset — false r411 on rows the user excluded via `if` |
| `trig_esfsup_ref322` | 322 → 0 | G10B `7ba1e16`: under default safety the unsupported-reference cohort made old abort opaquely; new drops it with a Caution |

#### I.8c Informative refusals: old ran silently wrong or failed cryptically/late, new refuses with a message (18 configs) — OLD BUGS, FIXED

Seventeen carried from the 08-23 §1.16 bucket B (attributions unchanged), plus one arrival:

| config | rc old → new | fixed defect |
|---|---|---|
| `scr_uncond_alone` | 322 → 198 | B15 `14ca983`+`657f6f0`: unconditionalse without robust/cluster wasted full estimation then bare r322 |
| `scr_cfxvar` | 0 → 198 | B17 `e4e3d1e`: a cfxplottypes token naming an existing variable was silently injected — rc 0 with a meaningless plot |
| `scr_nonintsg`, `ext_contract_fractional_subgroup` | 198 → 411 | B20 `7f7f73c`: non-integer subgroups died cryptically (`_joint_horiz_i_s1.5 invalid name`) |
| `scr_lat0` | 0 → 198 | latticenum `d1f7542`: old silently remapped `lattice(all 0)` to the gquantiles default |
| `scr_latnocon` | 0 → 198 | `4c8fb66`: `lattice()` without `contreat()` parsed then wholly ignored |
| `scr_fracttre`, `ext_contract10_fractional_time` | 0 → 459 | S5 commensurability screen (`905617d`..`f187d7d`): in-window fractional onsets ran rc 0 with the onset period mislabeled |
| `ext_val10_fractional_ttre` | 452 → 459 | same screen: old failed late/opaque |
| `ext_contract10_inconsistent_ttre`, `trig_mintsent_def` | 0 → 459 | F13 `c4d7e62`: ttre varying within `i` ran silently (+41% horizon bias measured) |
| `scr_sgsplit`, `adm_sgsplit` | 0 → 459 | F13 `c4d7e62` + `e5445c4`: units split across subgroups made SG contrasts meaningless (old returned the pooled estimate under every label) |
| `ext_contract10_if_outside_negative` | 411 → 322 | B23 `8d7bd01`: old's raw-domain assert false-fired on excluded rows; new runs the restricted sample, which fails r322 for the genuine reason |
| `ext_contract_bad_elasticity`, `ext_contract10_bad_safety`, `ext_contract10_poly_zero` | 0 → 198 | `b42b71f` VAL-G19/G21: silent-accept option domains closed (old rewrote an unknown elasticity type to `dydx`, ran `contreatpoly(0)` as order 1, treated garbage `safety()` as on) |
| **`ext_contract10_poisexp_without_poisson`** (arrival) | 0 → 198 | SI-9 `083ea89` (2026-08-25, author-approved): `poisexpresults` without `poisson` silently accepted on old; now an option-contract r198 (cert `test_si89`) |

#### I.8d `cov_ccc_absorb` (198 → 459) — BOTH-BUILDS DEFECT, rc moved by a message repair only

`ccc_absorb(x2)` absorbs the cohort-by-period treatment interactions themselves, so the
fitted model carries no treatment coefficient at any horizon. Old died with a bare
invalid-syntax r198 in the margins phase; since the BS-33 guard, `5c78560` (2026-09-07)
separates "no such period in the sample" from "the fitted model carries no coefficient" and
in the model case errors r459 pointing at absorbed fixed effects and collinear control sets.
The commit body is explicit that "the underlying defect, that ccc_absorb is unusable for event
studies, is unchanged and still open" (PARITY DEFECT TICKET in `configs.csv`; REPORT_2026-08-18
ticket 2). Classification: **both-builds-defect** — the rc change is a designed message
improvement on top of a fault neither build fixes.

#### I.8e Both-error parity pins (42 configs, RC_EXPLAINED as parity, not rc-by-design)

Identical rc on both builds, asserted by pin; the comparer checks the code only (message
texts are out of contract). Grouped by nature:
- **Input/option-contract refusals present on both builds (36):** the 32 `ext_val_*` /
  `ext_val10_*` / `ext_contract_fractional_lattice` / `ext_cov10_error_emptycells_restore`
  validation cases (r198/411/452 — too few/many vars, poisson+aweight, impvalues combos,
  lattice/contreat grammar, negative pre/post/weights, non-positive ids, fractional lattice)
  and `scr_cfxbad`, `scr_makecfx_nocon`, `scr_impsg`, `fcol_gate` (r198; `fcol_gate` is the
  fastcollapse-with-interactivecontrols gate on new vs unknown option on old — same code by
  different routes).
- **Design rejections (4, r459):** `patho_nonever`, `ext_val10_zero_coded_never` (no
  never-treated group; author ruling `cb199c8`), `patho_nocontrols`, `ext_val10_mismatched_clock`.
- **Both-build failure shapes (2), ticketed elsewhere:** `trig_pois_sep` (r322, poisson
  separation; A7 coverage lives in `trig_allest_drop`/`cov_*`) and `scr_strclus` (r109, string
  cluster variable — old a bare late type mismatch, new clean at parse; verified 2026-08-18).
Not a divergence class; listed because the task asked for the full rc taxonomy.

---

### Open questions (collected)

1. **PS-62 config count.** Ledger row `PS62-esblock` says "18554 cells over 360 configs on the
   2026-09-10 dual-build run (main 301b706 vs 8ea5a25)". This run: 18,554 cells — identical —
   over **335** configs. Either the 360 was a miscount on 09-10 or the 09-10 run had 25 more
   ES-carrying configs that contributed 0 cells; the cell identity argues for the former. The
   row is a wildcard and fires correctly either way; the text should be corrected to 335.
   **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E6]** — `4280f3b` corrected it: the reason now reads *"18554 cells over
   335 configs on the 2026-09-11 dual-build run (main 301b706 vs 14c30d3)"*.
2. **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E7]** (**moot — the row was removed**: `4280f3b` deleted
   `B21b-ic8-p` along with `CORE10-cmdline`, the unbounded row that had been hiding it; the
   comparer now normalises the `"wooldid "` prefix instead.) The question as written was:
   **`B21b-ic8-p` hedge is stale.** Old posts `[b V]` on `cov_icbad2` (properties identical on
   both builds); the row fires on `cmdline` only. Harmless (config-scoped wildcard) but the
   "may not post e(V)" text does not describe what old does on this config.
3. **Two failure shapes for one B21(b) class.** The solved-issues B21 record describes the IC
   garbage mode as "no e(V), missing SEs"; `cov_icbad2` shows garbage b with a posted garbage V.
   Same site, same cure, different symptom; worth a one-line addendum to the record.
4. **[CLOSED by `4280f3b`, 2026-09-13 — retained, not deleted; see `LEDGER_DELTA.md` §E8]** (**moot — `CORE10-cmdline` was removed at `4280f3b`** and the
   difference is no longer reported at all.) The note as written was:
   **CORE10 count.** 398 firings + `cov_icbad2` (claimed by `B21b-ic8-p`) = all 399 rc-0
   configs; no config lacks the prefix. Cosmetic only — noted so the 398 is not read as a gap.
5. **XIC retirement.** The 08-23 report's §1.14 headline magnitude (1.33e-04 IC-copy) is no
   longer a B21 statement: those cells are now owned by the SI-20 rows, which are a
   model-construction change and belong to another group's sections.
6. **`ext_contract10_nonpositive_time`** carries a blank `rc_new`; it is judged EXPLAINED
   (rc 0 both) so nothing is wrong, but an explicit `0` would make the pin auditable.
7. **A6 volatility evidence** rests on two old runs of `trig_cfx` producing different stripes
   (`run_old` vs my run); the verbose shows the comparer's normalized `__TMP1` spelling, so
   the raw dumps, not the verbose, are the evidence.

### Untested / not re-run
- The 31 carried rc pins were not re-executed here (four representative probes were); their
  08-23 attributions and the adversarial verification of the `b42b71f` trio are cited, not
  redone. The 33 new-only arrivals were not individually re-run either: the pin mechanism
  (old's unknown-option r198) is uniform and one probe (`scr_fastnose`) demonstrates it.
- Storage-noise tiers were reproduced on `allest_aw` and `impv` only; the `ext_*`, `adm_*`,
  `cov_cwm`, `sg_attitime` magnitudes are read from the 0911 verbose tail.


---


# Part 7 - this run's own new material: the 31 new configs, their first-pass degeneracy stems, and the re-verification of Parts 1-6


**Scope.** Parts 1-6 were drafted against the 2026-09-11 dumps of a 507-config universe. Sections
7.0-7.4 below carry what the **2026-09-13 adjudication run** produced: the 31 configs new in it,
the degeneracy stems they surfaced, and the re-verification of Parts 1-6 against it. That run is
where the report's evidence was gathered, so it is kept as written.

**Section 7.5 is the run this report actually documents.** After the 2026-09-13 rulings landed in
the package, the suite was re-run end to end on 2026-09-14 over a 559-config universe - 21 more
configs, and a changed package. The header of this report quotes THAT run, not the one in 7.0.
Read 7.0-7.4 as the adjudication, and 7.5 as the closing verification.

### 7.0 The adjudication run's result (2026-09-13), stated first

    VERDICTS  DIFFERS=14 EXPLAINED=382 MATCH=26 RC_EXPLAINED=116
    COMPARED  538 of 538 configured
    SUMMARY   unexplained=56 explained=4597 additions=16517 stale_ledger=0 unscoped=0 prov_fail=0
    pool.log  106 shards OK, 0 TIMEOUT / INCOMPLETE / FAIL.  0 MISSING dumps.

**Not one pre-existing config moved.** All 14 differing configs are among the 31 new in this run.
Subtract them and the other 507 reproduce their prior verdicts exactly - 374 EXPLAINED, 25 MATCH,
108 RC_EXPLAINED - the same triple the pre-run judgement produced. That is the preregistration's
P0 (the two documentation commits move nothing) and P4 (plausible regressions in HEAD about zero)
both satisfied, and both were written down before the run rather than after it.

After adjudication the same dumps judge clean: **unexplained=0, stale_ledger=0, unscoped=0 over
538 of 538.**

### 7.1 The 14, adjudicated

**"New config" is not a classification.** Each still carries one of the report's four classes.

| class | count |
|---|---|
| old-bug-fixed | 10 |
| intended-change | 2 |
| mixed (both, on different objects) | 2 |
| both-builds-defect | 0 |
| unclear | 0 |
| **regression in the current build** | **0** |

On every one of the 14, HEAD is the more correct side or the author-ruled side.

| config | class | mechanism, and the commit that fixed it |
|---|---|---|
| `real_tt`, `trig_tt_year_rob`, `trig_tt_pois_rob` | old-bug-fixed | v1.3's uncentered `i.i##c.t` trend columns trip Stata's near-singularity screen, which zeroes the stage-1 VCE wholesale: `trace(e(V))` 0 against 261.68 / 232.79 / 5.3988, at identical k, rank and df_a. **f7019c9** |
| `bien_aw`, `real_aw`, `real_pw_uncond`, `real_cont` | old-bug-fixed | the estimand weights and the `ringer` copy were generated at float, so an intended exact collinearity cancelled only to ~1e-8. **3cb0c41** |
| `trig_bigids_month` | old-bug-fixed | float `_i` pooled 40 consecutive 3.0e7 IDs onto 20 float codes; stage-1 k 107 against 197. **8c5bfd7** |
| `real_pois` | old-bug-fixed | v1.3 hard-coded the `e(pwsgs)` stripe as `t` even under poisson, where margins reports `z`. **154d5dc** |
| `real_sg` | old-bug-fixed | v1.3 published a subgroup-contrast p of 2.64e-12 from a covariance rank-deficient in the contrast direction (implied correlation .9993); HEAD withholds it. **cfac3f5** |
| `bien_pois` | intended-change | `expression(predict(xb))` (numeric differentiation) to `predict(xb)` (analytic Jacobian). **32517df**, author-approved |
| `real_esf_adv` | intended-change | the adversarial null recalibrated; the stage-1 fit is bit-identical, so the whole move is the test. **BS-12 ruling 2026-09-08** |
| `ext_lat_dense_pid` | mixed | `cfx*_it_*` shape 3x1 to 3x6 old-bug-fixed (**3df74db**); the cfx value shifts intended-change (SI-20 **1df8613**) |
| `real_impv` | mixed | `ll` old-bug-fixed (**4de841d**); `df_r` 11 to 8 intended-change (**92e083b**, BS-10) |

**The strongest evidence.** Five of the seven float-related configs were bisected to a SINGLE site
in v1.3: patching only that `gen` or token in a copy of the old build reproduces HEAD bit-exactly
in `%21x`. And the timetrends mechanism was reproduced outside wooldid altogether - a raw
`i.i##c.t` fit gives `trace(e(V)) = 0` with a swallowed warning where the centered copy gives
22.48, at identical rank and k and with bit-identical residual sums. The variance was identified
all along and v1.3 withheld it.

**45 ledger rows in 14 blocks**, every one scoped to its config and pattern. No wildcards: the
project's rule is to name a new cell class in the comparer rather than blanket a matrix in the
ledger. That none of the 45 went stale on re-judging is the proof that each fired.

### 7.2 Degeneracy stems, adjudicated

The component is preregistered to FAIL on a first run over new configs, and it did: 9 failures, 0
stale, 0 bad rows, all 9 on the two new fixtures.

- **Six `bien_*`, class `DATASET_DESIGN`.** The event-study window is counted in contiguous
  relative-period SLOTS, not realized horizons, so a two-year grid spends every second requested
  slot on a period no cell realizes. Settled by a decisive experiment: a rank-compressed twin with
  N held at 1320 returns a full window where the biennial returns half, and doubling the biennial
  request reproduces the twin exactly. v1.3 is byte-identical on every shared spec. Documented at
  `sthlp:716` and predicted by the fixture's own build notes - but a user meets "biennial" 647
  lines earlier, which is now cross-referenced (commit `813c414`).
- **`real_sg`, class `RULED`.** The BS-15 rank-deficiency screen firing as designed and documented
  at `sthlp:852`. What was missing was the ledger row, not the ruling.
- **`real_boot`, `real_boot_t`, class `HIDDEN` and UNRULED.** All ten draws of one horizon fall
  below the observed estimate, so the bias-correction quantile is undefined and the column is
  skipped at `wooldid.ado:15172`. **This is the one finding here that is not documented anywhere:**
  `z0` appears at two lines in the help and neither names the undefined case, and unlike every
  sibling withholding in the package this one sets no `e()` macro saying so. For the author.

### 7.3 Shards that did not finish

None. `pool.log` carries no TIMEOUT, INCOMPLETE or FAIL line across all 106 shards. That is the
required answer, because the degeneracy component refuses to run at all if a dump is missing, and
a silent shard timeout would cost the whole run's stem adjudication rather than one shard's.

### 7.4 Re-verification of Parts 1-6 against this run

The package did not change between the dumps Parts 1-6 were drafted on and this run, apart from a
nine-line comment block and a shortened caution string, so every magnitude they quote should
reproduce. The audit that matters is recorded in `LEDGER_DELTA.md`: commit `4280f3b` replaced the
ledger hours after those sections were drafted, retiring 7 cited row ids, rewriting 26 reason
texts and closing 10 of the sections' open questions. Each is reconciled there and marked in
place; none was deleted. One caveat carried forward: `4280f3b` also changed the dump format, so
the re-tie compares value cells, not dump text.

### 7.5 The closing run (2026-09-14) - the one this report's header quotes

Sections 7.0-7.4 describe the run on which the adjudication was DONE. This section describes the
run that verifies it against the package the author's 2026-09-13 rulings produced. Both are
reported because they are different questions: 7.0 asks "what differs and why", 7.5 asks "does
anything still differ once the rulings are in".

    prov[old] = old|301b706
    prov[new] = new|04cc333
    VERDICTS EXPLAINED=411 MATCH=26 RC_EXPLAINED=122
    COMPARED 559 of 559 configured
    SUMMARY unexplained=0 explained=4626 additions=17069 stale_ledger=0 unscoped=0 prov_fail=0
    pool.log  112 shards OK, 0 TIMEOUT / INCOMPLETE / FAIL.  0 MISSING dumps.

**There is no DIFFERS line because the count is zero.** 411 + 26 + 122 = 559: every configured
config compared, every one of them either explained by a scoped ledger row, matching outright, or
explained by an rc pair. The comparer prints `DIFFERS=` only when it is nonzero.

**What changed between the two runs**, so the zero is not read as a re-run of the same thing:

| | 2026-09-13 (7.0) | 2026-09-14 (7.5) |
|---|---|---|
| configs compared | 538 | **559** (+21) |
| new side | `4280f3b`-era tree | **`04cc333`**, carrying the 2026-09-13 rulings |
| DIFFERS at first judging | 14 | **0** |
| shards | 106 | 112 |

The +21 are the coverage wave (18) and the three concentration configs (`real_perm_t`,
`real_perm_sd`, `real_boot_alpha`) that moved the riskiest single-fixture options onto a second
host, per the author's "move the riskiest onto a second fixture" ruling. The new side carries the
degrees-of-freedom fix, the weighted Bayesian bootstrap `P :/ N` correction, the fastnose posting
change and six smaller rulings.

**The 21 new configs produced nothing to adjudicate.** No new DIFFERS, and the degeneracy
component passed **on its first pass over them** - `DEGEN_CHECK_OK`, `fail=0 stale=0 badrow=0
classmoved=0`, with 175 curated expectation rows still summing exactly
(`curated=[BY_CONSTRUCTION=131 DATASET_DESIGN=11 HIDDEN=11 PART_RULED=13 RULED=9]`). That is
unlike 7.2, where a first pass over new fixtures correctly failed on nine stems. The difference is
that the coverage wave targets option SURFACES on fixtures already characterised, where the
2026-09-13 wave introduced two new fixtures.

**The one thing this run did surface: 46 ledger rows went stale, and 7 of them had blamed the
wrong cause.** The degrees-of-freedom fix made 46 rows unnecessary - the differences they excused
no longer exist. They were retired in `81b3ce2` rather than left in place, because a ledger row
that no longer fires is a row that will excuse the next real difference in its pattern.

Retiring them exposed a mis-attribution worth recording on its own. Seven rows carried the reason
id `R4W-fitstats`, attributing a root-mean-square-error divergence to the estimand-weight float
family of Part 3. It was not that. With the dof fix in, the same configs' rmse divergence falls
from **6.6e-4 to 3.0e-14** - ten orders of magnitude, and squarely the accounting change, not the
weights. The weight family is real and Part 3's evidence for it stands; what was wrong was
attaching these seven rows to it. **A ledger reason is a claim about mechanism, and this one was
checkable and wrong.** It is recorded here because the same mistake is cheap to repeat: a row
written while a second bug is live can blame the bug the author happened to be holding.

**Component verdicts for the closing run.**

| component | verdict |
|---|---|
| oldnew | 559/559, `unexplained=0`, `stale_ledger=0`, 0 missing dumps |
| cert | 43 / 43 |
| newfeat | 13 / 13 (each confirmed against the harness's own anchored `^[A-Z0-9_]+_ALL_PASS`) |
| coverage | `COVERAGE_OK` |
| degeneracy | `DEGEN_CHECK_OK` |
| battery | 46 of 48 bit-identical; the two predicted movers are treated in 7.6 |

**Honest note on how the run was judged.** `run_full.sh` reported `oldnew` and `battery` as
FAILURES on this pass. Both were expected and both were adjudicated rather than waved through:
`oldnew` failed only on `stale_ledger=46` - the retirement above, after which the same dumps
re-judge `stale_ledger=0` - and `battery` failed on exactly the two configs the dof fix was
predicted to move. Neither failure was a difference in a number the package reports.

### 7.6 The two battery movers, and why they were blessed

The degrees-of-freedom fix changes how many parameters an absorbing fit counts. Two of the 48
battery configs, `cont_lat` and `perm_cfxnose`, carry cells that depend on that count, and the
author's ruling ("fix the model, absorb the finest set only") carried a baseline re-bless with it.
The question was not whether to re-bless but whether the movement is the predicted one, and the
answer is a census rather than an assurance.

**The discipline came first.** The project's standing rule for any battery movement is
pristine-twice before per-config judgement. The battery was run a **second** time into a separate
fresh directory and the two runs compared **to each other at rtol 0**: `TOTAL_DIFF_GROUPS=0`, all
**48 configs bit-identical across the two runs**, `cont_lat` and `perm_cfxnose` included. The
movement is therefore deterministic, not run-to-run noise. After blessing, `test/baseline` compares
**48 IDENTICAL at rtol 0** against both that run and the full run's own battery output.

**What actually moved, counted cell by cell** (30 changed cells on `cont_lat`, 421 on
`perm_cfxnose`; the large count is the 20-replication `permdist_b` matrix, where every cell shifts
in its last digits):

| | `cont_lat` | `perm_cfxnose` |
|---|---|---|
| changed cells | 30 | 421 |
| at or below the project's 1e-10 indifference threshold | 24 | 410 |
| above it | 5 | 10 |
| `rmse` (the designed effect) | 6.63e-4 | 6.40e-4 |
| `rss` | **unchanged** | **unchanged** |

**The cells above the threshold, named individually, because a threshold breach that is not named
is a threshold that is not enforced:**

- `r2adj` and `r2_withinadj` on both configs, moving 2.1e-9 to 4.8e-9. These are **the designed
  effect, not an exception**: the degrees of freedom enter the adjustment formula directly, so a
  correct dof count must move them. Together with `rmse` they are the whole point of the change.
- `cont_lat`: **one** coefficient, appearing in three objects (`b[1,6]`, `b_es_att[1,2]`,
  `wdidmainresults[6,1]` are the same number). Its value is **8.10e-07** and its absolute change is
  **7.11e-15** - the double-precision floor for arithmetic over order-one quantities. Divided by a
  near-zero denominator that reads as 8.77e-09 relative. A relative tolerance is the wrong
  instrument on a cell whose denominator is itself numerical noise; this is the same tail
  amplification the ledger's own `SI20-P` rows describe for the comparison harness.
- `perm_cfxnose`: **8 `permdist_b` cells** at 1.45e-10 relative - absolute 1.08e-11 on values near
  0.075. These are permutation *replication* coefficients: each replication is its own fit, so each
  inherits the same last-digit shift the observed fit gets.

**The check that settles whether any of that reaches a user.** On `perm_cfxnose`:

- **All 56 permutation p-value cells are bit-identical.** Not close - identical. A 1e-11 shift in
  the replication distribution cannot move a p-value that counts replications against an observed
  statistic, and the count confirms it rather than assuming it.
- Every joint and adversarial statistic: unchanged.
- Worst movement across every reported coefficient and variance cell (`b`, `V`, `b_es_*`,
  `V_es_*`): **3.27e-11**, below the threshold.
- Permutation confidence-interval endpoints move at worst **2.50e-11**, below the threshold.

So the honest summary is: `rss` is unchanged and the fit is the same fit; the dof-dependent fit
statistics move because they are supposed to; every p-value a user reads is bit-identical; and the
only cells that trip a mechanical relative gate are two adjusted R-squareds that the fix is meant
to move, one near-zero coefficient whose absolute change is 7e-15, and eight permutation
replication cells whose absolute change is 1e-11.

### 7.7 How the closing green was actually obtained - provenance, stated because it is not one clean pass

The final component states below are real and each was verified. But `run_full.sh` did **not**
produce them in a single pass, and a reader entitled to the verdicts is entitled to how they were
reached.

**Two shards timed out on the run that produced these dumps.** `old_shard_covars_01` at its 600s
limit and `new_shard_real_slow_00` at its 1800s limit, because the suite was run at six jobs
concurrently with a four-shard Monte Carlo campaign and the box has ten cores. That orphaned four
configs' dumps - `cov_aw_disp` on the old side, and `real_boot_alpha`, `real_perm_sd` and
`real_perm_t` on the new side, the three concentration configs whose replication counts were
deliberately raised and which are therefore the slowest in the deck. It also turned three
`cov_aw_disp` ledger rows stale, since a row cannot fire against a dump that does not exist.

**Nothing about that was a defect in the package**, and the first judging said so plainly:
`DIFFERS=1` (the genuine `sg_degpre` insertion) alongside `MISSING=4` and `stale_ledger=3` that
were purely an artifact of how the machine was scheduled.

**The repair, and the two places it departs from a clean run.** Both shards were re-run to
completion **outside the pool**, each reaching its own `SHARD_DONE` sentinel, and all four dumps
verified present on both sides. The comparison was then re-judged and is the one quoted below.
Two honest caveats follow from that:

1. **`cert`, `newfeat`, `battery` and `coverage` passed earlier in the same run**, against a tree
   whose oldnew dumps were still incomplete. The file times say so and can be checked:
   `cert.log` 13:38, `newfeat.log` 13:47, `battery_compare.txt` and `coverage.txt` both 13:53,
   against `out/run_old` completing its repair at 13:59 and `out/run_new` at 14:51. That does not weaken them - none of the four reads
   those dumps; `cert` and `battery` run against the working tree, `newfeat` against its own
   fixtures - but they were not re-run after the repair.
2. **`degeneracy` was run by invoking `degeneracy_scan.py` directly**, because the `run_full.sh`
   wrapper refuses to scan when the last run left dumps missing, and `pool.log` is append-only:
   its two `TIMEOUT` lines persist even after the shards are repaired, so the guard stays tripped
   while the comparison is clean. The scanner's own output answers the question the guard exists
   to ask - `nodump=0` over all 559 configs - which is why running it directly is sound rather
   than a way around a red light. **The wrapper's guard cannot see a manually repaired shard**, and
   that is worth someone's attention as a suite observation - **it was ruled on and fixed the
   same day; see 7.8**. (`pool.log`, frozen at 13:53,
   carries the two `TIMEOUT` lines at its lines 29 and 112; `oldnew_compare.txt` is 14:52 and
   `degeneracy.txt` 14:53, both after the repair.)

**The lesson, which is about scheduling and not about the package.** The chunk's instruction is to
iterate until all tests pass; the first iteration here failed for a reason I introduced by
saturating the machine, and reporting it as anything other than that would misattribute my own
resource contention to the code under test.

### 7.8 Two defects in the suite's own tooling, found while closing out and fixed (`23400f9`)

Both were found by the closeout rather than by a test, both were put to the author, and both
were ruled FIX. Neither is a defect in wooldid; both are defects in the apparatus that judges
wooldid, which is the harder kind to notice because a broken judge reports success.

#### 7.8.1 The degeneracy guard could not see a repaired shard

`run_full.sh`'s degeneracy precondition refused to scan when EITHER `pool.log` carried a
`TIMEOUT` line OR the comparer reported a missing dump. `pool.log` is **append-only**. So once
a shard had timed out, the first condition was tripped permanently - including after that shard
was re-run by hand to its own `SHARD_DONE` sentinel and every dump verified recovered. The
wrapper then refused a scan that the scanner itself passes at `nodump=0`, and the only way to
close the run was to invoke `degeneracy_scan.py` directly, which is precisely the bypass a
guard exists to make unnecessary. That is what happened on 2026-09-14 and it is why 7.7 has to
describe the closing scan as a direct call.

The comparer's own MISSING report is the correct condition **alone**: it is recomputed from the
dumps on disk at judging time, so it clears when a repair really lands and stays tripped when
one does not. The stale-dump case the `pool.log` condition was guarding against is ruled out by
that same re-judge. Dropped the `pool.log` clause. **Proven** by running `./run_full.sh
degeneracy` with the identical stale `pool.log` still on disk: it now passes the scan it
previously refused.

#### 7.8.2 The ledger join credited rows written about something else

This is the more interesting one, and it is a defect in a **classifier**, so its failure mode is
a confident wrong answer rather than an error.

`degeneracy_scan.py`'s `ledger_known` credited degeneracy coverage from any ledger row matching
the same e() **matrix**. A config whose dead objects were all covered was promoted to RULED.
The promotion therefore did not require anyone to have adjudicated the degeneracy - only that
some row, written for any reason at all, happened to touch the same matrix.

Two configs were being classified on that basis:

- **`sg_degpre`**: its `WINDOW_SHORT(att_s2:pre0/2)` - subgroup 2 has no pre-treatment row - was
  promoted to RULED by `F1SG-degpre-m-b`, a row about a **positional insertion in `b`** that
  says nothing whatever about subgroup 2's window.
- **`ext_cov10_rcs_pweight_ddd_c1`**: the same shape, via `DDG-pw-1`.

What makes this hard is that **nothing mechanical separates those from the honest cases**.
`BS15-pw-esptrip` really is about `esptrip`'s `pwsgs` deadness; it has the same `kind`, the same
shape of `pattern`, and names the same matrix. The difference lives only in the prose, so no
glob or kind filter can recover it - a human has to assert the tie.

**The fix, and a correction to the mechanism first proposed.** The first design was to require a
row's `pattern` to name the dead object's stem. That is not available: `pattern` is also
`compare_oldnew.py`'s cell-match key (`compare_oldnew.py:216`), so narrowing it would have
stopped those rows explaining the cells they exist to explain and turned the 559-config
comparison red. Checked before anything was edited.

What shipped instead uses a field that already exists for exactly this purpose.
`degeneracy_expected.csv`'s `citation` column for a RULED/PART_RULED row **is** a human saying
"these ledger ids rule this config's surfaces" - it was simply never used as anything but a
danglingness check. `ledger_known` now takes a `candidates` set and credits coverage only from
the cited ids. No ledger edits, no effect on the comparer, and no machine-parsing of a free-text
field. The all-ledger pass survives as a `DEGEN_NOTE` so "the ledger now covers this, consider
citing it" stays visible without being allowed to classify anything.

**What it moved, measured before the rule was changed.** 27 configs were RULED or PART_RULED and
22 of them were credited on a dead-row object, which made the change look enormous. Under the
citation gate exactly **three** moved - `cont_sg_tidose`, `cont_sg_gap`,
`ext_san_input_missing` - all RULED to PART_RULED, all on an uncovered `mat:b`. All three turned
out to be honest rows the citation list had simply omitted: `BS8-cont-sg-tidose-b`,
`BS8-cont-sg-gap-b` and `SI20-L-ext_san_input_missing-mat-b` each explicitly assert that
config's own `e(b)` is non-estimable. Citing them restored RULED, now as an explicit assertion
rather than a coincidence of matrix names.

**`sg_degpre` went back to `DATASET_DESIGN`**, which is what it has been all along: subgroup 2
has no pre-treatment row because the fixture was built that way (`data_classes.csv`), and the
ledger says so in `F1SG-degpre-n-sg` ("Subgroup 2 genuinely has no pre cells and is still
withheld"). It was only ever moved to RULED because `class_compatible()` rejects
`DATASET_DESIGN` against a computed RULED, so the curated file had to be edited to match a
promotion nobody had made. Census back to **11 DATASET_DESIGN / 9 RULED**.

**A correction to my own framing.** When this was put to the author I described it as fixing the
class of error behind "R4W-fitstats, DDG-pw-1 and sg_degpre" - three instances. R4W-fitstats is
not one of them. Those rows were retired from the ledger and were never cited in
`degeneracy_expected.csv`; they were a **comparer-side** mis-attribution, blaming the AMESE R4
refit removal for an rmse/r2adj divergence that was really a degrees-of-freedom over-count. The
same shape of error in a different tool, untouched by this fix. For the scanner join the count
is **two**.

**Verification.** `./run_full.sh degeneracy` returns `DEGEN_CHECK_OK` with `fail=0 stale=0
badrow=0 classmoved=0` and a curated census matching `DEGENERACY.md`. `./run_full.sh compare` is
unchanged at `EXPLAINED=411 MATCH=26 RC_EXPLAINED=122` over 559/559 with `unexplained=0
stale_ledger=0`, which is what proves the ledger and the comparison were not disturbed.



---


# Part 8 - external validation of the new features: what four external packages could and could not benchmark, and the issues that remain open


**Scope.** This is the external-validation part of the 2026-09-13 assessment. Four external
benchmarks were probed for their usability as oracles for `wooldid`: **jwdid** (Rios-Avila's
Wooldridge/Mundlak extended-TWFE estimator), **ritest** plus Stata's own **`permute`** prefix,
Stata's own **`bootstrap`** prefix plus **`boottest`**, and the **staggered-DiD references**
(honestdid, csdid, did_imputation, xtevent). All four probes returned the verdict **exact** —
not "qualitatively similar", not "agrees within Monte-Carlo error" — each under an option
mapping the probe had to find and state. The value of that is not a warm feeling about
correctness; it is **isolation**. When two independently written estimators agree on the point
estimates to 1e-14 on the same data with the same clustering, and then disagree on the variance,
the variance path is the only thing left to blame. That is exactly what happened in one place
(§C.6), and it is why an unruled candidate comes out of a feasibility study.

**Provenance.** `git diff 08ec0b4 HEAD -- wooldidinstall/` is empty and the working tree carries
no modification to `wooldidinstall/` (only an untracked `wooldid.ado.bak`), so every
`wooldid.sthlp` and `wooldid.ado` line number below holds at the probes' pin `08ec0b4`, at HEAD
`21db885`, and in the working tree alike. **Every wooldid line number cited here was re-opened
and read in this pass**, and so were the installed externals (`jwdid.ado`/`jwdid.sthlp` v1.41 /
v1.1, `ritest.ado` 1.1.7, `/Applications/Stata/ado/base/p/permute.ado`, `honestdid.ado` 1.1.2);
**six** citations needed correcting and are corrected in place, and **two** measured-only claims
are additionally pinned to their source; the consolidated list is §C.9.
Log-line citations (`run1.log:17977`, `bench_run2.log:674`, …) are **relayed** from the four
probe reports and were not re-opened here. **No Stata was run for this section**; where a claim
needs a run it is marked and the run is named.

**Sources.** `stageC_feas/jwdid/REPORT.md` (312 lines), `stageC_feas/ritest/REPORT.md` (506),
`stageC_feas/bootstrap_prefix/REPORT.md` (555; byte-identical to the frozen
`REPORT_agent3_final.md`, checked with `cmp` this pass), `stageC_feas/staggered_refs/REPORT.md`
(527), indexed by `stageC_probe_returns.json`.

---

### C.1 The four oracles

| benchmark | oracle **for** | mapping | measured agreement | the one thing it cannot benchmark |
|---|---|---|---|---|
| **jwdid** 1.41 + `jwdid_estat` 1.42 | the whole `g×t` cell matrix, `att`, every ES horizon, per-cohort and per-calendar aggregates, **and their clustered V** | give wooldid a per-onset cohort id as `i` (`icoh = cond(ttre<., ttre, i)`) and the same id as jwdid's `ivar()`, `cluster(i)` explicit on both; wooldid default ↔ jwdid default, `esfixedbaseperiod esrelativeto(-1)` ↔ jwdid `never` | b **and** V 1e-18 … 1.5e-13 over 10 legs, balanced or not, weighted or not, incl. `dta_s_admin` (`e(N)` = 29,067 on both sides) | wooldid's **unit-level-`i` clustered V** — INFO-only in every leg (SE ratios 0.62–1.24); jwdid is not an oracle for it anywhere |
| **ritest** 1.1.7 + Stata 19 **`permute`** | the permutation **assignment space** and every replication draw | ritest `cluster(i)` ≡ `permute(shuffle)`, `strata()` ≡ `permstrata()`, `samplingsourcefile()` replays `permsave()`; `permute` needs a one-row-per-cluster frame behind an `eclass sortpreserve` wrapper; only `permrawp` maps to an external p | `permute, enumerate` vs `permexact`: 20/20 draws **bitwise**, p_lower 1 / p_upper .05 / p_two .1000 equal digit for digit; ritest replaying wooldid's own draws: 199/199 **bitwise**, `r(p)` = .8592964824120602 = 171/199 = wooldid's `permrawp` share | wooldid's **default and `permabsp` p-values** (they carry an add-one; neither external engine has one nor offers one) and **`permtstat`**, never benchmarked |
| Stata **`bootstrap`** prefix (+`bstat`, `estat bootstrap`); **`boottest`** 4.4.10 | `bootstrap(standard)`: every replication draw, the bootstrap SE, the percentile and BC intervals, `z0`, and the `boottstat` re-identification leg | wooldid `rseed(S)` = prefix `seed(S)` puts both on **the same RNG stream** (one `bsample, cluster() idcluster()` per rep); drive `wooldid … fastnose`; `saving(, double)` mandatory | 796 draws at 8.10e-15; percentile/BC/`z0` **bit-identical** on common draws; `bootsave()` vs the prefix re-id leg at one seed: 200/200, `maxreldif_b` and `maxreldif_se` **exactly 0.000e+00** | `boottest`'s **p-value** — all six measured p's sat at the 1/999 floor; and wooldid's percentile p, which the prefix does not compute at all |
| **csdid** 1.72 / **did_imputation** Nov-2023 / **xtevent** 2.2.0 `impute()` / **honestdid** 1.1.2 | the `esfixedbaseperiod` ES path (csdid, xtevent), the pooled-reference post horizons (BJS), and the honestdid handoff on both routes | csdid `long2` + `estat event` ↔ `esfixedbaseperiod esrelativeto(-1)`; BJS `horizons(0/H) minn(0)` ↔ `esprelength(0) espostlength(H+1)`; xtevent `impute(stag) window(-P Q) vce(cluster i)` with binned endpoints; honestdid route 1 `b()/vcov()/numpre()` vs route 2 `pre()/post()` on the posted `e(b)`/`e(V)` | csdid 1.02e-14 (100 MC draws × 2 DGPs × 14 horizons); xtevent 3.8e-15 absolute / 1.1e-12 relative; BJS 1.16e-07 (two solvers' FE tolerance); honestdid route 1 ≡ route 2 `mreldif` 0 under `delta(rm)`, 2.71e-09 under `delta(sd)` | **standard errors below admin scale** — cross-tool SE ratios 0.80–1.74 on 40-unit panels against 0.995–1.023 on `dta_s_admin`; and p-values in the imbalanced regime (a 0.42-SE agreement moved p by 0.183) |

**Read the last column as the boundary of the claim.** Three of the four oracles are exact in the
points and silent on the variance. The fourth (jwdid) is exact in the variance too, but only in
the cohort-id parameterization. Nothing here validates wooldid's unit-level clustered variance
against an independent variance engine — see §C.7(ii) and §C.6, which are the same fact seen
from two sides.

---

### C.2 wooldid **help** defects (four)

#### C.2.1 The `e(b)` coefficient names are never stated, while the help twice tells users to use them

- **What is wrong.** `main_att` occurs **0 times** in `wooldid.sthlp`; `pretreatment_att` occurs
  **0 times**; `ES_att_pre#`/`ES_att_post#` occur exactly once, inside an `event_plot` stub
  example at **sthlp:1422**. Yet **sthlp:484** (the `permsave()` entry) writes "one column per
  coefficient of `e(b)` under its `e(b)` name; under `permtstat` also …", **sthlp:541** (the
  `bootsave()` entry) writes "the replication coefficients under their `e(b)` names". (The
  bootstrap probe also lists `estout`/`coefplot`/`honestdid` as pointed at those names; re-checked
  this pass, `coefplot` occurs **0 times** in the help and the `estout` passages at sthlp:283 and
  :438 are about `e(p)`/`e(ci_l)`, not coefficient names — so the instruction is line-cited
  **twice**, at 484 and 541, not three times.) The only place the help names
  anything coefficient-like is **sthlp:1281**, and what it names there are the *row labels of
  `e(wdidmainresults)`* — `main:att`, `pretreatment:att`, and the `ES_att`/`ES_att_it` equations
  "named `post#` and `pre#`". Those are a **different naming scheme** from the `e(b)` column
  names, so the one piece of naming the help does supply actively invites `_b[main:att]`.
- **File and line.** `wooldidinstall/wooldid.sthlp:484`, `:541`, `:1281`, `:1422`; the
  "display purposes only" caveat at `:1437-1439` does not cover it (it is about the covariance,
  and the same help keeps directing users to the names).
- **Measured.** The occurrence counts above were re-counted in this pass with
  `grep -c` / `grep -on 'ES_att[A-Za-z_#]*'` on the working tree (identical to `08ec0b4`). The
  user-facing consequence was measured by the bootstrap probe: `bootstrap b=_b[att]` fails
  **r(111)** "error in expression: `_b[att]`" — recorded at `run1.log:18279` — and it cost that
  probe a whole run.
- **Who it affects.** Every user who leaves wooldid's own tables: `test`, `lincom`, `estout
  keep()`, `coefplot`, `bootstrap b=_b[…]`, a `boottest` hypothesis, or any script that reads
  `e(bootdist)`/`e(permdist_b)` columns by name. Also every Stage C test in §C.8.
- **Fix.** Name the four stubs once, in the Stored Results entry for `e(b)`: `main_att`,
  `pretreatment_att`, `ES_att_pre#`, `ES_att_post#` (and the `ame` twins). One sentence.
- **Citation correction.** The bootstrap probe labels `sthlp:484` "permtstat"; line 484 is inside
  the **`permsave()`** option entry (which mentions `permtstat` in the same sentence). The line
  and the claim are right; the entry name is not.

#### C.2.2 The `did_imputation` equivalence conditions omit `esprelength(0)`

- **What is wrong.** `wooldid.sthlp:1461-1465` ("Comparing Output to did_imputation") lists four
  conditions for wooldid to reproduce `did_imputation`; the sentence carrying them is
  **sthlp:1464**: "given that one uses the `att` estimand-type, uses a pooled reference period
  (i.e., not using `esfixedbaseperiod`), does not use unconditional standard errors, and does not
  use controls". A pooled reference is **necessary but not sufficient**: the reference pool is
  the set of horizons the run asks for, so asking for pre horizons removes them from the pool and
  shifts every post coefficient.
- **File and line.** `wooldid.sthlp:1464` (condition list), entry `:1461-1465`.
- **Measured, and where recorded.** A dedicated confirmation run on the `cohom` panel
  (common onset t=7, 40 units, T=12): wooldid `esprelength(0)` vs BJS max |diff| **4.0738e-08**
  (`EQUIVALENCE_HOLDS`); wooldid `esprelength(4)` vs BJS **4.7619e-02 at every horizon**
  (`EQUIVALENCE_BROKEN`); and the clincher — `did_imputation`'s own `pretrends(4)` moves its post
  taus by **exactly 0.0000e+00** (`POST_TAUS_UNCHANGED`). Recorded at
  `stageC_feas/staggered_refs/espre0.txt:36-48`, log `espre0.log`. The shift is **not constant
  everywhere**: 4.76e-2 on cohom/cohet, 3.6e-4 … 1.17e-2 on costag, up to **7.05e-2** on
  `dta_binary` (`analysis.txt:267, 311, 359, 407`).
- **Who it affects.** Precisely the user who does the natural thing — runs both packages side by
  side **as an event study**, i.e. with pre horizons — and sees a uniform gap with nothing in the
  help to explain it. Asking the reference tool for pre-trends costs nothing; asking wooldid for
  pre horizons costs 0.0476.
- **Fix.** Add `esprelength(0)` to the condition list at sthlp:1464.

#### C.2.3 The `e(b_es_type)` entry claims `fastnose` posts a coefficient vector; it does not

- **What is wrong.** **sthlp:1295**, inside the `e(b_es_type)` / `e(V_es_type)` Stored Results
  entry that begins at sthlp:1288, reads: "When NO horizon has an identified variance
  (`impvalues`, `fastnose`, typically `timetrends`, and any run whose estimation returns an
  all-zero covariance), the covariance is not returned but the coefficient vector still is, under
  a Caution". Under `fastnose` the coefficient vector is **not** returned.
- **File and line.** `wooldid.sthlp:1295`. It contradicts the help's own text in three other
  places: **sthlp:276** (`impvalues`: "`e(b_es_)`, `e(V_es_)` and `e(numpre_)` are not returned
  and the `honestdid` and `pretrends` handoffs have no inputs"), **sthlp:284** (inside the
  **`fastnose` option entry itself**, sthlp:278-286: "`honestdid` and Roth-style diagnostics have
  no inputs"), and **sthlp:438** (the `permute()` "Fast-mode default" passage: "`e(V)`,
  `e(b_es_*)` and `e(V_es_*)` are then not posted").
- **Measured, and where recorded.** `wooldid y i t ttre, att cluster(i) esprelength(5)
  espostlength(6) esfixedbaseperiod fastnose` → `confirm matrix e(b_es_att)` returns **rc 111**,
  `e(numpre_att)` empty, `e(wdid_vcvstate)` = `absent`; the identical run without `fastnose`
  returns rc 0 with `numpre_att` 5 and `vcvstate` `analytic-blocks`. Recorded at
  `stageC_feas/staggered_refs/xtdiag.txt:175-181`. The pooled-reference `fastnose` run behaves
  the same way.
- **Who it affects.** Anyone building an MC or a test over `fastnose` runs — they must read the
  `ES_att_*` columns of `e(b)`, which **do** post. This defect has already cost real work: the
  staggered-refs probe's 100-draw MC captured the `i = cohort` arm from `e(b_es_att)` and
  therefore has **no 100-draw denominator** for that arm (§C.6).
- **Fix.** Either drop `fastnose` (and, on the evidence of sthlp:276, `impvalues`) from the
  sentence at 1295, or fix the three code paths. Three other help lines (276, 284, 438) and one
  measurement (`xtdiag.txt:175-181`) stand against line 1295, so the sentence is the likely error
  — but this is a **text ruling for the author**, not a mechanical fix.
- **Citation correction.** The staggered-refs probe calls sthlp:284 "the help's `permute()`
  fast-mode paragraph". Lines 278-286 are the **`fastnose` option entry**; the `permute()`
  fast-mode passage is at 437-438. Both lines say what the probe says they say; the labels were
  swapped. The correction strengthens the finding: the `fastnose` entry contradicts line 1295
  about `fastnose`.

#### C.2.4 The honestdid example omits `parallel(0)`

- **What is wrong.** **sthlp:1534** shows
  ``honestdid, numpre(`=e(numpre_att)') b(b) vcov(V) delta(rm) mvec(0.5(0.5)2)``. On any machine
  with the `parallel` package installed, honestdid defaults to **`parallel(4)`** and a user
  copying the example spawns four child Stata processes per call — unexpected, and multiplicative
  inside a loop.
- **File and line.** `wooldid.sthlp:1534` (example 10, the fit at :1528 and the matrix extraction
  at :1530/:1532).
- **Measured, and where recorded.** The default was measured by the probe (honestdid 1.1.2, the
  `parallel` package installed on this machine) and **pinned from source in this pass**:
  `honestdid.ado:128-130` reads ``if ("`parallel'" == "") { cap parallel numprocessors ; local
  parallel = cond(_rc, 0, 4) }`` — 4 whenever `parallel` is installed, 0 only when it is not; the
  package author's own `* TODO: xx Should default be no parallel?` sits at `:132`. Cost with
  `parallel(0)`: 7.2–8.1 s per `delta(rm)` call at k=8, `mvec 0.5(0.5)2`
  (`staggered_refs/results.txt:144, 151` for the two `delta(rm)` calls and `:189, 197, 216` for
  the sentinel `rm` runs; `:164, 171` are the `delta(sd)` timings, 0.084–0.219 s).
- **Who it affects.** Users copying the example; and it is the same omission that made the suite
  fan out. **Suite status, corrected:** the probe reported `test/suite/newfeat/honestdid_grid.do:126`
  and `honestdid_crosscheck.do:73` as fixed-but-uncommitted in the working tree. Re-checked this
  pass: both lines now carry `parallel(0)` **in commit `4280f3b`**, an ancestor of HEAD, and both
  omit it at `08ec0b4`. It is a landed fix, not an outstanding item, and it is **not counted as a
  defect** here.
- **Fix.** Add `parallel(0)` to the example at sthlp:1534, or one clause naming honestdid's
  default. The CI the example produces is correct; this is hygiene.

---

### C.3 wooldid **code** defects (three)

#### C.3.1 The standard r(459) trailer prints twice on a validation failure

- **What is wrong.** On a bootstrap validation failure the user sees wooldid's own `Error: …`
  line followed by **two** copies of Stata's standard trailer, "something that should be true of
  your data is not".
- **File and line — and the mechanism, which the probe hypothesised and this pass diagnosed from
  source.** `wooldid.ado:6217` prints the custom message and `:6218` raises `error 459`, inside
  `wooldid_bootctl`. The call site is `wooldid.ado:4769`, `cap noi wooldid_bootctl , …`, followed
  by **`wooldid.ado:4770`: `if _rc != 0 error _rc`**. The `noi` lets the subprogram's own error
  print the trailer once; `error _rc` then **reprints** it. The file documents this exact hazard
  itself at `wooldid.ado:61` — "rc and error messages pass through exactly (**exit, not error -
  error reprints**)" — which is why the shell uses ``exit `wdsh_rc'`` at `:100` and prints nothing
  twice. The probe's guess ("presumably the `error 459` surfacing through both the `wooldid` shell
  and `wooldid_worker`") is **wrong about the shell**: the shell is correct; the worker's call
  site is not.
- **Scope is wider than measured.** A survey of all 14 `cap noi wooldid_*` call sites in the file
  found **exactly two** that re-raise with `error` rather than storing the rc: `:4770`
  (`wooldid_bootctl`) and **`:4756` (`wooldid_permctl`, called at `:4755`)**. So every validation
  refusal raised inside either control program double-prints — including `permute()`'s own
  nesting refusal at `wooldid.ado:5844` and its `permsave()` r(198) at `:5698`. The `bootctl`
  case is **measured**; the `permctl` case is **reasoned from source, not measured**.
- **Measured, and where recorded.** `wooldid y cohid year ttre, att cluster(state) … boottstat`
  on `dta_s_admin`, where `cluster(state)` does not nest `i = cohid`:
  `bootstrap_prefix/run5.log:17818-17821`. Contrast the r(111) case, which prints once
  (`run1.log:18279`).
- **Who it affects.** Everyone who trips a bootstrap or permutation validation screen, and every
  batch log that records one. Cosmetic; nothing is wrong with the rc or the custom message.
- **Fix.** `exit _rc` in place of `error _rc` at `wooldid.ado:4770` and `:4756` — the idiom the
  file's own comment at `:61` prescribes and the shell already uses.

#### C.3.2 The low-`reps()` Caution prints **before** the validation error that aborts the run

- **What is wrong.** On the same failing call the user is advised about p-value resolution on a
  run that never executes.
- **File and line.** Both blocks live in `wooldid_bootctl`: the Caution at
  `wooldid.ado:6016-6018` (three branches: studentized, joint/equality/adversarial, and the
  coefficient-grid form), the nesting assertion at `:6208-6219`. The Caution is simply upstream.
- **Measured, and where recorded.** `bootstrap_prefix/run5.log:17810-17821` — the Caution at the
  top of the block, the refusal at the bottom. Confirmed by reading the two line ranges this pass.
- **Who it affects.** Anyone reading a failed run's log; the advice is not merely useless but
  misleading, since it implies the run got as far as choosing a replication count.
- **Fix.** Move the `6016-6018` Caution block below the validation checks that can abort
  (`:6205-6220` and the cluster-count sanity block that follows at `:6223`).

#### C.3.3 The `permsave()`/`bootsave()` name check rejects only `/`, `\` and empty, so `permsave(name, replace)` is swallowed whole as a filename

- **What is wrong.** `permsave()` and `bootsave()` take a **plain file name**, and `replace` is a
  separate top-level wooldid option — deliberately *not* Stata's `saving(f, replace)` idiom. The
  name validation rejects only a forward slash, a backslash and the empty string, so
  `permsave(wdB0, replace)` — the obvious analogy every Stata user reaches for — is accepted
  **whole**, and the run reports success writing `wdB0, replace.dta`.
- **File and line — four guard sites, all with the same three-condition test.** The probe cited
  two; this pass located all four:
  - `wooldid.ado:5697-5700` — `wooldid_permctl` pre-flight (``if strpos("`wdpfn'", "/") |
    strpos("`wdpfn'", "\") | "`wdpfn'" == ""`` → r(198));
  - `wooldid.ado:8513-8516` — `wooldid_permsave`, the writer;
  - `wooldid.ado:6042-6045` — `wooldid_bootctl` pre-flight, `bootsave()`;
  - `wooldid.ado:6288-6291` — `wooldid_bootsave`, the writer.
  Each is preceded by the trailing-`.dta` strip (`:5696`, `:8512`, `:6041`, `:6287`), which is why
  the probe's `5696-5700` / `8513-8517` ranges are off by one line at each end.
- **Measured, and where recorded.** The artefact `wdB0, replace.dta` is on disk in
  `stageC_feas/ritest/`; the run printed `permsave(): the permutation draws were written to
  wdB0, replace.dta`. It then **cost that probe the remainder of a run**: `runB` aborted at
  `use wdB0.dta` with r(601) and blocks B2 onward had to be rebuilt as `runD`.
- **Who it affects.** Any user who writes the standard idiom. Severity **low** — nothing is lost
  and the next `use` fails loudly, which is how it was found — but it is a silent surprise where
  an r(198) belongs.
- **The help is not at fault.** `wooldid.sthlp:485` says "it is still a plain file name, with no
  directory part", and `:482` says "(a plain file name; an existing file needs `replace`)". This
  is a validation gap, not a documentation one.
- **Fix.** Reject a name containing a comma at all four sites, or name the `, replace` suffix in
  the error text at `:5698`, `:8514`, `:6043`, `:6289`.

---

### C.4 A defect in the **external** tool, not in wooldid (one)

#### C.4.1 jwdid's help documents `[pw]`; the installed code returns r(101) for any weight

- **What is wrong.** `jwdid.sthlp` v1.1 declares `{weight}` in its syntax diagram (**line 16**)
  and states "`pweight`s are allowed; see `help weight`" (**line 35**). Installed `jwdid.ado`
  **v1.41** (file dated 2024-02-07) returns **r(101) "weights not allowed"** for both `[aw=]` and
  `[pw=]`.
- **Diagnosed from source in this pass.** The program's **first** parse is the replay guard at
  `jwdid.ado:41`, `syntax [ anything(everything)], [example *]` — which carries **no weight
  specification**, so Stata's parser refuses the weighted call there. The real estimation syntax
  at `jwdid.ado:53` **does** declare `[pw]`, and the estimation line at `:202-203` passes
  ``[`weight'`exp']`` straight through to `reghdfe`. So the code was written for weights and the
  replay guard shadows it. This is a jwdid code defect, not a help error.
- **Measured, and where recorded.** Reproduced twice: `jwdid_crosscheck.log:1695-1698`
  (`J6_JWDID_PW_RC=101`, on the finalize re-run) and by a standalone `_wtprobe` in `bench.log`.
- **Consequence for wooldid.** The weighted benchmark leg cannot be written against `jwdid` +
  `estat` directly; it must emulate jwdid's own linear path (`reghdfe y <cells> [pw=w],
  abs(icoh t) cluster(i) keepsingletons`, `jwdid.ado:202-203`, plus `estat`'s aggregation as
  explicit weight vectors). Under that emulation wooldid `[aw=w]` and `[pw=w]` are exact
  (b ≤ 3.3e-15, V ≤ 1.3e-12 on `dta_binary`; bitwise V on the common-onset panel) and are bitwise
  identical to each other.
- **Fix.** None available in this repo. The adopted Stage C script already prints the return code
  as `J6_JWDID_PW_RC` rather than judging it, so a later jwdid release that accepts weights
  announces itself and the leg can be switched to a direct `estat` comparison.

---

### C.5 An undiagnosed **performance observation** — not a defect

- **What was measured.** On `dta_binary` (960 obs, 80 clusters, 12 periods), wooldid's own
  `permute()` loop costs **1.86 s per replication** (369.5 s / 199 reps) where an external prefix
  whose per-draw command is `wooldid … fastnose` costs **0.95–0.97 s/rep** (ritest 189.3 s;
  `permute` through the cluster-frame wrapper 193.3 s), and the ritest replay — which draws no
  assignment — costs **0.66 s/rep**. That is a ~2× per-replication gap against engines doing the
  same work.
- **What was NOT measured.** The split between fixed setup and per-replication overhead. An
  earlier draft's claim that "180 s of the 369 s is not in the replication loop" was an inference
  and the probe **withdrew** it. A standalone `wooldid … fastnose` call on this fixture timed at
  1.44 s and a full call at 1.80 s (`ritest/A_results.txt`, `PROBE_*`), but those are standalone
  calls rather than warm loop iterations and are not a like-for-like third point.
- **A cross-probe narrowing, available only by putting the two reports together.** On the *same
  dataset*, the same prefix-vs-wooldid comparison runs the **other way** for the bootstrap engine:
  wooldid `bootstrap(standard)` 127.8 s at R=200 (**0.639 s/rep**) against the prefix over
  `wooldid … fastnose` 145.5 s (**0.727 s/rep**) — `bootstrap_prefix` run5 timers 8 and 9. So
  this is **permute-loop-specific**, not a general wooldid loop overhead, which is a useful
  narrowing for whoever profiles it. Reasoned from two measured pairs, not itself measured; the
  two legs also differ in their option strings (the bootstrap pair carried
  `esprelength(1) espostlength(1)`).
- **Caveat that binds all of it.** Every timing in both probes was taken on a machine running
  three to five other Stata jobs, and `run5`'s timers 8/9 overlapped `run6`. The absolute seconds
  are upper bounds; only the ratios are meaningful. **Label: one profiling pass, not a defect.**

---

### C.6 CANDIDATE — wooldid standard errors undercover in two arms, with the external tools on identical points as the control

> **STATUS AT CLOSE (2026-09-14): closed by RULING, not by measurement, and the distinction
> matters.** This section had two halves. The **help half** - `wooldid.sthlp:67` steering users
> toward the `i` geometry that these draws show undercovering - was put to the author on
> 2026-09-13 and he ruled: *"Let's remove this statement from the help file and call it a day; I
> am skeptical of it."* The claim was removed, so the help no longer points anyone at the arm in
> question, and that half is closed.
>
> The **statistical half is not closed and is not measured.** No 500-draw MC was run. Everything
> below still rests on 30 draws at a single G, and the shape the section itself specifies as the
> minimum before anyone rules (500+ draws, G in {20, 40, 200}, a coarse-`i` arm, csdid and
> did_imputation as controls in every cell) remains unbuilt. Recorded plainly so that the removal
> of a help sentence is never later read as evidence that the SEs were checked and found
> correct. **This is one of the items going to the author as open at the close of this
> assessment.** Cost estimate stands at roughly half an hour of machine time; what it lacked was
> not feasibility but a slot.


**Why this is here at all.** csdid produces the *same point estimates* as wooldid to 1.02e-14 and
did_imputation to 1.16e-07, from the same data with the same `cluster(i)`. Same data, same
points, same clustering, different variance formula: the benchmarks isolate the variance path
exactly. That is the entire value of an exactness oracle, and it is what turned a feasibility
probe into a number-mover.

**Exactly what was measured.** `se_mc.do`: **30 draws × 2 DGPs**, both 40 units × 12 periods,
`cluster(i)`, **G = 40**, full analytic SEs, five estimators per draw, dumping `e(b)` and the
`e(V)` diagonal. DGPs: `hom` (common onset t=7, τ = 0.5) and `stag` (two onsets 5/8 plus
never-treated, τ heterogeneous with a trend). Empirical SD denominators come from **two
independent MCs** — the 30-draw one (`semc_out.txt`) and the independent **100-draw** point-parity
MC (`mc.txt`, re-derived this pass into `mc_analyze_out.txt`). Coverage is against the DGP's true
horizon effects and does not depend on the SD estimate at all.

| arm (40 units, `cluster(i)`, G = 40) | DGP | meanSE / SD₁₀₀ (min/med/max) | mean 95% coverage |
|---|---|---|---|
| **wooldid `esfixedbaseperiod`, `i` = unit id** | **stag** | **0.57 / 0.79 / 0.88** | **0.838** |
| wooldid `esfixedbaseperiod`, `i` = cohort id | stag | (no 100-draw arm; 30-draw 0.88/0.98/1.12) | 0.939 |
| csdid `long2` — *identical points* | stag | 0.90 / 1.01 / 1.23 | 0.928 |
| **wooldid pooled reference (`esprelength(0)`), `i` = unit** | **stag** | **0.58 / 0.62 / 0.74** | **0.800** |
| did_imputation — *same points to 1e-7* | stag | 0.89 / 0.93 / 1.08 | 0.929 |
| **wooldid pooled reference, `i` = unit** | **hom** | **0.77 / 0.82 / 0.86** | **0.903** |
| did_imputation | hom | 0.92 / 0.97 / 1.00 | 0.947 |
| wooldid `esfixedbaseperiod`, `i` = unit | hom | 0.86 / 1.00 / 1.09 | 0.947 |

(**S1**) is the staggered `esfixedbaseperiod` unit-`i` row: SEs 20–40% below the sampling SD,
coverage 0.838 at a nominal 0.95, `main_att` worst at 0.77–0.79. The same run with `i` = the
cohort id is clean, csdid on the identical points is clean, and the **common-onset** unit-`i` run
is clean — so it is staggered-specific. (**S2**) is the pooled-reference path, low in **both**
DGPs, with did_imputation clean in the same draws.

**Why it is not noise, stated as the probe stated it.** The SD is measured twice independently
and the shortfall survives both denominators; the ratios are 3+ MC standard errors from 1; and
the control estimators sit at 1 in the same draws. On `dta_s_admin` with a coarse `i` wooldid's
SEs match csdid's to 0.7–2.3% and BJS's to 0.5%, so it is **not visible at admin scale**.

**Why it is explicitly NOT ruling-grade.** One DGP family, one panel length, **G = 40 only**,
**30 draws**, `cluster(i)` only, no pooled/coarse-`i` cross arm, no G sweep. The `i = cohort` arm
has no 100-draw denominator at all, because `mc.do` captured that arm from `e(b_es_att)` — which
`fastnose` does not post (§C.2.3): the help defect and the measurement gap are the same fact.

**Exactly what MC is needed before anyone rules.** The `se_mc.do` shape, run as a
`review_process` MC and **not** as a suite test: **500+ draws**, arms
{`esfixedbaseperiod`, pooled `esprelength(0)`} × {`i` = unit id, `i` = cohort id} × {common
onset, staggered} × **G ∈ {20, 40, 200}**, dumping `e(b)` and the `e(V)` diagonal each draw, with
**csdid `long2` and did_imputation as controls in every cell** (they give the same points to
1e-14/1e-7, so they isolate the variance path in every arm). Flag rule: any arm whose mean
reported SE / empirical SD falls below 0.90, or whose 95% coverage falls below 0.93, **while both
controls in the same draws sit inside those bands**. Add a coarse-cell arm to test the mechanism
directly. Cost: **~30 min** at 3.3 s per draw-pair.

**The mechanism on offer is a hypothesis, not a diagnosis.** `did_imputation`'s own help warns
that SE groups which are too small carry "downward bias from overfitting", and its default pools
by cohort-year; with `i` = the unit id wooldid's treated cells are singletons and the saturated
regression's treated-cell residuals are mechanically zero. That is consistent with the `i` =
cohort arm and `dta_s_admin` both being clean — and **it does not explain why the common-onset
unit-`i` arm is clean**. Do not repeat it as a finding.

**Status: number-mover, held for the author.** Nothing in this section should be treated as a
defect until the 500-draw MC exists.

**Adjacent, and deliberately not counted as a defect: honestdid does not refuse a degenerate
covariance.** On a synthetically corrupted covariance, `honestdid, delta(rm)` with a zeroed pre
row returns **rc 0 with no warning at all** and a silently different robust CI (lb 0.66881 at
M=0.5 against the clean run's 0.64679, same upper endpoint) — `staggered_refs/results.txt:198-203`
vs `:147`, console at `refs.log:17470-17560`. With the `-1e100` sentinel on the same row it
returns rc 0 but flags "(*) CI is open at both endpoints" with bounds ±2.7496; a zeroed **post**
row is silently accepted by both deltas. Treating a zero-variance coefficient as known exactly is
defensible for a genuine zero, so this is not a bug report against honestdid; the hazard is that
wooldid's structural zero means "not identified". **Consequence, which is a test-design rule:**
every honestdid handoff in the suite and in user code must assert `e(es_ndegen_att)` is **absent**
first — honestdid's return code is not a degeneracy guard. `wooldid.sthlp:1291` already requires
"a covariance with no structurally degenerate rows" for the handoff, so the help is not wrong; it
does not say that violating the condition **fails silently downstream**, and one clause there
would be worth adding. **Caveat: every probe was a synthetic corruption of a clean covariance.**
Producing a genuinely degenerate wooldid fit (`e(es_ndegen_att)` > 0) and feeding it through both
routes was not done (§C.7 iv).

---

### C.7 What could **not** be benchmarked, and why

The four probes' remaining-gaps lists, consolidated and grouped by cause. This is the honest
answer to "what is still unvalidated".

#### (i) Blocked by the external tool — nothing in this repo can close these

| gap | the block | what would close it |
|---|---|---|
| **Direct weighted jwdid comparison** | installed jwdid v1.41 returns r(101) for `[aw]` and `[pw]` (§C.4.1). All four weighted legs (S5, S6, D5, J6) run against an *emulation* of jwdid's own `reghdfe` line, not through `jwdid` + `estat` | a jwdid release whose replay guard (`jwdid.ado:41`) admits weights; the adopted script already prints `J6_JWDID_PW_RC` as the tripwire |
| **`boottest` p-value calibration** | all six measured wild-bootstrap p's (3 fixtures × equal-tail/symmetric/Webb, `reps(999)`) printed **0**, below the 1/999 resolution, because every available fixture's effect is **10–35 analytic SEs from zero**; wooldid's own comparators are equally saturated (8-state panel: percentile p = 0 at resolution 2/299, Wald p = 2.49e-06, analytic p = 2.0e-06) | a fixture built to sit near the 5% boundary — see (iv) |
| **xtevent's Sun-Abraham mode** (`cohort()`/`control_cohort()`) | requires `avar`, **not installed on this machine**; xtevent's own aggregation and its `overlay`/`proxy` machinery were untouched | installing `avar`, then one run |
| **`bootstrap(bayes)`** | **no Dirichlet mode exists anywhere**: Stata's prefix has none, and `bayesboot` (SSC) needs `iweight` and cannot drive wooldid | nothing external; only internal consistency checks are possible |

#### (ii) No external analogue exists — different by construction, not merely untested

- **Continuous treatment (`contreat()`) and the whole `ame` family**, including the `dta_cont`
  fixture: csdid and did_imputation have no continuous-treatment analogue; jwdid, ritest and the
  `bootstrap` prefix were never pointed at it. **No external oracle at all.**
- **`subgroups()`** and the pairwise contrasts: no counterpart in any of the four.
- **The adversarial line test** and `pretrendpower()`: no counterpart.
- **`att_it` / `att_i` / `att_itime`**, `impvalues`, `unconditionalse`, `customvce`,
  `oostreatedcontrols()`, `timetrends`, `esrelativeto(k ≠ −1)`: no jwdid, csdid or BJS counterpart.
- **`esprelength(k>0)` under the pooled reference**: it dummies the k nearest pre-periods and
  shrinks the pooled baseline; jwdid has no analogue of that geometry, so the pooled legs must
  stay `espost`-only and **pre horizons are benchmarked only under `esfixedbaseperiod`**. The same
  geometry is what breaks the BJS equivalence (§C.2.2).
- **Controls**: `controls()` / `coarsecohortcontrols()` / `interactivecontrols()`. csdid's DR
  machinery and wooldid's regression part company with covariates by construction; BJS matches
  only in the "fine cohorts" case; `coarsecohortcontrols()` has no analogue anywhere.
- **wooldid's percentile bootstrap p**: the Stata prefix has no percentile p at all, so only
  `bootdisplay(wald)` + `bootdf(normal)` is comparable. Likewise **wooldid's default and
  `permabsp` permutation p**, which carry an add-one that neither ritest nor `permute` has or
  offers (verified from source this pass: `ritest.ado:861-870` is `c/n` via `cii`;
  `permute.ado:1163-1168` `SetTwoSided` is `min(n, 2·min(c_lo,c_hi))/n`). Only `permrawp` maps.
- **BCa**: wooldid jackknifes by **cluster**, the prefix by **observation** unless
  `jackknifeopts(cluster())` — a genuine convention difference needing its own measurement.
- **`pretrends(k)` vs wooldid's pre coefficients** outside common onset: a *conditional* identity
  (1e-13 under common onset; 0.012–0.029 absolute on `costag`, 0.006–0.049 on `dta_binary`, i.e.
  at most **0.22 SE** and **0.46 SE** respectively). Do not write an equality test there. The
  joint pre-trend tests are different statistics outright (`e(joint_horiz)` 0.858 vs BJS `pre_p`
  0.850 on cohom, but 0.011 vs 0.556 on costag).
- **No-never-treated designs**: wooldid errors r459 while jwdid drops t ≥ max(g) and continues.
  Out of scope by construction — a leg there would test the divergence.
- **Multi-way clustering**: refused by wooldid's own bootstrap engine; not probed.
- **Poisson with `i` = unit or `icoh`**: jwdid's poisson path gives all never-treated **one**
  cohort FE, so only the `g2` cohort-FE formulation was compared (and there it is exact to solver
  precision: cells 6.0e-9 … 1.1e-8, aggregates 1.0e-9 … 2.1e-9, SEs 1.8e-12 … 9.3e-10).
- **wooldid's unit-level-`i` clustered V**: no **independent-engine** benchmark exists anywhere in
  this work. jwdid's V is INFO-only on every unit-level leg (SE ratios 0.62–1.24) and
  `test/suite/newfeat/reghdfe_oracle.do` remains the only check. *Partial exception, worth
  recording:* the bootstrap probe's hand-rebuilt single-absorb stage 1 plus a count-weighted
  `lincom` reproduced wooldid's `main_att` **and its clustered analytic SE** at
  reldif 6.07e-16 / **1.22e-14** on `dta_binary` at `i` = unit (`run6.log:17836`). That is a
  `reghdfe`-based reconstruction of the aggregation, covering `att` only and not the ES horizons
  — not an independent variance engine.

#### (iii) Simply not run yet, and cheap

| gap | cost to close |
|---|---|
| **xtevent `impute(stag)` on a staggered panel** — measured on ONE dataset (`cohom`, common onset, balanced); the mode's whole point is staggered adoption, and until this runs the Tier-1 xtevent claim is **common-onset only and an extrapolation elsewhere** | **0.03 s** per call; should lead the follow-up |
| **jwdid J11 on `test/suite/data/dta_s_biennial.dta`** — jwdid's reference is hard-wired to g − gap (gap = min spacing over the pooled unique g>0 and t values), so a biennial panel silently moves the reference. The highest-value single addition; may resolve to a **documented non-mapping** rather than a test | a few seconds inside the existing ~30 s script |
| **csdid with `notyet` controls** — csdid was run never-treated-only throughout; wooldid's control group is every not-yet-treated cohort *plus* never-treated, so `notyet` may be the better map on panels with no never-treated units | 0.18–0.53 s per call |
| **X6 convention identities** (recompute c_lo/c_hi/c_abs offline and pin every p convention) | **free** — no Stata run beyond one `permute()` run per convention |
| **A direct `permrawp` engine run on `dta_binary`** — today the external p's are compared to shares recomputed *offline* from `e(permdist_b)`; the engine-vs-offline identity is pinned only on the tiny panel | free inside X6 (one rerun at the same `rseed`) |
| **The full `e(b)` column set** — only the pooled `att` column was compared under permutation, and only the four-coefficient `att` vector under bootstrap. `e(permdist_b)` and `e(bootdist)` carry every coefficient (pre-treatment, ES horizons, subgroup rows, `gteffects` cells) | "nearly free" / "free in design, costs only runtime" |
| **Residual-saturation on `dta_binary`** — `predict double rsd, res` failed after `reghdfe`, so `SENT_P1_SAT_zeroresid_inwindow=.` and `..._base=325` (`run6.log:17757, 17769`) are **stale and must not be cited** | one re-run with `reghdfe …, resid` |
| **`permtstat`** — not benchmarked at all. ritest can rank `_b[main_att]/_se[main_att]`, but only off a per-draw command computing analytic SEs (no `fastnose`) and wooldid's studentized observed/replication statistic must be matched term for term | ~1.8 s/rep; feasible, unmeasured |
| **honestdid beyond `delta(rm)`/`delta(sd)` at method defaults** — FLCI, Conditional, C-F, `l_vec()` | seconds each, but `delta(rm)` at k=8 is 7.2–8.1 s and scales badly in k |
| **ritest `strict`, `eps()`, `left`, `right`; `permute`'s `standardize`** — read from source, never exercised, so the eps mapping is verified only at the default 1e-7 | minutes |
| **`bootdf()` values other than default and `normal`; `ties`** | minutes |
| **`dta_s_admin` under `permute()`** — untested; two route constraints known in advance: `permute(probabilities)` needs a numeric cluster variable, and the `permsave`→ritest replay keys on numeric cluster values through `onset_<clusterval>` column names, so a string-ish id fixture may need an `encode` first | one run |
| **A K>2-level probability table** — the cumulative-probability walk in ritest's `samplingprogram` must mirror `wooldid_permxp` exactly; today's generalisation is **reasoning in the report, not a measured result** (the one measured design is 6 clusters, two levels) | one run |
| **The unbalanced-panel identity boundary** — measured only at `dta_s_admin`; *where* the exact identity breaks (missing periods? unbalanced onset? multiple obs per i-t cell?) was never isolated, so Stage C has no principled predicate for choosing Tier 1 vs Tier 2 on an arbitrary dataset | a small sweep |
| **Data classes** — all benchmark data is `wellbehaved`. Not run: `dta_s_bigcodes` (the id-collision hazard this work hit by accident with FIPS codes), `dta_s_extmiss`, `dta_s_singleton`, `dta_s_biennial`; no stress or pathological class anywhere | one run each |
| **Inference variants** — multi-way clustering (jwdid side), `robust`, `fe()`, `level()` | minutes |
| **`reps(999)` permutation pair** — *all* measurement is at 199 reps (`dta_binary`) and 20–1000 (tiny panels); every `reps(999)` figure in the ritest report is a **linear extrapolation** from measured per-rep wall time on a loaded machine | **one unattended ~45-minute run** |
| **The S1/S2 coverage MC** (§C.6) | **~30 min** at 3.3 s per draw-pair |

#### (iv) Needs a fixture that does not exist yet

- **A failed or shrunk bootstrap replication.** Every run in the bootstrap probe had
  `N_misreps` = 0, so the documented divergence — the prefix **drops** a replication whose `e(b)`
  changes size (`e(N_misreps)`), wooldid **keeps** the row with missing cells (`e(boot_N_V)` per
  column) — is taken from code and help and has **never been measured**. Test C-BOOT-7(c) cannot
  be written until a fixture engineered to lose a horizon in some replication exists.
- **A deliberate `permsave()` column collision.** wooldid documents (sthlp:485) that an
  `onset_<clusterval>` column can be silently omitted on a name collision with **nothing inside
  the file marking the gap**. The probe checked only row/column counts (80 × 199+1). X3's
  proposed guard — "`onset_` column count == cluster count" — is therefore **untried**.
- **A near-5%-boundary fixture for `boottest`** (e.g. `dta_binary` with the effect scaled down
  and/or the noise raised, at G = 8–12). Without it no numeric pass rule for a p-value comparison
  can be calibrated; see (i).
- **A genuinely degenerate wooldid fit** with `e(es_ndegen_att)` > 0, fed through both honestdid
  routes. Every sentinel probe was a **synthetic corruption of a clean covariance**, so §C.6's
  adjacent finding establishes how honestdid behaves on those inputs, **not** that wooldid's own
  degenerate output reaches honestdid in that shape. That is the missing half.

---

### C.8 The tests Stage C should implement

Three classes, and the distinction matters because it sets the tolerance:

- **[I] deterministic identity** — two estimators, no RNG anywhere. Gate at machine precision
  (`atol 1e-12 + rtol 1e-10`), except where two ML/FE solvers meet, which floors at 1e-6/1e-8.
- **[S] same-seed exactness** — a common seed puts both engines **on the same RNG stream**, so
  the correct assertion is `rtol 1e-12`, **not** a Monte-Carlo band. At R=199 the MC band
  `SE/sqrt(2R)` is 0.0039 — three orders of magnitude too loose, and it would pass through real
  breakage.
- **[M] Monte-Carlo tolerance** — independent draw sets; compare one-sided tails at 4 MC standard
  errors, never a two-sided p.

**Write every Tier-1 tolerance in the mixed form `|a − b| ≤ atol + rtol·|b|`, never as a relative
difference alone**: pre-period coefficients on random draws pass through zero (costag pre7 =
0.000617 cleared a pure relative gate only by luck, and xtevent's worst case is a relative
1.1e-12 on an absolute 3.8e-15).

#### jwdid

1. **[I] `newfeat/jwdid_crosscheck.do` — ADOPT AS IS.** A ready-to-adopt **718-line** prototype
   exists at `stageC_feas/jwdid/jwdid_crosscheck.do`, self-contained (helpers inlined), arg 1 =
   path to `wooldid.ado`, pre-flight `cap which jwdid` printing `JWDID_XCHECK_SKIPPED` when jwdid
   is absent (as `reghdfe_oracle.do` does for ppmlhdfe), success line `JWDID_XCHECK_ALL_PASS`.
   Verified end to end: **exit 0, 25.9 s wall, 86 judged comparisons, 0 FAIL**, and all 97 `RES|`
   lines byte-identical to the previous run. **Two edits needed**: the default `wooldid.ado` path
   in the preamble, and the two absolute data paths (`test/data/dta_binary.dta`,
   `test/suite/data/dta_s_admin.dta`). Legs **J1–J10**: the synthetic common-onset DGP (40×12,
   onset 7, seed 20260824), `dta_binary`, a derived unbalanced panel, and `dta_s_admin`; wooldid
   `y icoh t ttre, att cluster(i) [esfixedbaseperiod] esprelength(P) espostlength(Q) gteffects`
   with `icoh = cond(ttre<., ttre, i)` against `jwdid y, ivar(icoh) tvar(t) gvar(g) [never]
   cluster(i)` plus `estat simple/event/group/calendar`. **Pass rule:** `e(N)` equal; cells b
   ≤ 1e-8 and cells se ≤ 1e-6; att b ≤ 1e-8, att se ≤ 1e-6; ES b ≤ 1e-8, ES V ≤ 1e-6; group and
   calendar b ≤ 1e-8 (mreldif = max|w−j|/(|j|+1)). Poisson (J7) at ≤ 1e-6 both sides. **J5 is the
   unit-level-`i` leg and judges V as INFO, never as a verdict.** Measured worst per leg: J1
   1.9e-16 … J10 1.5e-13 — four orders of headroom. Budget **~30 s**.
2. **[I] J11, `dta_s_biennial` — not yet written or run.** The J1 mapping on a biennial panel, at
   the same tolerances, to probe jwdid's `gap` logic. Expect either agreement (the mapping is
   wider than documented) or a clean **documented non-mapping** that belongs in the notes rather
   than the suite.

#### Staggered-DiD references

3. **[I] `newfeat/csdid_identity.do` — the single highest-value new test.** Datasets: `cohom`,
   `cohet`, `costag` (generator `staggered_refs/refs.do:91-113`) and `test/data/dta_binary.dta`.
   wooldid `y i t ttre, att cluster(i) esprelength(P) espostlength(Q) esfixedbaseperiod` with
   (P,Q) = (5,6)/(5,6)/(6,8)/(7,7); external `csdid y, ivar(i) time(t) gvar(g) long2` with
   `g = cond(ttre<., ttre, 0)` then `estat event` / `estat simple`; map `Tm<h>`→`pre<h>`,
   `Tp<h>`→`post<h>`, simple ATT→`main_att`. **Pass:** every mapped horizon
   `|b_w − b_cs| ≤ 1e-12 + 1e-10·|b_cs|`; ATT `≤ 1e-12 + 1e-12·|b_cs|`. Sentinel
   `CSDID_IDENTITY_ALL_PASS`. Measured worst **1.02e-14** over 100 MC draws × 2 DGPs. ~3 s/dataset.
4. **[I] `newfeat/bjs_identity.do`.** Same four datasets; wooldid `… att cluster(i)
   esprelength(0) espostlength(H+1)` against `did_imputation y i t ttre, horizons(0/H) minn(0)
   cluster(i)`; map `tau<h>`→`ES_att_post<h>`, BJS aggregate `tau`→`main_att`. **Pass:**
   `|b_w − b_bjs| ≤ 1e-6 + 1e-6·|b_bjs|` on every post horizon and the ATT (measured worst
   1.16e-07 — the two solvers' FE tolerance, tightenable with `tol()` on both sides).
   **Negative control in the same file**, guarding the §C.2.2 condition against silent drift:
   rerun with `esprelength(4)` and assert `max_h |b_w4 − b_bjs| > 1e-3`, recording the shift.
   **Third leg:** `anticipation(1)` vs `did_imputation … shift(1)` at the same 1e-6 gate, and vs
   wooldid on `ttre−1` at **exact 0, SEs included**.
5. **[I] `newfeat/xtevent_identity.do`.** `cohom` **first, then `costag`** — the staggered panel
   is the gap (§C.7 iii). `gen byte z = (ttre<. & t>=ttre)` then `xtevent y, policyvar(z)
   panelvar(i) timevar(t) window(-P Q) impute(stag) vce(cluster i)` against the same
   `esfixedbaseperiod` wooldid fit; map `_k_eq_m<h>`→`pre<h>`, `_k_eq_p<h>`→`post<h>`, and
   **exclude the outermost lead/lag on each side** — the window endpoints are binned and are not
   wooldid horizons. **Pass:** `|b_w − b_xt| ≤ 1e-12 + 1e-10·|b_xt|` (measured worst 3.8e-15
   absolute, 1.1e-12 relative on a 0.0036 coefficient — **the relative gate is the binding one**).
   Skip with a sentinel if `which xtevent` fails. 0.03–0.24 s.
6. **[I] `newfeat/honestdid_routes.do`** (or fold into `honestdid_grid.do`). For 2–3 existing grid
   cases, fit once with `esfixedbaseperiod` and the default `esrelativeto(-1)`, then call
   honestdid twice: route 1 ``honestdid, numpre(`=e(numpre_att)') b(b) vcov(V) delta(`d')
   mvec(0.5(0.5)2) parallel(0)`` and route 2 ``honestdid, pre(<ES_att_pre* indices>)
   post(<ES_att_post* indices>) delta(`d') mvec(0.5(0.5)2) parallel(0)`` on the posted `e(b)`/`e(V)`
   (find the indices by scanning `colnames e(b)` for the `ES_att_` prefix, excluding
   `ES_att_it_`). **Pass:** `mreldif` over CI columns 2–3 ≤ **1e-12** for `delta(rm)` (measured 0)
   and ≤ **1e-7** for `delta(sd)` (measured 2.71e-09, honestdid's own FLCI grid search); assert
   the `omit` variant changes nothing (measured 0). **`parallel(0)` is mandatory.** **Cheap leg
   that covers PS-62 without paying honestdid's 8 s:** assert `mreldif(e(V) ES block,
   e(V_es_att)) ≤ 1e-15` and `mreldif(b) == 0` directly — measured 2.73e-17 and 0, **0.1 s**.
7. **[I] `newfeat/honestdid_degenerate.do`.** One clean `esfixedbaseperiod` fit; corrupt a single
   ES row **and** column of the covariance three ways (diagonal → `-1e100`, → 0, → missing) on a
   pre row and on a post row, and pin the measured table (honestdid 1.1.2): `-1e100` pre → `rm`
   rc 0 with the "(*) CI is open at both endpoints" flag and CI = the grid bounds ±2.7496, `sd`
   rc 3598; `-1e100` post → `rm` rc 198 (ECOS), `sd` rc 0 with a CI identical to clean; 0 on a pre
   row → `rm` **rc 0 with no warning** and lb 0.6688 against the clean 0.6468 at M=0.5, `sd`
   rc 3498; 0 on a post row → both deltas rc 0 and identical to clean; missing on either row →
   rc 198 / 3598 / 3498. **The real assertion is the doctrine: rc 0 is not a safety signal**, so
   every honestdid handoff must first assert `e(es_ndegen_att)` is absent.

#### Permutation (ritest + `permute`)

8. **[I] X1 — `permute, enumerate` vs `permexact`.** Fixture: synthetic tiny panel, 48 obs,
   6 clusters, onsets {4,4,4,.,.,.}, 8 periods (built by `ritest/runA.do`). wooldid
   `y i t ttre, att cluster(i) permute(shuffle) permexact`. External: a cluster-level frame (one
   row per `i`, never-treated recoded to 0) plus an **eclass `sortpreserve`** wrapper `wd_cl`
   that `frlink`s the onset back onto the panel frame and calls `wooldid y i t ttre_p, att
   cluster(i) fastnose`, driven by `permute ons _b[main_att], enumerate nodrop nodots nowarning
   saving(f, replace double): wd_cl ons`. **Pass:** sorted `_pm_1` vs sorted `e(permdist_b)` at
   reldif ≤ 1e-10 (report the bitwise count; 20/20 measured, max reldif 0); `r(n_reps)` ==
   `e(perm_npossible)` == `e(perm_completed)` == 20; `r(c_lower)/20` == `e(perm_p_lo)`,
   `r(c_upper)/20` == `e(perm_p_hi)`, `r(p_twosided)` == `e(perm_p_att)`, each to 1e-12;
   `|r(b) − e(perm_comparator_b)|/|e(perm_comparator_b)| ≤ 1e-13` — **never bitwise**. ~2.2 s.
9. **[S] X2 — sampled assignment space.** Same panel and wrapper, `reps(500) rseed(7)`. **Pass:**
   every draw within 1e-10 relative of some member of the exact 20-value set (500/500 bitwise
   measured); exactly 20 distinct cluster-level patterns; every pattern carries 3 onset-4 and
   3 never (0 multiset violations); min pattern frequency ≥ 1 and chi2(19) < 43.8; tails compared
   against exact-set counts **recomputed under `permute`'s own eps rule** (`#{s ≤ s_obs+1e-7}` =
   20, `#{s ≥ s_obs−1e-7}` = 2), and **any tail whose exact share is 0 or 1 asserted for equality
   rather than z-tested** (the lower share here is exactly 1, so a `4·sqrt(p(1−p)/n)` rule would
   be a zero-tolerance test). Informative upper tail .04 vs exact .05, |z| < 4 at MC se .0097.
   ~44 s.
10. **[S] X3 — ritest replay of wooldid's own draws. The highest-value permutation test.**
    Fixture `dta_binary` with the **null** outcome `y0 = y − eff*(t>=ttre)`, reps 199. wooldid
    `y0 i t ttre, att cluster(i) permute(shuffle) permreps(199) rseed(#) permsave(f) replace` —
    **note `permsave` takes a PLAIN name (§C.3.3)**. Transpose `f.dta` (`keep _rep onset_*`;
    `xpose, clear varname`; `i = real(subinstr(_varname,"onset_","",.))`; drop the `_rep` row and
    the rep-0 column; rename to `ttre1..ttre199`) and **assert the `onset_` column count equals
    the cluster count before transposing**. External: `ritest ttre _b[main_att], reps(199)
    samplingsourcefile(rsB.dta) samplingmatchvar(i) saving(h, replace double) nodots: wooldid y0
    i t ttre, att cluster(i) fastnose`. **Pass:** rep-by-rep `h._pm_1` vs column `main_att` of
    `e(permdist_b)` at reldif ≤ 1e-10 (report the bitwise count; **199/199, max reldif 0**
    measured); `sd(h)` == `e(perm_sd)` to 1e-10; `r(p)` == `c_abs/n` recomputed from wooldid's own
    draws to 1e-12; `r(b)` == `e(perm_comparator_b)` at reldif 1e-13; and **assert no draw lies
    within 1e-7 of |s_obs|** (0 of 199 measured) — if one does, compare against
    `#{|s| ≥ |s_obs| − 1e-7}/n` instead. ~9 min at 199 reps, ~4.5 min at 99.
11. **[M] X4 — independent-draw p agreement.** Same null outcome. wooldid `permute(shuffle)
    permreps(199) rseed(11)` with c_lo/c_hi/c_abs recomputed offline; external A
    `ritest ttre _b[main_att], reps(199) cluster(i) seed(22) saving(f, replace double) nodots:
    wooldid y0 i t ttre, att cluster(i) fastnose`; external B the `wd_cl` wrapper under
    `permute … reps(199) rseed(33)`. **Pass:** `|r(p)_ritest − c_abs/n| ≤ 4·se_diff` with
    `se_diff = sqrt(p̄(1−p̄)(1/n1+1/n2))` (z = +0.59 measured); the same for each `permute` tail
    (z = ±1.40); `|sd_ext − e(perm_sd)| ≤ 4·sqrt(sd1²/(2(n1−1)) + sd2²/(2(n2−1)))` (z = −0.41,
    −1.04). **Compare one-sided tails, never the two-sided p** (min-and-double amplifies a tail
    difference near ½), **and use the null outcome** — on the effect-carrying outcome every engine
    sits at its floor (wooldid .04 = 2/(n+1) vs ritest 0 at 49 reps) and the test is empty.
    ~10 min.
12. **[M] X5 — strata and probabilities.** *Strata*: 8-cluster/2-strata tiny panel (strata
    {4,4,.,.} | {4,.,.,.}, 24 assignments); wooldid `… permute(shuffle) permstrata(s) permexact`
    against `ritest ttre _b[main_att], reps(720) cluster(i) strata(s) saving(f, replace double)
    nodots: wooldid … fastnose`. **Pass:** 24 distinct patterns; per-stratum level counts fixed on
    every draw (0 violations); membership in the exact 24-value set at reldif 1e-10; `r(p)` within
    4 MC se of `c_abs/24` = .0417 (measured .0375, z = −0.56). *Probabilities*: 6-cluster panel,
    per-cluster onset probabilities .8/.6/.4/.2/0/0; wooldid `permute(probabilities) permprobs(…)
    permreps(400) permrawp` against `ritest … samplingprogram(permprob_sp)
    samplingprogramoptions("clus(i) lev(4) plist(1 .8 2 .6 3 .4 4 .2 5 0 6 0)")`. **Pass:** both
    p's within 4 MC se of the exact design tail from the 15-point estimable support (.4093;
    measured .3842 z = −1.00 and .3802 z = −1.16); per-cluster treated frequencies within 5
    binomial sd of each cluster's conditional marginal; probability-0 clusters never treated;
    ritest's dropped-draw count equals its all-never count (the role `e(perm_failed)` plays for
    wooldid). ~90 s + ~45 s.
13. **[I] X6 — convention identities; free.** From one `wooldid … permute(shuffle) permreps(R)
    rseed(#)` run on any fixture, recompute c_lo, c_hi, c_abs offline from `e(permdist_b)` and
    `e(perm_comparator_b)` and assert to **1e-12**: `e(perm_p_att)` == `min(1,
    2(1+min(c_lo,c_hi))/(R+1))`; a rerun at the same `rseed` under `permabsp` gives
    `(1+c_abs)/(R+1)` **and a bitwise-identical `e(permdist_b)`**; under `permrawp`, exactly
    `c_abs/R` — the number ritest and `permute` report; under `permside(left)`/`(right)`,
    `e(perm_p_lo)`/`e(perm_p_hi)`; under `permexact` every one of these loses the add-one. **This
    is the block that makes X3 and X4 compare like with like**, and it closes the one provenance
    gap in the permutation study for free.

#### Bootstrap

14. **[S] C-BOOT-1 — core oracle, count leg.** `test/data/dta_binary.dta`. A: `wooldid y i t ttre,
    att cluster(i) esprelength(1) espostlength(1) bootstrap(standard) reps(199) rseed(20260912)`.
    B: `bootstrap _b, reps(199) seed(20260912) cluster(i) idcluster(ni) nowarn saving(pf, double
    every(1) replace): wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
    fastnose`. **Pass at max reldif ≤ 1e-12** on (i) every cell of `e(bootdist)` vs `pf.dta`'s
    `_b_*` columns (199×4 = 796 cells), (ii) the bootstrap SE vs `e(se)`, (iii)
    `e(boot_ci_percentile)` vs `estat bootstrap, all`'s `e(ci_percentile)`, (iv) `e(boot_ci_bc)`
    vs `e(ci_bc)`, (v) `e(boot_z0)` vs `e(z0)`. **Must NOT assert `e(boot_ci_normal)` here** — the
    default `bootdf(clusters)` = t(79) against the prefix's z differs by 2.5e-3 **by design**.
    Measured 2.43e-16 / 8.10e-15 / 2.47e-16 / 2.88e-15 / 2.41e-15 / exactly 0. ~128 s + 145 s
    (timed at R=200; R=199 on this fixture costs the same).
15. **[S] C-BOOT-2 — convention pin, Wald leg.** The same pair with `bootdf(normal)` added to the
    wooldid call. **Pass:** max reldif ≤ 1e-12 between `e(boot_ci_normal)` and the prefix's
    `e(ci_normal)`. Companion guard in the same file: **without** `bootdf(normal)`, assert the
    wooldid half-width ÷ prefix half-width == `invttail(79,.025)/invnormal(.975)` to rtol 1e-10 —
    i.e. the gap is exactly the critical-value ratio and nothing else.
16. **[I] C-BOOT-3 — quantile rule, no second estimator; the cheapest guard in the set.** Take any
    bootstrap run's `e(bootdist)`, write it one row per replication, then `bstat, stat(b_obs)
    n(N)` and `estat bootstrap, all`. **Pass:** `e(boot_ci_percentile)` vs `e(ci_percentile)`,
    `e(boot_ci_bc)` vs `e(ci_bc)`, `e(boot_z0)` vs `e(z0)` at max reldif ≤ 1e-12 (**measured
    exactly 0**). Run at **both R=199 and R=200**: only an R divisible by 40 at `level(95)`
    reaches the integer-`n·alpha` branch that averages adjacent order statistics (measured
    1.39e-15 at R=200); 199/499/999 exercise only the `ceil` branch.
17. **[S] C-BOOT-4 — re-identification / `boottstat` oracle; highest value per second.** Fixture
    `use test/data/dta_binary.dta; keep if inlist(i,1,2,3,26,27,28) | inrange(i,51,56)` (G=12,
    N=144, att = 1.00571, se_obs = .092284). A: `wooldid y i t ttre, att cluster(i)
    esprelength(1) espostlength(1) bootstrap(standard) reps(200) rseed(31415) boottstat
    bootsave(w12ts) replace`. B: `bootstrap b=_b[main_att] se=_se[main_att], reps(200) seed(31415)
    cluster(i) idcluster(ni) nowarn nodots saving(r12fresh, double every(1) replace): wooldid y ni
    t ttre, att cluster(ni) esprelength(1) espostlength(1)`. Merge `w12ts.dta` rows `_rep` 1..200
    (columns `main_att`, `se_main_att`) against `r12fresh.dta` rows 1..200 (`b`, `se`). **Pass:**
    max reldif ≤ 1e-12 on **both** columns — **measured exactly 0.000e+00 on both, 200/200
    matched.** This is the only external check of the re-id leg that exists and it simultaneously
    certifies `bootsave()`'s layout and its `se_<name>` column naming. ~46 s + 52 s.
18. **[S] C-BOOT-5 — count leg == physical resample for the observation-weighted estimands.** Same
    12-cluster fixture and seed; drive the prefix twice, once as `wooldid y i t ttre, att
    cluster(i)` (copies pooled into the original cohort) and once as `wooldid y ni t ttre, att
    cluster(ni)` (each copy its own cohort). **Pass:** max reldif of the per-replication **point
    estimates** ≤ 1e-12 (measured 2.14e-15 over 200 draws; 6.82e-15 over 199×4 on `dta_binary`).
    This pins the design claim that justifies the cheap count leg for `att`/`ame`; no existing
    test would notice if a change to the margins weights broke it. **Do not extend the rule to the
    `se` column** — the two legs differ there by design (that is test 17's territory). Site it on
    the 12-cluster fixture (~97 s), not `dta_binary` (~640 s at R=199).
19. **[S] C-BOOT-6 — second dataset, realistic panel.** `test/suite/data/dta_s_admin.dta` (29,667
    rows, 4 cohorts, 48 states, `e(N)` = 29,067 after 600 drop). A: `wooldid y cohid year ttre,
    att cluster(state) esprelength(1) espostlength(1) bootstrap(standard) reps(99) rseed(4242)`.
    B: the prefix at `seed(4242) cluster(state) idcluster(ns)` over the same call with `fastnose`.
    Same five assertions as test 14. **Pass at max reldif ≤ 1e-12 — not 1e-15**: the 30× larger
    panel loses a decade to summation (measured `e(b)` exactly 0, draws 4.26e-14, SE 5.84e-15,
    percentile 2.82e-14, BC 3.13e-14, z0 exactly 0). Covers a panel where **the cluster variable
    is not the `i` variable**, cluster sizes are unequal, and rows leave the sample. ~52 s + 59 s.
20. **[negative control] C-BOOT-7 — scope, three asserts, no second engine.** (a) `bootstrap _b,
    reps(20) cluster(i): wooldid y i t ttre [aw=w], att cluster(i)` must fail nonzero (the prefix
    refuses weighted commands) while the same wooldid call with `bootstrap(standard) reps(20)`
    must exit 0. (b) On `dta_s_admin`, `wooldid y cohid year ttre, att cluster(state)
    bootstrap(standard) reps(5) rseed(1) boottstat` must exit **exactly r(459)** with error text
    containing "requires the cluster() variable to nest i" (measured). (c) On a fixture engineered
    so some replication loses a horizon, assert the prefix reports `e(N_misreps)` > 0 while
    wooldid reports `e(boot_N_V)` < `reps()` **with the row retained** — the two engines are
    *expected* to disagree there. **(c) still needs a fixture that does not exist (§C.7 iv).**
    These are scope statements that stop a future agent "fixing" a deliberate difference.
21. **[I + recorded] C-BOOT-8 — `boottest` / hand-rebuild identity, two halves with different
    status.** *Asserted half*: on `dta_s_admin` and again on `keep if state <= 8`, build
    `gen relidx = year - ttre + 2 if ttre<. & year >= ttre-1; gen long cell = 0; replace cell =
    cohid*100 + relidx if relidx<.`, fit `reghdfe y i.cell i.year, absorb(cohid) cluster(state)
    tol(1e-9)` (**match wooldid's own tol** — the probe's rebuild omitted it, which is the likely
    reason the admin SE agrees only to 8.4e-12), then `lincom` the post cells
    (`mod(cell,100) >= 2`) weighted by their in-sample observation counts ÷ the total. **Pass:**
    the lincom estimate and its SE vs `_b[main_att]` / `_se[main_att]` at max reldif ≤ **1e-10**
    (measured 6.26e-14 / 8.40e-12 at G=48; ~1.4e-13 / 5.9e-13 at G=8; 6.07e-16 / 1.22e-14 on
    `dta_binary`). **This is the only external check of wooldid's analytic clustered SE in the
    whole study.** *Recorded-not-asserted half*: run `boottest ((<that lincomb>) = 0),
    cluster(state) reps(999) seed(11) nograph ptype(equaltail)` and log `r(p)`, `r(t)`, `r(CI)`
    beside wooldid's `e(boot_p_percentile)` / `e(boot_ci_percentile)`. **Assert only** that
    boottest exits 0 and that `r(t)` equals wooldid's analytic t at rtol 1e-6 (measured 14.446776
    / 36.677022 / 13.963961 against 14.446776159 / 36.677021625 / 13.963961283). **Do not assert
    the p** until a near-boundary fixture exists. Under 1 s on top of one wooldid run.

#### Not a suite test

22. **[M] `review_process` coverage MC for the S1/S2 candidates.** Template `se_mc.do`; full
    specification and flag rule in §C.6. ~30 min. **Do not site this in `test/suite/`.**

**Two numerical traps that must be written into every permutation and bootstrap test.**
(a) `saving()` is **float** unless `double` is passed — in *both* ritest and `permute`, and in the
`bootstrap` prefix. Without it, bitwise membership in the exact reference set falls from 1000/1000
to **0/1000** with a 5.55e-08 spread. (b) Both permutation prefixes report the observed statistic
`r(b)` through a **~15-significant-digit decimal round trip** — both re-read this pass:
``ritest.ado:607-608`` is ``local x = `name'[1]`` then ``char `name'[permute] `x'``, and
`permute.ado` carries the identical pair at ``:800``/``:802`` with the read-back at ``:942``
(``local value : char `name'[permute]``) — so `r(b)` sits one ulp off
`e(perm_comparator_b)` (reldif ~1.04e-15). Compare observed statistics at **reldif ≤ 1e-13, never
bitwise** — on one fixture the round trip happened to be exact, so a bitwise assertion would pass
on some fixtures and fail on others: the worst kind of flake. The same round trip bit the
*analyzer* during the study (a superseded `analyzeD.out` block prints c_lo = 19 and z = +5.13;
`analyzeD_D2.out` and `postcheck.out` supersede it with c_lo = 20).

---

### C.9 Citation corrections, and two source diagnoses

Every wooldid line number in the four probe reports was re-opened and read. **Six** citations
needed correction (rows 1–6); the last two rows are **not** corrections but source pins that
upgrade a measured-only claim to a diagnosed one. **No claim was found to rest on a line that does
not say what the probe said it says.**

| probe's citation | what the file actually has | effect on the claim |
|---|---|---|
| `sthlp:484` described as the **`permtstat`** instruction | line 484 is inside the **`permsave()`** option entry, which mentions `permtstat` in the same sentence | none — the line and the instruction are as quoted |
| `sthlp:284` described as the **`permute()` fast-mode** paragraph | lines **278-286 are the `fastnose` option entry**; the `permute()` fast-mode passage is at **437-438** | **strengthens** C.2.3: the `fastnose` entry itself contradicts sthlp:1295 about `fastnose` |
| `ado:5696-5700` and `:8513-8517` for the name check | the guard proper is **`:5697-5700`** and **`:8513-8516`**; `:5696`/`:8512` are the trailing-`.dta` strip | none; and the `bootsave()` twins the probe named without lines are **`:6042-6045`** and **`:6288-6291`** |
| D2 "presumably `error 459` surfacing through both the `wooldid` shell and `wooldid_worker`" | the **shell is correct** (`exit` at `:100`, and its own comment at `:61` says why). The reprint is **`:4770 error _rc`** after `cap noi wooldid_bootctl` at `:4769`; the identical pattern is at **`:4756`** for `wooldid_permctl`, and those are the only 2 of 14 `cap noi wooldid_*` sites that re-raise with `error` | **diagnoses** the defect and names a one-line fix at two sites |
| suite `parallel(0)` fix "uncommitted in the working tree" | landed in commit **`4280f3b`**, an ancestor of HEAD `21db885`; absent at `08ec0b4` | closed, **not counted** as a defect |
| `permute.ado:955` for the lower-tail count | the `gen double count = (name <= value + eps)` is at **`:954`**; `:967` (upper tail) is correct | none |
| "jwdid help documents `[pw]`, installed code returns r(101)" (measured only) | `jwdid.sthlp:16` and `:35` document it; **`jwdid.ado:41`**, the replay guard `syntax [anything(everything)], [example *]`, is parsed first and carries no weight spec, while the real syntax at **`:53` does declare `[pw]`** and `:202-203` passes the weight through | **diagnoses** it as a jwdid code defect (replay syntax shadowing), not a help error |
| honestdid "defaults to `parallel(4)` when `parallel` is installed" (measured only) | **`honestdid.ado:128-130`**: `cap parallel numprocessors` then `local parallel = cond(_rc, 0, 4)`; the author's own `* TODO: xx Should default be no parallel?` at `:132` | pins the default from source |

**Numbers this section deliberately does not cite, and why.**
`SENT_P1_SAT_zeroresid_inwindow` / `..._base` (`run6.log:17757, 17769`) — **stale**, the `predict`
that should have set them failed. The `analyzeD.out` D2 block — **superseded**. All `reps(999)`
costs — **linear extrapolations**, no 999-rep run was ever executed. "~7 minutes for 100
`boottstat` replications at G=80" — a **projection**; that run was never completed. The `2 − 1/G`
replication-SE inflation constant — the meat-matrix prediction is 1.92–1.99 and the **measurement
is 1.37–1.48**, because a resample holds only ~63% of the distinct clusters, so quote the
direction and rough size, never the constant. The `jwdid` probe return's `measured_agreement`
field — a self-referential placeholder; the numbers here come from the report's tables and its
`bench_run2.log` / `jwdid_crosscheck.log` citations instead.


---


# Part 9 - every bug fixed during this assessment, and what was left open


The chunk asks me to "end by reporting on all bugs fixed". This is that list, kept as the work
goes so nothing is reconstructed from memory at the end. It is split by what was broken, because
the split is itself the headline.

### 0. The headline, in two parts

**Stage A needed no change to the package at all.** `git diff --stat HEAD -- wooldidinstall/` was
empty across the whole Stage A commit. Every Stage A defect was in the ASSESSMENT SUITE - the
thing that is supposed to catch defects - or in documentation.

**Stage B then found no regression either.** Over 538 configs, 524 were explained before triage
and the other 14 were all configs meeting v1.3 for the first time. Not one PRE-EXISTING config
moved: strip the new ones out and the other 507 reproduce their prior verdicts exactly, 374
explained / 25 matching / 108 by return code. Adjudicating the 14 gave **10 old-bug-fixed, 2
intended-change, 2 mixed, 0 both-builds-defect, 0 unclear, and 0 regressions** - on every one of
them the current build is the more correct or the author-ruled side.

**And the closing run, over a larger universe and a changed package, found nothing further.**
After the author's 2026-09-13 rulings landed, the suite was re-run end to end over **559** configs
(21 more than the adjudication run) against the same v1.3 old side. **Zero DIFFERS at first
judging** - `EXPLAINED=411 MATCH=26 RC_EXPLAINED=122`, `unexplained=0 stale_ledger=0 unscoped=0`
over 559 of 559, with cert 43/43, newfeat 13/13, coverage OK and the degeneracy component clean on
its first pass over all 21 new configs. Details and the one thing that run did surface - 46 stale
ledger rows, 7 of which had blamed the wrong cause - are in Part 7.5.

The four package fixes below were found by the assessment, not by the comparison: none of them
moves a number, which is why the dual build could never have caught them.

### 0b. The original Stage A headline (kept for the record)

`git diff --stat HEAD -- wooldidinstall/` is empty across the whole Stage A commit. Every defect
below was in the ASSESSMENT SUITE - the thing that is supposed to catch defects - or in
documentation. A seven-area audit of the suite, each area re-checked by an adversarial critic,
found no code defect in wooldid that Stage A had to fix before the comparison could run. Whether
that survives Stage B is the open question; this line gets updated, not deleted, if it does not.

### 1. Defects that made the suite BLIND (it reported success it had not earned)

| # | what was broken | consequence | fix |
|---|---|---|---|
| 1.1 | Batch 3 re-posted `joint_horiz*`/`joint_Choriz*` as SCALARS under their old MACRO names, and the comparer keys its kinds separately | 361 values on 103 configs were blessed as "a removal plus an informational addition" and never value-compared at all | the comparer compares across the kind change; six blanket ledger rows retired; MATCH became reachable again on 25 configs |
| 1.2 | CORE-10's uniform `wooldid ` prefix on `e(cmdline)` was absorbed by an unbounded all-config ledger row | that row also hid `B21b-ic8-p`, which had only ever fired on the prefix | the comparer normalises the prefix; both rows gone |
| 1.3 | `coverage_check.py` counted an option as covered when its only clean-rc config was one the OLD build refuses | which is every `permute()` and `bootstrap()` option: 30 options credited as covered while the dual build only ever judged their return-code pair | a new COVERAGE_WARN rc-only class that names them, and says where their bit-exact tracking actually comes from |
| 1.4 | no component asserted that a config's target surface was LIVE | a surface dead on BOTH builds compares equal and passes, and the static coverage check still credits the option | the new `degeneracy` component, with a curated expectation file where every accepted degeneracy carries a citation |
| 1.5 | `run_full.sh`'s data-regeneration threshold was a bare count | a stale threshold makes the data step skip, and every config on a missing dataset dies r(601) on BOTH builds, which reads as a parity pass | threshold raised with the two new fixtures, and the comment now states the real builder split (12 + 12 + 1 = 25) rather than a miscount |

### 2. Defects that made a measurement MEASURE THE WRONG THING

| # | what was broken | consequence | fix |
|---|---|---|---|
| 2.1 | `mc_ri2.do` never passed `permstepdown`, which PS-45 made opt-in | the entire family-wise error phase wrote missing values and the analyzer printed an empty table | the option is passed; the analyzer now fails loudly on absent columns instead of reporting n=0 |
| 2.2 | `mc_ri.do` and `analyze_mc.py` were preregistered for the `|s|` convention now called `permabsp`, but the shipped default is equal-tailed | check R6 could not pass at any seed - a structurally impossible test | the bands certify the default users actually get, with the deviation recorded in the preregistration; the convention question is item 3 of the rulings document |
| 2.3 | `mc_bootse2.do` read only the undisplayed coefficient-draw p-values on its studentized arms | the displayed p was never checked | it now reads the displayed one too |

### 3. Defects that would have ABORTED or CORRUPTED a run

| # | what was broken | consequence | fix |
|---|---|---|---|
| 3.1 | a new sup-t reconstruction leg in `test_esband.do` wrote a one-line Mata `if (c) { x++; continue }` | Mata rejects it with "'}' found where nothing expected", r(3000), killing the entire certification suite | split across lines, and every `test_*.do` and newfeat script swept for the same shape (no other instance) |
| 3.2 | all three suite `honestdid` call sites omitted `parallel(0)` while the `parallel` package is installed | every `delta(rm)` call silently spawned four extra Stata processes, so the suite breached its own instance cap from the inside | `parallel(0)` added at all three sites; the same omission in the help's own example is item 6 of the external-benchmark defects |

### 4. Broken links: tracked things pointing at things that were not there

| # | what was broken | fix |
|---|---|---|
| 4.1 | `run_verify2.sh` was a tracked runner whose target `smoke2/verify2.do` was untracked, so a fresh checkout had nothing to run | both `smoke2/*.do` tracked; the runner keeps a guard that fails loudly if the target goes missing again |
| 4.2 | three permutation artifacts cited by the TRACKED F058 conclusion documents as their evidence of record were themselves untracked | all three tracked |
| 4.3 | the three new Monte Carlo drivers had no runner and no documented invocation | one shared cap-honouring runner `run_stagec.sh`, a usage header in each driver, and an inventory entry that says in plain words that they are scaffolds with no campaign and no analyzer |
| 4.4 | artifact classes regenerated by every certification run were being offered for staging | `test/.gitignore` covers them: `*.smcl`, `*.gph`, `perm3_*`, `l6probe_*`, and the shard scheduling and log scratch under the tracked `out2` tree |

### 5. Records that had gone stale or self-contradictory

- ledger reasons still describing the AMESE R4 refit removed on 2026-09-08; the SI-7 rows
  describing a horizon drop the default path no longer takes; PS-62's config count (360 -> 335);
  fourteen rows whose unquoted commas split their own reason field across columns.
- `DEGENERACY.md` cited a config count that could not be reconciled with the file it documents;
  it now states the shipped composition and the rule that you count UNRULED citations, not the
  HIDDEN class. The same off-by-one was in my own working notes and in the brief I gave an agent.

### 6. Still open, restated at the close of 2026-09-14

The nine issues the external-benchmark probes found are no longer the open list. Six are fixed,
one was never ours, and two are closed by ruling rather than by repair. What is actually open is
below, and it is shorter and sharper than what it replaced.

**Where these stand after the 2026-09-14 rulings.** Two items that stood here on the morning of
2026-09-14 - the subgroup pre-treatment over-drop and its invisibility - were ruled and are now
built; they moved to section 9b. Four remain below, and they are not all the same kind of thing.
Items 1 and 4 have been RULED and are kept here anyway, because a ruling to take no action does
not make the finding go away and a reader looking for known limitations deserves to meet them.
Item 2 is a live call still awaiting the author - a disclosure surface, not a computation.
Item 3 is not ours to fix.

1. **The bootstrapped adversarial statistic is INERT on a containment geometry - measured both
   ways.** Under the null `AB12/bt` rejected 0 of 200 at 5% and 0 of 200 at 10%, which alongside
   the analytic arm's 72.5% over-rejection looked like the bootstrap repairing it. The author
   ruled "measure power first, then decide", and power is now measured: with a genuine violation
   planted (a step, which is not a line through the reference period and so is outside BS-12's
   null), **the analytic arm rejects 35 of 40 with p reaching 1.7e-09 while the bootstrapped arm
   rejects 0 of 40, minimum p 0.080** - paired on the same 40 datasets, never once agreeing on a
   rejection. So it is not conservative, it is inert: it rejects nothing under the null and
   nothing under a violation the analytic arm detects at one part in a billion. A user reading a
   non-significant `AB12/bt` p learns nothing from it.
   **RULED 2026-09-14, after the measurement: LEAVE IT UNTOUCHED.** The options put to the
   author were to withhold the statistic under the existing containment screens, to withhold and
   then measure whether the inertness extends beyond containment, to ship it documented as a
   known limitation, or to leave it. He chose to leave it: no code change, no help-file change.
   The measurement stands here as a recorded fact about the shipped statistic rather than as a
   pending repair. Stating the limit of the evidence alongside the ruling, because the ruling
   does not extend it: **the ladder is all-containment geometry**, so nothing here establishes
   how `bt` behaves elsewhere, and the 0-of-40 result should not be read as a general claim
   about the bootstrapped adversarial test.
2. **Numerically degenerate pairwise-subgroup inference is DISCLOSED rather than withheld, on
   both builds.** On `ext_cov10_rcs_pweight_ddd_c1` the four `att`-family pairwise rows all carry
   one point estimate - exactly 1 - and the inference published beside it is **se 2.33e-10, t
   4.29e9, p 2.1e-207, CI [1, 1]**. That is a contrast the fit determines exactly, reported as
   though it were estimated. v1.3 published inference on `att_i` alone; HEAD publishes it on
   `att`, `att_it` and `att_itime` and withholds `att_i`. `att` and `att_itime` carry v1.3's
   `att_i` numbers to the bit; `att_it`'s are its own and equally degenerate (se 5.21e-10,
   t 1.92e9). The package's own degeneracy scan sees part of
   this - it records `PWSGS_DEAD(0v+2i/8)` - but not these rows. **Found while checking a ledger
   reason, not by any test**, which is the point: the suite explained this config as a recovery
   for months. The reason text is corrected; the surface itself is the author's call, and it
   belongs with items 1 and 2 as the same kind of problem.
3. **One external defect that is not ours**: the `jwdid` help documents weights its code refuses.
4. **The continuous `ame` line test is about twice nominal size AT the cell size the help
   recommends** - disclosed rather than repaired, and listed here because a reader looking for
   known limitations should not have to find it inside a campaign section. Measured by varying
   observations per cohort-by-period cell and nothing else: rejection 0.402 at six observations
   per cell, **0.110 at fifty**, against a nominal 0.05, with the `att` family from the same fits
   inside its band at both (0.049 and 0.080). So the help's precondition is **necessary but not
   sufficient**: "on the order of fifty" is a large improvement, not a cure. The author ruled the
   cell-size question to user judgement, and the help now carries the measured magnitudes
   (`6ccceec`), so this is a documented limitation rather than an owed repair - but it is a
   limitation, and 10.2 has the mechanism, including the part that is NOT explained (a uniform
   20% understatement of the diagonal predicts at most ~0.22 joint rejection at eight degrees of
   freedom against 0.402 observed, so `e(V_es_ame)` is mis-shaped rather than merely scaled, and
   the measurement covers the diagonal only).

**Worth putting items 1 and 4 side by side, because they were treated differently.** Both are
known limitations of shipped statistics that this assessment measured for the first time. Item 4
is disclosed in the help; item 1 is not disclosed anywhere, by ruling. That is the author's call
and it is recorded as his, not presented as a conclusion of the assessment.

**Closed by ruling rather than by repair, which is a different thing and worth the distinction.**
The efficiency claim about where to set the cohort identifier was REMOVED from the help rather
than measured, at the author's instruction ("I am skeptical of it"). The standard-error
undercoverage candidate therefore has no pending simulation: it stays in Part 8 as a measured
observation, not as an owed task. The honestdid degenerate-covariance guard was ruled to stay as
it is, on the grounds that our handoff need not protect against cases honestdid does not protect
against itself. The silent window truncation was ruled to stay silent.

**Also closed**: the ticket once listed here as an impossible p-value above 1. It was a misreading
of a `%18.12g` field with the exponent dropped by eye; the value is 3.34e-07, it reproduces at the
original seed, no code path can emit a value above 1, and the old build is bit-identical there.

### 7. Fixed during the Stage B triage, 2026-09-13 (each verified in BOTH directions)

| # | what was broken | how it was proven fixed |
|---|---|---|
| 7.1 | The standard `r(459)` trailer printed TWICE on a validation failure. Two `cap noi` call sites used `error _rc`, which reprints a message the callee has already shown, against the file's own rule at `wooldid.ado:61` ("exit, not error - error reprints"). | Same reproduction on the pre-fix and post-fix builds: **2 trailers before, 1 after, return code 459 in both**. The third `error _rc` site was deliberately LEFT ALONE: it follows a `qui cap`, where the message was suppressed, so it is the only disclosure a ppmlhdfe failure gets. |
| 7.2 | The low-`reps()` caution printed BEFORE the `boottstat` nesting check and the two-cluster check, both of which abort - so a user whose run never executed was advised about its p-value resolution. | Moved below every aborting validation. Verified both ways: the aborting run now shows the error and no caution; a healthy run still prints the caution exactly once. |
| 7.3 | `bootsave(name, replace)` - the obvious analogy to Stata's own `saving(f, replace)`, which this option deliberately does not follow - was swallowed whole and wrote a file literally called `name, replace.dta` at **rc 0**. The name check rejected only a directory part. | Now `rc 198` with a message naming `replace` as a separate option, while `bootsave(name) replace` still returns rc 0. No comma-named file is written. |
| 7.4 | The `degeneracy` component detected a curated class moving and emitted only a WARNING, then exited 0 - so genuine drift would print and pass. Its own contract is that a curated expectation still describes the run. | Promoted to a failure, along with a dangling ledger citation. Proven by deliberately corrupting a row: `DEGEN_CHECK_FAIL ... 1 class-moved/dangling`, and clean again when restored. **This hole was in a component written during Stage A of this same assessment**, and was found by an agent auditing it. |
| 7.5 | `DEGEN_SUMMARY`'s class bracket printed SCANNER classes only, so curated `DATASET_DESIGN` rows were counted under `HIDDEN` and the line could not be reconciled against the file it summarises. This is what made a spurious "22 hidden versus 9 expected" discrepancy look like drift. | The line now prints the file's own tally beside the scanner's; the curated counts sum to the row count exactly. |

#### 7.6 A ledger reason that was factually wrong at HEAD (fixed 2026-09-14)

`DDG-pw-1` said v1.3 returns missing pairwise-subgroup cells "where HEAD returns valid values".
A cell-by-cell census of the config's 48 `pwsgs` cells says otherwise: **15 missing to valued, 5
valued to missing, 5 missing on both, 23 valued on both** (8 of those changed value). The five
going the other way are `att_i`'s entire inference row. The move is a relabelling of which
degenerate contrast gets disclosed, not a recovery. `DDG-pw-2` and `-pw-3` had been expanded to
read "as DDG-pw-1 ... missing-cells-recovered class", propagating the same wrong framing to two
more configs; both now say only that the disclosure surface moves, and say explicitly that their
direction was not characterised. Re-judging with the corrected reasons gives verdicts identical to
the run's own and `stale_ledger=0`, so all three rows still fire exactly as before - the text was
wrong, the scoping was not.

### 8. Ledger and expectation rows written by the triage

45 ledger rows across 14 blocks, every one scoped to its config and pattern rather than
blanketing a matrix - the project's own rule, which Stage A had to enforce by retiring six
blanket rows that were swallowing 361 unexamined values. That none of the 45 went stale is the
proof each actually fired. Nine degeneracy expectation rows: six `DATASET_DESIGN` on the biennial
fixture, one `RULED`, and two `HIDDEN`+`UNRULED` for the bias-corrected bootstrap case.


### 9. Fixed by the 2026-09-13/14 rulings (the author ruled; these are the repairs)

| # | what was wrong | how it was proven |
|---|---|---|
| 9.1 | **Absorbed degrees of freedom were over-counted.** Duplicate dose terms are bit-exact copies of columns already in the equation - their job is to steer the collinearity drop, not to be parameters - yet they consumed a degree of freedom wherever a finer set was also absorbed. | Three ways: a rank computation is unchanged; `regress` with explicit dummies agrees with every spelling except the shipped one; and across the suite the residual count moves while the residual sum of squares differs by 0 or 3e-15. After the fix the absorbing standard error equals the non-absorbing one for the same model to sixteen digits, where it was 0.2% off. Battery 46 of 48 identical; no point estimate moved. |
| 9.2 | **The weighted Bayesian bootstrap sampled the wrong estimand.** The Dirichlet parameter came from the cluster's row count and its mass was divided by the cluster's weight total; equal only when the base weight is 1, which every call site happened to be. | Replications sat 9.36 bootstrap deviations off centre with every draw on one side, and the printed interval excluded the user's own estimate. Now 0.10. Unweighted is an exact no-op on draws, coefficients and covariance. The convention the code cited was measured, not assumed, and the comment naming it was false. A new pin FAILS on the pre-fix build and passes after. |
| 9.3 | **`fastnose` withheld the event-study coefficient vector** the help promised it returned. | Code changed to match the help, per ruling; the covariance stays absent, which is the other half of that claim. One battery config moved, additively: the added cells are byte-identical to cells already in the same result. |
| 9.4 | **A caution prescribed a remedy that cannot work.** It advised raising the replication count whenever every draw fell on one side - correct for sparsity, useless when the distribution is displaced. | Qualified so it no longer prescribes unconditionally. |
| 9.5 | **`bootalpha(size)` is now disclosed as unweighted-only** where it claims parity with the reference package, because that package refuses user weights outright. | Measured against the package itself. |
| 9.6 | **Three help claims were false**: `impvalues` listed as returning a vector it does not, the `did_imputation` equivalence conditions missing the one option that makes them true, and the honestdid example missing the option that stops it spawning four child processes. | Each verified against the code or the external package before editing. The suite's own three honestdid call sites were fixed hours earlier; the line a user copies was not. |
| 9.8 | **The two battery baselines the dof fix moved were re-blessed**, `cont_lat` and `perm_cfxnose`. | Not blessed on assertion. Pristine-twice first: the battery was run a second time into a fresh directory and the two runs compared to each other at rtol 0 - `TOTAL_DIFF_GROUPS=0`, all 48 configs bit-identical including both movers, so the movement is deterministic. Then a cell census: 30 and 421 changed cells, of which 24 and 410 sit at or below the 1e-10 indifference threshold. `rss` is unchanged; `rmse` moves 6.6e-4 and the two adjusted R-squareds 2-5e-9, which is the designed effect. The only other above-threshold cells are one `cont_lat` coefficient of value 8.1e-07 whose absolute change is 7.11e-15, and eight `permdist_b` replication cells whose absolute change is 1.08e-11. **All 56 permutation p-value cells are bit-identical**, as are every joint and adversarial statistic; the worst movement across all reported `b`/`V`/`b_es_*`/`V_es_*` cells is 3.27e-11. After blessing, baseline compares 48 IDENTICAL at rtol 0. See Part 7.6. |
| 9.7 | **The coefficient names in `e(b)` were undocumented** - two of them appeared zero times - in a help that tells users eleven times to address coefficients by them. | Now documented where the stored results are. |

### 9b. Fixed by the 2026-09-14 rulings

| # | what was wrong | how it was proven |
|---|---|---|
| 9b.1 | **One degenerate subgroup deleted EVERY subgroup's pre-treatment row**, on both builds. When only one subgroup realises in the pre-period subpopulation, `margins` drops the grouping stripe and returns the target-flag tempvar name; no name carries a subgroup token, so the name match that decides the row found nothing for anyone. The guard written for this shape keyed on the POST-period subgroup count - which is 2 in exactly the case that breaks - so it never fired. | The guard is re-keyed on which subgroups the PRE period actually realises, captured on the same subpopulation as the call that produces the collapsed stripe. On `sg_degpre`, subgroup 1 regains a row that is live rather than residue: estimate 2.08e-4, se 0.067, p 0.9975. Subgroup 2, which genuinely has no pre cells, is still withheld. The three-subgroup case where only the third is degenerate is bit-unchanged - the stripe survives there and the name match does the work. Cert 43/43, battery `TOTAL_DIFF_GROUPS=0`. |
| 9b.2 | **A second defect the same change repairs.** The old guard's first term is vestigial - always true after two `subinstr` strips - so it reduced to "one subgroup in the list implies the row is the LAST one", which would have handed subgroup 1 the wrong row had a one-subgroup list ever met a multi-level pre stripe. | Repaired by the same re-keying. The vestigial term is kept as the file's idiom rather than widening the diff. |
| 9b.3 | **The drop was invisible in both directions** - kept or deleted, nothing posted, no caution, no marker, so neither the user nor the degeneracy scanner nor the ledger could see it. | Three self-gating hidden macros on the existing `e(wdid_sgdegen)` shape record which subgroups lost the row, why, and WHICH row; a per-subgroup caution prints beside the table and survives replay. The third macro exists because the first caution written was wrong: on `cont_sg_gap` the withheld row is the `ame` while `pretreatment:att` sits in the same table with a live estimate, so "the pre-treatment average is not estimable" would have been false. **A disclosure that misstates what it discloses is worse than silence**, and it was caught before it shipped. |
| 9b.4 | **The 105 rewritten ledger reasons reproduced the very defect they repaired.** The honest mechanism is ~1800 characters and near-identical across all 105, so the file went 159KB to 300KB - about 190KB of duplicated prose, which is the shape that let one sentence go stale across 138 rows in the first place. | Stated once as a named `SI20-L-MECHANISM` note and cited by the 105, following the file's own convention where `PS62-esblock` covers a class with one row. Measured first: the 116 `SI20-L` rows carried only **16 distinct texts**, 105 of them beginning with the same 1837-character block. Verified by the comparer's own loader - 437 rows, same id set, **zero rows whose config, kind, pattern or rtol changed**. 300KB to 124KB. |
| 9b.5 | **The help's adversarial-size warning made four claims and none had been measured** on the statistic that ships. | All four now carry numbers, and two were modestly wrong: the "threefold" degradation at two clusters per cell is 4.0x, and the single-cluster case is 16x rather than the "order of magnitude" claimed. See 10.1. |
| 9b.6 | **Was the F1 repair protected against regression?** The agent that built 9b.1 recommended a dedicated cert leg and left it for a ruling. | **RULED 2026-09-14: no cert leg; the battery pin is enough.** `sg_degpre`'s blessed baseline pins the 15-name shape, so an F1 regression turns the battery red today. Recording the contingency the ruling rests on: that protection exists only while `sg_degpre` stays in the battery deck, so if it is ever removed this ruling lapses with it and the cert leg becomes the open question again. |


### 10. The Stage C campaigns: what was measured, 2026-09-14

Three campaigns were preregistered at the start of this assessment - C1 the BS-12 adversarial
size ladder (`mc_adv_bs12.do`), C2 simultaneous sup-t band coverage (`mc_bands.do`), C3 the
randomization-inference extensions (`mc_ri_ext.do`) - and all three shipped as scaffolds that
had never been executed. They are reported here because the chunk asks whether the new features
are correctly sized, and until they ran the honest answer was that nobody knew.

C1 is below at 10.1. **C2 and C3 are in 10.4 and 10.5.** Three further measurements were not
preregistered and were added because the work turned them up: the continuous-treatment
cell-size dose-response (10.2), which exists because C1's continuous rung violates a
precondition the help states; the standard-error coverage MC (10.3), which C.6 asked for; and
the AB12/bt power run, which the author ordered before ruling on the inert statistic and which
is reported with the defect it measures in section 6 rather than here.

#### 10.1 The BS-12 adversarial size ladder - the shipped statistic finally has a number

**Why it was owed.** Every adversarial size figure on record - both Monte Carlo reports and the
help - was measured on the free-intercept statistic that BS-12 replaced on 2026-09-08. The
ruling's own prerequisite simulation was never run. So the package's documented size claims
described a statistic that no longer shipped.

200 replications per cell, B = 199, slope 0, four geometries, 8 shards, about seven hours.
Integrity was clean throughout: `rc = 0` on all 1200 rows, `advabsent` 0, `nfail` 0, `ncomp` =
199 on every resampling row, every p inside [0,1], and `advw_big` 0 on the permutation arm.

| cell | what it is | size at nominal .05 | |
|---|---|---|---|
| `AB1` / an | 6 clusters per cell - the compliant geometry | 9/200 = **0.045** | **PASS** |
| `AB1C` / an | continuous treatment, att family | 10/200 = **0.050** | **PASS** |
| `PM` / pm | permutation p, exact by construction | 7/200 = **0.035** | **PASS** |
| `AB13` / an | 2 clusters per cell | 36/200 = 0.180 | FAIL |
| `AB12` / an | containment, 1 cluster per cell | 145/200 = 0.725 | not scored |
| `AB12` / bt | bootstrapped, containment cell | 0/200 = 0.000 | FAIL |

**The headline is the three passes.** The statistic the package actually ships is correctly
sized where its help says it should be, including under a continuous treatment on the att
family, and the permutation arm is exact. That is the first evidence for BS-12 rather than for
its predecessor.

**The two failures are both regimes the help already describes, and they put numbers on claims
that had none.** The help said size rises "roughly threefold when cells were spanned by only two
clusters": measured 4.0x. It said "an order of magnitude when each cell lay entirely within a
single cluster": measured 16x, worse than its own estimate. Both corrections are now in the help.

**`AB12/bt` is new and is open.** See section 6 item 3.

#### 10.2 The continuous cell-size precondition - a caveat confirmed, and a correction to my own reading

The help states one condition specific to the continuous case: cells "on the order of fifty"
observations, because each cell's dose slope is otherwise too imprecisely estimated for the
test's covariance to be trustworthy. It had never been measured, and - this is the part worth
recording - **the rung preregistered as the continuous half of the ask violates it**. `AB1C` is
30 cohorts x 6 units x 12 periods, so a treated cohort-by-period cell holds six observations
against the help's fifty.

That was found the hard way. A probe measured the `ame` line test rejecting at 40% against a
nominal 5% and it looked like an undocumented defect in a shipped feature. It is not: the help
predicts it, and the probe was measuring a design the help warns against. **A measured
over-rejection is not a finding until the help has been read.**

Measured properly, by varying observations per cell and nothing else:

| observations per cell | `ame` p at nominal .05 | uniformity | `att` p from the same fits |
|---|---|---|---|
| 6 | **0.402** | KS 0.520, rejected | 0.049, inside band |
| 50 | **0.110** | KS 0.109, not rejected | 0.080, inside band |

So the caveat is real, the direction is right, and fifty observations per cell is a large
improvement rather than exactness - 0.110 is still about twice nominal. The help now says that.

The mechanism is only partly diagnosed, and the report says so rather than overclaiming. Over
200 draws the reported `ame` standard errors run at 0.800 of the sampling standard deviation,
every one of eight horizons below 0.85, while the `att` family sits at 1.036 with every horizon
in [0.96, 1.14]; `e(V_es_ame)` is posted on every draw, so the documented cell-containment screen
never fires here. But a uniform 20 percent understatement of the diagonal predicts at most about
0.22 joint rejection at eight degrees of freedom, and 0.402 was observed - so `e(V_es_ame)` is
**mis-shaped, not merely scaled down**. The measurement records the diagonal only and therefore
establishes that the variance is understated, not by how much in which direction.

**What this cost the assessment, stated plainly.** A preregistered sizing rung was built on a
design its own help declares untrustworthy, and neither the preregistration nor the driver
header recorded the precondition. That is a design error in the assessment, not in the package.
Both are amended, and the compliant measurement now ships as `mc_ame_cellsize.do` so the help's
numbers are reproducible from the repository.

#### 10.3 The standard-error coverage MC - C.6's candidate, confirmed at 500 draws

Part 8 C.6 flagged wooldid's analytic standard errors running 20-40% below the sampling standard
deviation on two arms, with csdid and did_imputation covering correctly on IDENTICAL point
estimates in the same draws. It rested on 30 draws at one cluster count and C.6 itself said that
was not ruling-grade. The author's help-side ruling of 2026-09-13 closed the documentation half by
removing the sentence that steered users toward the geometry; it did not measure anything. On
2026-09-14 he ruled to run the measurement at G = 40. Medians across horizons, 500 draws:

| arm | hom SE/SD | hom cov | stag SE/SD | stag cov |
|---|---|---|---|---|
| `esfixedbaseperiod`, **i = unit id** | 0.949 | 0.934 | **0.742** | **0.838** |
| `esfixedbaseperiod`, i = cohort id | 1.019 | 0.953 | 0.972 | 0.934 |
| pooled `esprelength(0)`, **i = unit id** | **0.799** | **0.876** | **0.661** | **0.794** |
| pooled `esprelength(0)`, i = cohort id | 1.010 | 0.950 | 1.006 | 0.946 |
| csdid `long2` *(control)* | 0.961 | 0.939 | 0.927 | 0.926 |
| did_imputation *(control)* | 0.963 | 0.939 | 0.943 | 0.936 |

**On the common-onset DGP C.6's own pre-written flag rule fires cleanly**: the pooled arm at
`i` = unit reports SEs at 0.799 of the sampling SD and covers at 0.876 against a nominal 0.95
while both controls cover at 0.939. **S2 is confirmed at 500 draws.**

**On the staggered DGP the rule abstains, on a technicality that is stated rather than hidden**:
csdid's coverage is 0.926, four thousandths below the rule's 0.93 gate, so "both controls inside
band" fails. The magnitudes are not close - wooldid's two `i` = unit arms cover at 0.838 and 0.794.
C.6 predicted 0.838 on that arm from thirty draws and five hundred draws return 0.838. Read the
gate as abstaining, not exonerating.

**The decisive detail is per horizon.** On the staggered DGP the empirical SD of wooldid's
estimates and csdid's are identical to four decimals, horizon by horizon - 0.1563, 0.1645, 0.1489,
0.2339 - because they ARE the same estimates. What differs is the reported standard error: 0.1219
against 0.1530 at the first pre horizon, 0.0845 against 0.1389 at the fourth. Same estimates, same
sampling variability, different reported uncertainty, and no DGP explanation is available for it.

**`i` = the cohort id is clean everywhere** (0.95-1.02, 0.934-0.953), across both DGPs and both
estimators. The effect is specific to `i` = the unit id, which is what the removed help sentence
used to recommend - so the ruling that removed it is validated by this run.

Scope: one cluster count, one panel length, `cluster(i)` only. C.6 asked for a G sweep and this run
was scoped to G = 40, so it establishes the effect is real and large there and says nothing about
how it behaves in G. **That is a ruling, not an omission.** A G = 200 follow-up (one shard,
about four hours) was offered on 2026-09-14 on the explicit ground that the original concern was
whether the shortfall is a small-G artifact and that G = 40 is itself small; the author ruled
SKIP. The G = 40 result therefore stands as the scoped answer and the sentence above is the
honest boundary of it. Two staggered horizons are far worse than their neighbours (0.567, 0.591) and
that is unexplained. Report and driver: `review_process/se_coverage_2026-09-14/`.

#### 10.4 C2 - the simultaneous sup-t band, measured for the first time (`558f18d`)

**Why it needed measuring.** The band shipped on 2026-09-10, *after* every Monte Carlo in the
project had already run. `test_esband.do` pins its structure - that the matrix exists and has the
right shape - and nothing pinned its **coverage**. No MC had ever used a cluster strictly below
`i`, which is both the geometry the help recommends and the only one on which the analytic
`e(V_es_*)` the band is simulated from is the right object. So the headline claim of the feature
was unverified.

**Design.** 150 reps per arm, 6 shards, B = 199, tau = 0.5, three engines: analytic-simulated
(`an`), bootstrap (`bs`), bootstrap-studentized (`bt`). 450 rows, no `level()` passed, so
wooldid's default 95 is in force.

**The result: every preregistered check passes, 0 failures.**

**P1, the headline.** Simultaneous coverage against nominal 0.95, preregistered count band
[137,148]:

| arm | engine | `band_all` | rate | verdict |
|---|---|---|---|---|
| an | analytic-simulated | 139/150 | 0.9267 [0.8735,0.9586] | PASS |
| bs | bootstrap | 139/150 | 0.9267 [0.8735,0.9586] | PASS |
| bt | bootstrap-studentized | 146/150 | 0.9733 [0.9334,0.9896] | PASS |

The margin is worth stating rather than hiding behind the verdict: `an` and `bs` sit **two counts
above the band's lower edge**. That is a pass and not a close one in any statistical sense - the
Wilson interval comfortably contains 0.95 - but they are not central either, and if this band is
ever re-measured at higher n the honest prior is that analytic and bootstrap sit a little under
nominal while the studentized arm runs slightly conservative.

**P2, the control.** All 8 horizons of all 3 arms lie inside their pointwise binomial regions.

**P3, the directional check, and the sign IS the claim.** A simultaneous band that is not wider
than the pointwise intervals is not doing its job. With 8 horizons, intervals claiming only
pointwise coverage should cover *all* of them well below 95%. They do, and the gap is large and
consistent:

| arm | `pw_all` | `band_all` | gap |
|---|---|---|---|
| an | 0.7000 | 0.9267 | +0.2267 |
| bs | 0.7067 | 0.9267 | +0.2200 |
| bt | 0.7133 | 0.9733 | +0.2600 |

**P0 integrity, clean on all 450 rows.** `rc = 0`; a band present on every rep
(`band_absent = 0`, so no rep was silently dropped by the kernel refusing a degenerate horizon);
`nh = 8`; `pw_n = nh`, so the control is measured on the same horizon set as the band;
`ncomp = 199` and `nfail = 0` on both bootstrap arms; and the summary columns arithmetically
consistent with the per-horizon ones on every row. That last family exists because **a driver
that mis-writes its summary columns turns a coverage number into fiction without failing
anything**, and it is the check that would have caught it.

**One secondary fact that carries the mechanism.** Mean critical value is 2.644 (`an`), 2.658
(`bs`) and 5.049 (`bt`), against a pointwise 1.96. The critical value both exceeds the pointwise
one and **moves with the engine** - if it did not, the simultaneity would be nominal only, a
relabelled pointwise interval. It is real.

**Scope.** One geometry (G = 30, U = 6, T = 12, cluster strictly below `i`), one window
(`esprelength(4) espostlength(4) esfixedbaseperiod`), one level (95), one tau. It establishes
that the band covers at its nominal level where the help says to use it; it says nothing about
other windows or levels. Data and report: `test/suite/mc/out3/bands/`.


---
