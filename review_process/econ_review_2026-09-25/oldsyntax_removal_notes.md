# Removal of the `oldsyntax` option

Worktree: /Users/hegland/Documents/Research/wooldid_dev/.claude/worktrees/agent-a35e83036e5ce789e
I ran `git reset --hard 925d0e8` first, and `git log --oneline -1` then printed
`925d0e8 Record the docs audit, its post-apply check and the verification of the round-2 tree`.
Nothing is committed. The diff is in `removal.diff` in this directory: 6 files, 28 insertions, 74 deletions.

## Sites changed

### wooldidinstall/wooldid.ado
Line numbers are from 925d0e8.
- **753 and 5388 (the two syntax lines):** removed the `oldsyntax` token with one substitution (` update oldsyntax PRETRENDPOWer(numlist) ` -> ` update PRETRENDPOWer(numlist) `) applied to both lines, so they are still byte-identical.
- **5405-5411 (validate):** deleted the normalization block and the comment above it (`* This controls the names we use for ES plots; ...`), plus the blank line that followed it.
- **3565 and 4049:** removed ` old(`oldsyntax')` from both `wooldid_sglbl` calls.
- **14215 (`wooldid_sglbl` syntax):** removed the required option `old(integer)`.
- **14284-14286 (comment):** removed the three oldsyntax sentences. The membership-list rule is kept in one line: `* lbl() is a MEMBERSHIP list consumed by a strpos() test.` The existing "Never index lbl() positionally ..." lines follow it unchanged.
- **14293:** deleted `if `r' == 0 & `old' == 1 local lbls `lbls' contemp`.
- **4327, 4339, 4352, 4364:** deleted the four `if `relyr' == 0 & `oldsyntax' == 1 local name "contemp"` lines. The `post`relyr'` / `pre`arelyr'` lines are unchanged.
- **4646 and 4732:** deleted the `else if `f8rv' == 0 & `oldsyntax' == 1 ...` branches. Each chain is now `if <0 pre... / else post...`. When oldsyntax was 0 the deleted branch was always false, so every case still takes the same path.
- **5763 and 5776:** removed ` os(`oldsyntax')` from the `wooldid_permctl` and `wooldid_bootctl` calls.
- **6623 and 7129:** removed `os(string)` from the syntax of both programs. I checked: `os'` was read in exactly one place in each program, the refusal.
- **6660-6663 (permctl):** deleted the refusal "Error: permute() is not compatible with oldsyntax.".
- **7401-7405 (bootctl):** deleted the two comment lines and the refusal "Error: bootstrap() is not available under oldsyntax.".
- **6322:** removed `oldsyntax` from the `foreach wdv_ret in ...` c_local list.
- **Mata `wooldid_esrelx` (18686-18696):**
  - Deleted `if (s == "contemp") return(0)`.
  - Reflowed the comment above it to drop the "contemp under oldsyntax" clause. The "it_post0 tail returns missing, which makes the scan self-disambiguating" sentence is kept.
- **Mata `wooldid_vspl` (19012-19031):**
  - Declaration `string scalar want, want2` -> `string scalar want`.
  - Deleted the two comment lines, the `want2 = ...` line, and the `else if (want2 != "" & ...)` branch.
  - The inner loop is now just `if (vs[i, 2] == want) idx[j] = i`. Without oldsyntax no `px+"contemp"` row existed, so the else-if never fired.

### wooldidinstall/wooldid.sthlp
- **110:** `{opt safety(on|off)}   {opt oldsyntax}  {p_end}` -> `{opt safety(on|off)} {p_end}`.
- **358 (bootstrap refusal list):** removed `{it:oldsyntax},`. The line now ends `{it:customvce()},` and line 359 continues with `{it:suppresssgprimaryeffects} and multi-way ...`.
- **436 (permute refusal list):** removed `{it:oldsyntax}, `.
- **924-927:** removed the whole `{phang}{opt oldsyntax}` entry and one blank line. The safety() `{p_end}`, one blank line and `{marker runtime}` now follow each other.
- Braces balance on every touched line (5/5, 2/2, 7/7).

### wooldidinstall/wooldid_agent_guide.md (chapter 14 not touched)
- **164:** removed `` `oldsyntax` `` from the Technical Features list.
- **277:** removed `` `oldsyntax`, `` from the bootstrap refusal list.
- **283 (NOT in the brief):** removed `` `oldsyntax`, `` from the incompatibility list in the `boot:strap(mode)` entry.
- **375:** removed `` `oldsyntax`, `` from the permute refusal list.
- **383 (NOT in the brief):** removed `*oldsyntax*, ` from the "Not compatible with" list in the permute entry.
- **663:** the sentence now ends "... when `reghdfe` leaves the session in a bad state." The clause "`oldsyntax` restores the old `contemp` coefficient suffix and is refused by both replication engines" is gone.
- **769 (engine-compatibility table row):** `impvalues`, `unconditionalse`, `oldsyntax` -> `impvalues`, `unconditionalse`.
- **799:** removed the whole `oldsyntax` entry and its trailing blank line.

### test/test_bootstrap_gates.do
- **112-113:** deleted gate 5f (the `oldsyntax` call and `chk 5f 198 _rc`). 5g and 5h keep their numbers.
- **287-288:** reworded the comment to "(a normalization bug once rode under matching 198s)", with no mention of oldsyntax.
- **366-368:** deleted the m4 leg (the `log using boot_gates_m4.smcl`, the oldsyntax call, and `logh m4 not available under oldsyntax`). m5 and later keep their numbers.
- The file keeps no gate count. Its only tally is `BOOTG_FAILS` and the `BOOTGATES_ALL_PASS` sentinel, which matches run_full.sh's `'^[A-Z0-9_]+_ALL_PASS'`, and neither changed.
- No file ever erased `boot_gates_m4.smcl`, so there was no cleanup line to remove.
- Neither deleted leg drew random numbers: both were refused at parse, before any draw. So the later legs' RNG streams do not move.

### test/suite
- **configs_ext.csv:321:** deleted the `ext_cov10_oldsyntax` row. File otherwise byte-identical (LF endings kept).
- **README.md:182:** removed `` `oldsyntax`, `` from the "display/compat switches" list.
- After the edits, a grep of test/ (excluding test/suite/out and *.log) finds no `oldsyntax` and no `ext_cov10_oldsyntax`. ledger.csv, degeneracy_expected.csv, data_classes.csv, preps.json and skipped_cases.txt never referenced either.

## check_validate_sync
The script is at `test/check_validate_sync.py`, not `test/suite/`. I ran it from the worktree root:

```
OK: worker and validate syntax lines are byte-identical
OK: c_local closure holds (43 assigned, 97 returned)
OK: every syntax option is covered by the c_local blanket
OK: the permutation draw-matrix cap is 10000 at all three sites
OK: the bootstrap draw-matrix cap is 10000 in the e(bootdist) test, the e(bootdist_pwse) test, the shell caution's test and its text, and the caution prints below part(boot)
exit=0
```

The same script run against the unedited ado printed `44 assigned, 98 returned`. The drop of exactly one on each side is `oldsyntax`.

## Greps
| Command | Count |
|---|---|
| `grep -n -i oldsyntax wooldidinstall/wooldid.ado` | **0** |
| `grep -n contemp wooldidinstall/wooldid.ado` | **0** (also 0 case-insensitive) |
| `grep -n -i 'oldsyntax\|contemp' wooldidinstall/wooldid.sthlp` | 0 |
| `grep -n -i 'oldsyntax\|contemp' wooldidinstall/wooldid_agent_guide.md` | **2**, both in GENERATED chapter 14 |

The two guide hits are at lines 1181 and 1239 (1183/1241 at 925d0e8): the permctl and bootctl rows of the error table. I left them as instructed. Regenerating chapter 14 from the edited ado should drop both rows.

## Smoke check
Scripts: `runs.do` (with wrappers `smoke_new.do` and `smoke_old.do`) and `cmp.do`.

**Setup:**
- Stata ran as `stata-mp -e do ...` from this directory, with at most 2 instances of mine at once. `pgrep` showed 0 before launch.
- Data: the 12-unit panel copied from test/test_bootstrap_gates.do (`y unit t ttre`, `set seed 11`), plus `sgv = 1 + mod(unit,2)`. Leg 2d uses test/suite/data/continuous_control_dense.dta.
- Each session puts its build first with `adopath ++`, and `which wooldid` confirmed the build:
  - new: `<worktree>/wooldidinstall/wooldid.ado`
  - old: `old/wooldid.ado`, from `git show 925d0e8:...`
- The old and new builds ran in separate sessions.

**New build (`smoke_new.log`), ending in `SMOKE_new_ALL_PASS`:**
- (1) `wooldid y unit t ttre, att`: ran.
- (2) `wooldid y unit t ttre, att esprelength(2) espostlength(2)`:
  - ASSERT: `strpos(colfullnames e(b), "post0") > 0` and `strpos(..., "contemp") == 0`, which passed ("r2 ASSERT OK").
  - The names are `main_att pretreatment_att ES_att_pre2 ES_att_pre1 ES_att_post0 ES_att_post1`, with equation-qualified forms `ES_att:post0` and so on.
- Extra legs, to exercise the edited code paths:
  - (2b) `att att_it esprelength(2) espostlength(2) simultaneousbands cluster(unit)`: `wooldid_esrelx`, `wooldid_vspl`, e(esband_*).
  - (2c) `att esprelength(2) espostlength(2) subgroups(sgv)`: `wooldid_sglbl` binary call, the f8 subgroup block, e(escomp_s*).
  - (2d) `y_cont_single cohort t ttre_single, att att_it contreat(dose_common) contreatcontrols(basic) subgroups(sg_three) esprelength(1) espostlength(1) robust safety(off)`: ES_ame rows and the `wooldid_sglbl` continuous-treatment call.
- (3) `att esprelength(2) espostlength(2) cluster(unit) bootstrap(standard) reps(20) rseed(1) simultaneousbands`: rc 0.
- (4) `att esprelength(2) espostlength(2) cluster(unit) permute(shuffle) permreps(20) rseed(1)`: rc 0.
- After every run, ASSERT: no row or column name of any e() matrix contains "contemp". Passed for all 7 runs.
- (5) `wooldid y unit t ttre, att oldsyntax`:
  - ASSERT: `_rc == 198` and the log contains "option oldsyntax not allowed".
  - Passed: rc 198, message `option oldsyntax not allowed`.

**Old build (`smoke_old.log`):** the same runs 1, 2, 2b, 2c, 2d, 3 and 4 all ran with rc 0, ending in `SMOKE_old_ALL_PASS`.

**Comparison of old and new (`cmp.log`), ending in `CMP_ALL_PASS`:**
- `cmp.do` loads each pair of `.ster` files and checks that both builds post the same sets of e() matrix, scalar and macro names.
- It then compares every e() matrix (mreldif == 0 plus identical colfullnames and rowfullnames), every scalar (bit-exact in `%21x`) and every macro (exact string compare in Mata).

| Run | mreldif(b) | mreldif(V) | b colnames identical | V rownames identical | Matrices / scalars / macros | Differ |
|---|---|---|---|---|---|---|
| r1 | 0 | 0 | 1 | 1 | 3 / 14 / 14 | 0 |
| r2 | 0 | 0 | 1 | 1 | 6 / 22 / 15 | 0 |
| r2b | 0 | 0 | 1 | 1 | 10 / 28 / 20 | 0 |
| r2c | 0 | 0 | 1 | 1 | 15 / 26 / 17 | 0 |
| r2d | 0 | 0 | 1 | 1 | 43 / 45 / 20 | 0 |
| r3 | 0 | 0 | 1 | 1 | 29 / 36 / 25 | 0 |
| r4 | 0 | e(V) absent in both* | 1 | n/a | 23 / 36 / 22 | 0 |

\* permute() fast mode posts no e(V) in either build. The r4 macros compared include `perm_rngstate`.

**Printed output:**
- I compared each run's text log between the builds after stripping only the log-header lines (log path, opened/closed timestamps).
- All seven are IDENTICAL: r1 31 lines, r2 49, r2b 66, r2c 134, r2d 403, r3 95, r4 108.
- So every caution, note and table is unchanged.

**Not run:** the cert suite and test_bootstrap_gates.do itself, as the brief said. The test-file edits are pure deletions plus a reworded comment.

## Things I was unsure about / left for the coordinator
1. **The brief missed two guide sites.** Lines 283 and 383 of the agent guide also named oldsyntax. I removed the word there too, so the guide is clean outside chapter 14.
2. **`test/suite/compare_oldnew.py:103`** still has `ES_ROW = re.compile(r"^(ES_.+?)_(pre\d+|post\d+|contemp)$")`. It does not say `oldsyntax`, so it was outside the grep scope. It is now dead but harmless (it could still matter when comparing against old dumps that ran ext_cov10_oldsyntax). I left it.
3. **`test/suite/coverage_check.py:54`** has the comment "Measured on the current syntax line: 85 options". Its own parser counts 83 options at 925d0e8 and 82 now. The number was already stale before this change, so I left it.
4. **The README phrase "display/compat switches"** no longer has a compat switch in its list (verbose, regtol, cleanmatrices, emptycellsoverride, suppresssgprimaryeffects). I only removed the word, as instructed.
5. **Suite baselines in test/suite/out/** (not touched) presumably still hold ext_cov10_oldsyntax results. The next dual-build or compare run may report the config as gone, or need a baseline refresh.
6. **The external bundle `wooldid_dev/wooldid_external/`** (untracked, main checkout only) still defines the cov10_oldsyntax case: coverage_extensions.do, manifest cases.csv, requirements.csv and the b003 expected manifest. `mine_external.py` writes only the `.GEN` files, never the curated configs_ext.csv, so a regeneration would bring the row back only if someone merges it. I did not touch the bundle.
7. **Hand-written comment text:** the new sglbl comment line, the esrelx comment reflow, and the test_bootstrap_gates.do comment ("a normalization bug once rode under matching 198s") are the only new prose. Everything else is deletion.
8. **Mata robustness:** `wooldid_vspl`, `wooldid_esbanda` and `wooldid_bootband` are called under `cap mata:`, so a broken function would fail silently. r2b and r3 exercise them: e(V) off-diagonals, e(esband_att), e(esband_att_it) and the bootstrap band. They match the old build bit for bit, which rules that out.
   Matching alone would also hold if the splice had failed in both builds, so `vsplchk.do` checks the new build directly:
   - r2: e(V)[ES_att_post0, ES_att_post1] = .18781937
   - r2b: e(V)[ES_att_it_post0, ES_att_it_pre1] = .18722086

   Both are non-zero cross-horizon covariances, so `wooldid_vspl` did splice the family blocks in the new build.
