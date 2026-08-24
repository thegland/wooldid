# Old-wooldid bug attribution report — full assessment run, 2026-08-23

**Claim proven by this report:** every old-vs-new difference observed in the 2026-08-23
full assessment run traces to a documented, verified bugfix in old wooldid or to a
documented intended change. No unexplained divergence remains; no regression in current
wooldid was found.

**2026-08-24 extension (§5):** the new-feature statistical-validation phase. Clean
bills of health for randomization inference (exactly sized, 7/7 preregistered MC
configs), pretrendpower's power claims (externally corroborated by Roth's `pretrends`
package), the F2 event-study b/V returns (reghdfe oracle at machine precision), the
honestdid handoff, and fastnose parity — plus one headline documented issue class
(NF-1/NF-2/NF-3): V_es-based joint tests (adversarial, jointtests, roth pre-trend
Wald, and all continuous-treatment ame ES inference) are severely oversized when the
cluster VCE is blind to within-cell variation, a long-standing property shared with
old wooldid, now mechanism-proven and quantified. Fixes deferred to the bug-fixing
chunk by work-file order.

Builds: **old** = published wooldid v1.3 (`git main`, prov `old|301b706`) vs
**new** = git HEAD (prov `new|e53671f`). Environment: Stata 19 MP, macOS, gtools
present. Sources of record: divergence ledger `test/suite/ledger.csv` (row ids cited
throughout), per-config firing log `test/suite/out/oldnew_compare_verbose.txt`, proof
narratives `wooldid_solved_issues.md`, prior certification
`test/suite/REPORT_2026-08-18.md`, preregistration `test/suite/PREREG_2026-08-23.md`,
run record `test/suite/out/SUITE_REPORT.md`.

---

## 0. Executive summary of the run

**Result: 450/450 configurations compared, 0 unexplained** (after three rc-pin
updates and one deliberate coverage addition, below). Final verdicts: **376
EXPLAINED** (every value/name/rc difference matched a cited ledger row: 3,878
explained difference groups, 11,783 informational new-only additions), **74
RC_EXPLAINED** (42 both-error parity pins + 32 designed return-code changes), **0
stale ledger rows, 0 unscoped rows, 0 provenance failures**.

Companion components, all green: **cert 21/21** (the 20 targeted `test/test_*.do`
suites plus `check_validate_sync.py`), **battery IDENTICAL** (29/29 configs bit-exact
vs the committed baseline at `--rtol 0` — the forward-regression guard), **coverage
OK** across the option surface. `newfeat` was deliberately excluded this run per the
preregistration (author holds new-feature validation for a later chunk).

### Run mechanics and the three rc-pin updates

- **21:20 full run** (`data oldnew cert battery coverage`): 375 EXPLAINED,
  70 RC_EXPLAINED, **3 RC_UNEXPLAINED**; cert/battery/coverage all green.
- The three failures were exactly the stale-pin class the preregistration's P1
  predicted: configs whose rc pins were frozen at `59b9183`, **predating the
  option-domain screens of commit `b42b71f`** (VAL-G19/G21A/B/C/E, landed 2026-08-18),
  which converted five silent-accept paths from rc 0 to informative r(198). Pin
  updates (rc_new `` → `198` in `configs_ext.csv`; dumps unchanged):
  1. `ext_contract_bad_elasticity` — `contreatelasticitytype(notAType)` silently
     accepted by old and by pre-campaign HEAD; now r(198) (`b42b71f`, VAL-G21).
  2. `ext_contract10_poly_zero` — `contreatpoly(0)` ran identically to order 1 on
     old; now r(198) (`b42b71f`, VAL-G21B — confirmed in that commit's message and
     test diff).
  3. `ext_contract10_bad_safety` — `safety(nonsense)` silently accepted; now r(198)
     (`b42b71f`, VAL-G21).
- **21:44 compare re-judge** of the same dumps: 375 EXPLAINED, 73 RC_EXPLAINED,
  **0 unexplained** — PASS.
- **Adversarial pin verification (independent reviewer, refute-the-citation
  mandate): all three attributions CONFIRMED.** Each config was re-run on both
  staged builds: the new-side error text is the `b42b71f` screen verbatim (git
  `log -S` dates each message uniquely to that commit; the text is absent from the
  old build). The verification also characterized what old silently DID with each
  garbage value — see sections 1.16/2 — and swept both config CSVs for any other
  configs exercising these screens: none exist, so no stale pins remain.
- **Coverage closure (the preregistration's P3 gap):** an independent prober
  confirmed NO config reached the `7ba1e16` esf support-first corner (the suite's
  clean run against a ledger with zero `7ba1e16` rows is itself the proof), then
  designed and empirically validated two new `trig`-family configs demonstrating it
  (section 1.17): `trig_esfsup_ref` (safety(off): old silently wrong, four new
  LARGE ledger rows `G10B-esf-s/m/e/n`) and `trig_esfsup_ref322` (default safety:
  old r322 opaque abort → new graceful rc 0, pinned). After the dual-build
  `ONLY=trig` re-run, the **final full re-judge: 450/450, 0 unexplained, 0 stale,
  0 unscoped — PASS.**

### Preregistration scorecard (PREREG_2026-08-23.md)

- **P0** (health claims): all satisfied — battery 29/29 at rtol 0, cert 21/21,
  coverage OK, zero diffs attributable to the tranche-17 commits.
- **P1** (stale rc pins): confirmed — exactly 3, all in the `b42b71f` screen family.
  No stale pins surfaced from the other screen-changing commit families P1 listed
  (270854d, 03cb89e, 970efc8, cb1996c, df70e81) — their reachable configs carry
  current pins.
- **P2** (pin↔ledger STALE fallout): none occurred (`stale_ledger=0`).
- **P3** (possible uncovered surfaces): confirmed for the `7ba1e16` esf
  support-first corner (no config fired it) — and CLOSED in-run by the two new
  demonstrating configs of section 1.17. The NFADV pooled-reference/degenerate-q
  unpost paths remain ledgered only for `ext_*` configs (section 4).
- **P4** (headline prediction): confirmed — the unexplained pile was entirely pin
  staleness; **plausible real regressions found: zero**.

### Why MATCH is 0 this run (it was 133 on 2026-08-18)

Not a regression signal. Commit `aec904e` (CORE-10) makes `e(cmdline)` begin with the
command name, so every rc-0 config now differs from old in exactly that macro. The 133
formerly-MATCH configs all moved to EXPLAINED via the single ledger row
`CORE10-cmdline` (fired on all 375 rc-0 configs), with **zero numeric differences** on
those configs beyond their previously-documented classes. See section 3.1.

---

## 1. Fixed old-wooldid bugs that explain mismatches

Each entry: what old wooldid did wrong and its user impact; the mechanism; the fix
commit(s) with the decisive verification evidence from `wooldid_solved_issues.md`; and
the suite configs demonstrating it in THIS run with observed magnitudes and ledger row
ids. Old-code line numbers refer to v1.3 (`301b706`).

### 1.1 Summary statistics: "control group" rows were pre-treatment-treated rows (B1) and sub-sample IQRs were the full-sample IQR (B2)

- **Symptom / impact:** every `y_control`/`c_control` row of `e(summarystats)`
  duplicated the `*_pretreat` row (the `if` condition at old line 1379 was identical
  to the pre-treatment block's — `ttre <. & t < ttre` — instead of selecting
  never-treated units), and every pre/post/control IQR cell plus the derived
  `e(*_iqr)` scalars carried the full-sample IQR (lines 1359/1372/1385 read the
  in-sample `p75/p25` locals). Users citing these descriptive tables published wrong
  control-group and dispersion statistics.
- **Fix + proof:** `015ae88` (B1) and `0188b8f` (B2), 2026-07-19. Parent probes:
  reported control mean = pretreat mean = 2.1575 while the true control mean was
  2.1257; reported pre-IQR = post-IQR = in-sample IQR = 1.5322 vs true 1.3513 (pre)
  and 1.3982 (post). Both verified by 3-assessor + 3-verifier waves; only the
  affected summarystats/IQR cells changed, estimation bit-identical.
- **This run:** ledger rows `A1-B1a..e`, `A2-B2a..e` fired on `basic` (worst
  summarystats reldiff 2.92e-01; `ypre_iqr` 1.5684 → 1.3347), `pw_ss` (2.92e-01),
  `cont_ss` (3.33e-01, contreat rows included), `adm_base` (1.68e-01 on the realistic
  admin panel), and the bundle configs `ext_bin_ext_support`, `ext_cont_ext_support`,
  `ext_contract10_summary_control`. Also bundled into `cov_kitchen` (row `A7-kit-m`:
  summarystats worst 1.48e-01).

### 1.2 `impvalues` corrupted the returned `e(ll)` (B4)

- **Symptom / impact:** with `impvalues`, `e(ll)` silently returned the last
  Ibragimov–Müller CI lower bound instead of the log-likelihood (local `ll` reused at
  old lines 1568 et al., then ereturned at 3070). Any model-comparison or IC built on
  `e(ll)` was garbage.
- **Fix + proof:** `4de841d` (2026-07-20), IM locals renamed. Parent probe: `e(ll)` =
  −1005.89 without `impvalues`, = 0.584 (a CI bound) with it. Battery showed exactly
  one changed value across 18 configs — `e(ll)` under impvalues, now the true ll.
- **This run:** rows `A3-B4/B4b/B4c`, `B4-ext-1/2` on `impv` (0.612 → −400.91),
  `adm_impv` (0.457 → −38745.48), `trig_imdegen` (1.875 → 4106.64; positive ll is
  correct on that near-noiseless panel), `ext_bin_ext_impvalues`,
  `ext_cov10_cont_impvalues`.

### 1.3 cfx matrices truncated to width 1 on every estimand after the first (B7)

- **Symptom / impact:** with `makecfxplots` + multiple estimands, old overwrote the
  `cfxplotnose` option flag with `marginsdof` contents mid-loop, truncating all cfx
  matrices after the first estimand to the estimate column only — SEs/CIs silently
  lost (and plot titles corrupted by a tempvar name). Assessment proved the truncation
  was cluster-independent and positional: every estimand after the first was hit.
- **Fix + proof:** `3df74db` (2026-07-20), per-iteration copies. Parent probe:
  `e(cfxatt_cfxpost)` had 6 columns while `e(cfxatt_it_cfxpost)` from the same run had
  1. Verified 3/3: estimand-2+ matrices restored full-width, bit-identical to
  single-estimand references.
- **This run:** rows `A5-cfx-m/n` on `trig_cfx`: `cfxatt_it_*`/`cfxame_it_*` matrices
  3x1 on old → 3x6 on new (the 5 inference columns ABSENT on old). Row stripes also
  differ by design — see NFFAST-02, section 3.2.

### 1.4 Extended missing values (`.a`–`.z`) escaped the missingness screens (B9)

- **Symptom / impact:** old tested `== .` (lines 421/443), so `.a`-style missings —
  routine in survey data (MEPS, CPS) — stayed in the screened sample: they inflated
  the `w_per_it`/`w_per_i`/`t_per_i` weight denominators (diluting real observations'
  weights) and `t = .a` on a treated unit polluted the post-treatment classification
  (`t >= ttre` evaluates true).
- **Fix + proof:** `2adf447` (2026-07-20), `>= .` at both sites. Verified 3+3:
  corrupted-data runs bit-match cleaned-reference runs across every injection point,
  estimand, weight, and estimator; clean-data behavior provably bit-identical.
- **This run:** row `A6-ext-m` on `trig_extmiss` (accepted blind-spot config, built to
  fire this class): b worst reldiff 5.07e-02, `wdidmainresults` worst 1.36e-01 —
  old's weights were diluted by phantom rows; new matches the cleaned-data answer.

### 1.5 Estimand weights built on the pre-regression sample, not `e(sample)` (B10)

- **Symptom / impact:** reghdfe singleton drops (and B9-class rows) made the
  regression sample strictly smaller than the screened sample, but the att_it/att_i/
  att_itime cell weights normalized over the larger set — cells containing dropped
  observations were silently underweighted, violating the documented equal-cell-weight
  semantics. Assessment proved old's att_it collapsed exactly onto plain att whenever
  a whole cell dropped.
- **Fix + proof:** `40e1913` (2026-07-20), weight construction moved after estimation
  onto `e(sample)`. Verified 3+3: under singleton/separation drops, full-data runs
  match cleaned-reference runs bit-for-bit; plain att provably invariant; bit-identical
  on no-drop data across 18+ configs.
- **This run:** row `A7-drop-m` on `trig_allest_drop` (blind-spot config with forced
  fe-singleton drops): b worst 8.22e-02, `wdidmainresults` worst 2.50e-01, while the
  plain-att cells agree to ~4e-09 (the documented signature: only the cell-weighted
  estimands move). Rows `A7-kit-s/m` on `cov_kitchen` combine this esamp weight base
  with B1/B2 (b worst 1.41e-07 plus the summarystats class).

### 1.6 Sentinel-coded onsets silently biased estimates under `safety(off)` (mint staleness)

- **Symptom / impact:** with onset variables carrying sentinel codes (e.g. 0/9999 for
  never-treated) and `safety(off)` — exactly what the F13 error message tells such
  users to specify — old's row-level no-pre-period drop could strand surviving units
  without their reference period; stale `mint`/`mintp1` then made a treated dummy
  silently drop: measured 18–29% biased att with fabricated pre-trends up to sign
  flips, rc 0, alarms disabled by the same `safety(off)`.
- **Fix + proof:** `0c03010` (2026-07-27): unit-level drop + post-drop recompute.
  Parent probe: att 1.4234 vs truth 1.7436. Verifier: partial contamination bench-exact
  where parent was biased mreldif .29; user-if twin closed (parent mreldif .43).
- **This run:** rows `A15-mint-s/m/e` on `trig_mintsent_off` (blind-spot config): new
  correctly drops the sentinel-coded units (N 296 → 290, Nclus 30 → 29, e(sample)
  digest shifts accordingly); every fit scalar moves. The default-safety companion
  `trig_mintsent_def` (0 → 459) is the F13 refusal — see table 1.16.

### 1.7 `timetrends` numerically inert at large calendar codes (tcen)

- **Symptom / impact:** with raw large time codes (onset measured as low as t ~ 5e5),
  the demeaned `i.i#c.t` trend columns lost all precision and the requested unit
  trends went **numerically inert** — rc 0, no caution, results exactly equal to the
  no-timetrends run (bisect probe: reldif exactly 0). Users believed they had
  controlled for trends when they had not.
- **Fix + proof:** `f7019c9` (2026-07-27): mean-centered `gen double` trend injector,
  exactly shift-invariant. Suite F8/F9; battery 28/28; two verifiers (48/48
  non-timetrends configs BYTE-identical; timetrends configs numeric-only within the
  pre-registered footprint).
- **This run:** rows `A12-ttb-*` on `trig_tt_big` (blind-spot config at t ~ 5e5): old
  reports the no-trend fit (r2 0.8887, ll −116.40) while new fits working trends
  (r2 0.9069, ll −70.05, b worst 2.26e-01); the documented signature also fires —
  working trends omit e(V) on this design (`properties` b V → b) where old's inert
  trends kept it. Row `A12-ttp-s` on `trig_tt_pois` bounds the poisson-side IRLS
  reparameterization noise (this run worst 2.06e-09, bound 1e-4; run-varying class).

### 1.8 The adversarial pre-trend family: wrong fits, wrong reference, crashes — then a recalibrated statistic (B16 + 4483758 + NFADV)

Three coupled old-build defect layers plus one author-approved intended recalibration
stack on the same e() objects, so they are documented together (the ledger separates
them by row).

- **Old-build defects (B16 main half, commits `c0e17f8`, `05a243b`, hardening
  `44fcf47`/`291090c`, tranche 12):**
  1. *Cut off-by-one/two:* under `esfixedbaseperiod` with the reference below the
     reported window, the relyr −1 coefficient fell out of the pre-trend line fit and
     was tested as post-treatment (hand-verified from e() exports: wrong 2-point-line
     p 1.67e-25 vs correct 1.51e-22 on the esrel(-3)+espre(2) cell).
  2. *Reference-gap hunt:* old located the reference hole by scanning consecutive-t
     jumps — aborting r(198) on gap-year data, and with one wrong jump silently
     planting the reference at a DATA hole (parent probe on `adv_gapb`'s design:
     p 0.00892 where the correct reference gives 0.00184 — crossing 0.05 on variants).
  3. *Degenerate fits:* 1-row fits aborted r(2001); silent `finaladvtest`
     cancellations printed nothing.
- **SG half (`4483758`):** the two subgroup adversarial + altmethod + esf blocks
  crashed r(2000) whenever `esrelativeto()` placed the reference below the reported
  pre-window (empty `skiplocation` degenerated the fit restriction to zero rows); the
  identical non-SG run returned rc 0.
- **Intended recalibration on top (NFADV-01/02/05, commit `4fcf127`, author-approved
  2026-08-19):** the primary plug-in statistic was measured **over-rejecting ~1.4x
  nominal** under realistic heteroskedastic V (not reliably conservative; null Monte
  Carlo rejection 0.04 recalibrated vs 0.90 old), and the altmethod forecast statistic
  was anti-conservative (corrected covariance propagates A·V_PP·A′ + cross terms).
  New statistic: calibrated minimum-distance composite-null p on the event-study
  covariance; the line is republished as `slope`/`intercept_rt0` with the **slope
  bitwise-invariant** and the intercept re-expressed at relative time 0 (NFADV-05).
  Pooled-reference and degenerate-covariance p-values are no longer posted (loud
  refusals replace manufactured values — the parent manufactured p=7.2e-4 on a q=0
  cell). Blind post-commit verification: AGREE (three-route reconstruction ≤ 2.5e-12).
- **This run:** rows `A13-*` and `NFADV-*`:
  - `adv_altb` (`A13-altb-p/l/ln`): p **6.07e-23 → 3.69e-06** — the pre-registered
    B16+NFADV-02 movement exactly.
  - `adv_gapb` (`A13-gapb-p/l/ln`): p **0.00892 → 0.00274**. Note: the B16 planting
    fix alone gave 0.00184; the shipped value differs because the NFADV-01
    recalibrated statistic replaced the plug-in test on top of the corrected
    reference — both movements are pre-registered in the ledger row.
  - `esf_adv` (`A13-esfadv-*`): p 3.53e-10 → 8.29e-11 (B16 + svmat `6267522` +
    NFADV-01 stack).
  - `adm_esf_adv` (`NFADV-adm-*`): p 2.52e-11 → 2.98e-12 on the realistic panel.
  - `NFADV-ext-p/l/ln` on **145 ext configs**: recalibrated p-values, plus REMOVALS
    where old posted pooled-reference/degenerate-covariance p-values that are no
    longer defensible (e.g. `ext_b1_h11_f1`: `adversarial_pval_att` removed, was
    2.86e-04; `ext_b1_h22_p`: removed, was 1.13e-40).
  - Line matrices everywhere: slope cells bit-identical or storage-noise; intercept
    re-expression moves the second cell by design (e.g. `esf_adv` −1.0243 → 0.4876).
  - Crash cures pinned as rc changes: `adv_gapa` 198 → 0, `adv_r2001` 2001 → 0,
    `trig_sgadv_deeppre` 2000 → 0 (table 1.16).

### 1.9 IM degeneracy guard: absurdly confident p on degenerate per-cluster effects (16263da)

- **Symptom / impact:** with `impvalues` on noiseless/near-deterministic outcomes, old
  reported a collapsed ~0 SE and an absurdly confident IM p (se ≈ 2e-16 → p ≈ 1e-227
  class) — a storage artifact, dangerous to placebo/simulation users.
- **Fix + proof:** `16263da` (2026-07-22), author-requested guard: when the IM ttest
  SE is negligible relative to the estimate (`se < 1e-10*|b|`), SE/t/p/CI are set to
  missing (point estimate kept) with one Caution; because `ereturn post` rejects a
  VCE with a missing diagonal, e(V) is dropped for such runs (`properties` b V → b).
  Verified 3/3 with 22/22 realistic configs bit-identical and zero false fires.
- **This run:** rows `A16-imd-*` on `trig_imdegen` (fully degenerate design): new
  drops e(V) entirely (old posted a zero-diag V), `wdidmainresults` se/t/p/CI cells
  0.0 → missing, b differs only at svmat storage noise (worst 3.01e-08, row
  `A16-imd-b`).

### 1.10 Silent b-only posting when some SEs are missing: the e(V)-preserve fix (47defa0)

- **Symptom / impact:** whenever any estimand's SE came back missing — a far
  event-study horizon reached by a single cluster under `impvalues`, or plain
  continuous runs with missing `*_ame` SEs — old silently posted **b-only** (no
  e(V), `e(properties)` = "b") with **no caution**: lincom/test/esttab downstream
  simply lost the VCV for every healthy estimand too.
- **Fix + proof:** `47defa0` + doc `c53e25b` (2026-07-26): degenerate rows get a
  negative sentinel (−1e100) V diagonal that renders as missing on every consumer
  surface tested (esttab/lincom/nlcom/test), healthy rows keep their true VCV, and a
  Caution reports the affected count. The wave proved zero-diagonal and
  identified-rows-only alternatives actively harmful. Verified 3+3 (32/32 clean
  configs bit-identical).
- **This run:** rows `A17-*`: `trig_imdeg_part` (partial IM degeneracy: `properties`
  b → b V, `wdidmainresults` rows restored, b at amplified svmat noise worst 2.41e-06
  within its 1e-4 bound), `cont_sg_tidose`/`cont_sg_gap` (missing per-SG VCV
  diagonals), `ext_cont_ext_untreatedzero` and `ext_path10_cont_constant_dose`
  (plain-continuous missing *_ame SEs) — on all of them old silently returned b-only,
  new posts e(V) with the sentinel plus Caution.

### 1.11 Per-subgroup stacked results: silent label shift via Stata's rowname padding (dupcol)

- **Symptom / impact:** when an early per-SG "pre" round appended a row but a later
  one failed (e.g. a time-invariant dose making the continuous pre-margins all-
  omitted), old's namelist ended up shorter than the matrix and **Stata silently pads
  a short namelist by repeating the last name**: every later column wore the next
  quantity's label — `_b[ES_att_SG0_post1]` returned the post0 estimate; up to 13
  trailing duplicate labels; the displayed SG tables contradicted their own relyr
  column. Values were positionally correct; the labels lied.
- **Fix + proof:** `9d2e4df` (pins `0b12ee1`, tranche 12): per-(round x object)
  append record; names can never exceed rows. Verifier: parity strictly
  matnames-only on 7 dumps; 130+ ES entries matched the independently-labeled
  `e(b_es_*_s<sg>)` at reldif 0.
- **This run:** rows `C-SGlab-b/s/b2/s2` on `cont_sg_tidose` and `cont_sg_gap`:
  e(b)/e(V)/`wdidsgresults_*` colnames show old's trailing `ES_ame_SG2_post1`
  repeats replaced by the correct `pretreatment_SG#_att` labels; **zero value cells
  differ** in these groups.

### 1.12 Whole-sample joint pre-trend tests ignored the requested unconditional VCE (E1b)

- **Symptom / impact:** the two binary whole-sample pre-trend joint tests hardcoded
  `df()` and silently used the conditional VCE even when the user specified
  `unconditionalse` — the reported joint p was computed under the wrong VCE.
- **Fix + proof:** `7ea1533` (2026-07-26, separable fold-in of the E1 chain): threads
  `vce(unconditional)` via a dedicated `jointmarginsdof` local; acceptance suite
  `test/test_dffix.do`.
- **This run:** row `A18-jt` on `trig_uncond_jt`: `joint_horiz` 0.0994 → 0.1270 —
  the test now uses the unconditional VCE the user asked for.

### 1.13 v1.3 returned missing pairwise-subgroup cells that HEAD recovers (DDG)

- **Symptom / impact:** on RCS+pweight/twoway-cluster+subgroups designs, old's
  double-margins pairwise machinery returned MISSING `e(pwsgs)` cells and missing
  `joint_*_allsgssame` scalars (degenerate/rank-deficient contrast tests), plus one
  garbage negative SE (−8.9e-16). HEAD's E1 posted-pwcompare + rank guard (commits
  `cd8ce95`/`2a1c9ca`/`7ea1533`/`184a6ff`/`ed780ac`, 2026-07-26) computes valid
  values. (Framing corrected in the 2026-08-18 report, ticket 1: the old cells are
  the %21x encoding of Stata missing, not maxdouble artifacts.)
- **This run:** rows `DDG-pw-1..4`, `DDG-jt-1..3` on the four
  `ext_cov10_rcs_*` configs: e.g. `ext_cov10_rcs_pweight_ddd_c1`
  `joint_att_allsgssame` **. → 6.03e-204** and 41 pwsgs cells missing-on-old →
  valid-on-new (e.g. [1,2] . → 3.29e-10); the garbage SE cell is a clean 0 on new.

### 1.14 Float-precision storage defects: `gen` float truncation (B21) and `svmat` float truncation (6267522) — the bounded storage-noise tiers

- **Symptom / impact:** old stored engineered-collinearity copies, interactive-control
  copies, estimand weights (B21), and every `svmat`-materialized intermediate (IM
  machinery, adversarial fit inputs, plot frames) as **float**. Consequences ranged
  from benign last-digit noise to one proven silent-garbage mode: a double-typed
  interactive control with level large relative to its within-unit variation
  **silently produced absurd coefficients with no e(V) at rc 0** on old, because
  float quantization defeated the engineered exact-collinearity drop (B21 assessment,
  two independent verifier reproductions).
- **Fix + proof:** `3cb0c41` (+ docs `7672719`) — seven `gen double` conversions;
  `6267522` — all 16 svmat sites and 4 predict calls to double. New values are the
  correct ones; movement is confined to documented classes (~1e-8 rel gen-double,
  ~4e-8 rel svmat estimates with tail-amplified p; IC regressor copies up to the
  measured 1.33e-4 ceiling).
- **This run — every bounded tier fired inside its bound** (bounds are ~10x measured
  ceilings; full list in the verbose log's tail):
  - `XIC-*` (IC-copy class, bound 1e-3): worst 1.33e-04 on
    `ext_bin_ext_interactive`/`ext_cont_ext_interactive`.
  - `A9-svmat-*` on `impv` (worst 2.12e-07 vs 1e-6 bounds), `A9-adm-*` on `adm_impv`
    (IM t/p tail-amplified to 1.66e-05 vs 1e-4), `A16-imd-b`, `A17-part-b`.
  - `A8-B21-*` on `allest_aw`, `A8-sgit-*` on `sg_attitime` (joint scalar
    tail-amplified to 3.30e-06 vs 1e-4), `B21-cwm-*` on `cov_cwm`, `B21-adm-*` on
    `adm_aw` (worst 3.13e-08).
  - `X-est-*`/`X-fit-*`/`X-joint-*` across the ext estimation family (106 configs):
    worst 9.46e-05 on `ext_c1_h22_p` b (elasticity-transform amplification, bound
    2e-4); fit scalars at ≤ 2.12e-07.

### 1.15 `_ppmlhdfe_d` left in the user's dataset on every clean hdfe-poisson run (B12 + S6)

- **Symptom / impact:** the cleanup condition was inverted (old line 3440), so every
  clean poisson run left ppmlhdfe's internal `_ppmlhdfe_d` variable in the user's
  data — and the old condition **actively deleted a user's own pre-existing variable
  of that name** on the plain-poisson path.
- **Fix + proof:** `9d9d26f` (2026-07-20; condition mirrors the creation gate) plus
  S6 `edebab1` (2026-07-27; frame-qualified epilogue drop). Verified 3+3 including
  the rc-3200 fallback path. Under `verbose` the variable is deliberately retained
  (author directive, B18/22 entry).
- **This run:** rows `E-B12`/`E-B12e` on 6 configs (`pois`, `trig_tt_pois`,
  `adm_pois_pw`, `ext_contract10_poisson_converged`,
  `ext_contract10_ppml_d_collision`, `ext_path10_ppml_separation`): the leaked
  variable is absent from new's post-run varlist and e(sample) digest.

### 1.16 Error-contract repairs: crash cures and informative refusals (the 31 designed rc changes)

All 31 differing-rc configs are pinned (`rc-by-design`) and each cites its fix. Two
buckets here; the third (new-only options, 198 → 0) is intended-change material in
section 3.4. One pinned rc change that is NOT a fix — `ext_contract10_rcs_factor_fe`
(0 → 198) — is a known naming-fragility defect and is documented in section 2.2, not
here.

**Bucket A — old crashed; new completes correctly:**

| config | rc old → new | fixed defect |
|---|---|---|
| `cont_uz_ss` | 198 → 0 | B5 `4979abf`: contreat+untreatedZero+summarystats crashed r198 ereturning an undefined local after estimation completed |
| `adv_gapa` | 198 → 0 | B16 `c0e17f8`: adversarial reference-gap hunt aborted on gap-year data |
| `adv_r2001` | 2001 → 0 | B16 `05a243b`: degenerate 1-row pre-trend fit aborted; now clean skip + Caution |
| `trig_sgadv_deeppre` | 2000 → 0 | `4483758`: SG adversarial skiplocation crash (empty-string fallback) |
| `trig_bigcodes` | 461 → 0 | ttre double storage `5a97f70`: YYYYMMDD onset grids aborted r461 (the silent cohort-merging variant, att bias −11%/−64%, was proven in the `8c5bfd7` campaign work) |
| `scr_badif` | 411 → 0 | B23 `8d7bd01`: input asserts ran over the whole dataset — old false-fired r411 on bad rows the user excluded via `if` |
| `trig_esfsup_ref322` | 322 → 0 | G10B `7ba1e16`: under DEFAULT safety, the unsupported-reference cohort made old abort with an opaque r322 after full estimation setup; new drops the refless cohort gracefully with a Caution (the safety(off) twin is section 1.17) |

**Bucket B — old ran silently wrong or failed cryptically/late; new refuses informatively:**

| config | rc old → new | fixed defect |
|---|---|---|
| `scr_uncond_alone` | 322 → 198 | B15 `14ca983`+`657f6f0`: unconditionalse without robust/cluster wasted full estimation then bare r322; now early clear r198 |
| `scr_cfxvar` | 0 → 198 | B17 `e4e3d1e`: a cfxplottypes token naming an existing variable was silently injected — old completed rc 0 with a mislabeled, meaningless plot |
| `scr_nonintsg`, `ext_contract_fractional_subgroup` | 198 → 411 | B20 `7f7f73c`: non-integer subgroups died cryptically (`_joint_horiz_i_s1.5 invalid name`); now clear r411 |
| `scr_lat0` | 0 → 198 | latticenum `d1f7542`: old silently remapped `lattice(all 0)` to the gquantiles default (bit-identical to `lattice(all 2)`) |
| `scr_latnocon` | 0 → 198 | `4c8fb66`: `lattice()` without `contreat()` was parsed then wholly ignored, no message |
| `scr_fracttre`, `ext_contract10_fractional_time` | 0 → 459 | S5 commensurability screen (`905617d`..`f187d7d`, band `3151f6e`): in-window fractional onsets ran rc 0 with the onset period mislabeled as pre-treatment (fractional_time attributed by screen class; its config note is unpinned) |
| `ext_val10_fractional_ttre` | 452 → 459 | same screen: old failed late/opaque ("factor variables may not contain noninteger values"); new errors early with the fractional-grid hint |
| `ext_contract10_inconsistent_ttre`, `trig_mintsent_def` | 0 → 459 | F13 `c4d7e62`: ttre varying within `i` (incl. onset-recorded-only-post coding) ran silently — measured +41% horizon bias, sign-flipped manufactured pre-trends; now informative r459 naming the offending unit |
| `scr_sgsplit`, `adm_sgsplit` | 0 → 459 | F13 `c4d7e62` + `e5445c4`: units split across subgroups made SG contrasts structurally meaningless (old returned the pooled estimate under every label); now refused (T-B in REPORT_2026-08-18) |
| `ext_contract10_if_outside_negative` | 411 → 322 | B23 `8d7bd01`: old's raw-domain assert false-fired r411 on negative-y rows excluded via `if`; new correctly runs the restricted sample, which then fails r322 for the genuine downstream reason |
| `ext_contract_bad_elasticity`, `ext_contract10_bad_safety`, `ext_contract10_poly_zero` | 0 → 198 | `b42b71f` VAL-G19/G21: silent-accept option domains closed (the three rc-pin updates of this run, section 0). Verified old-side semantics: an unrecognized `contreatelasticitytype()` was silently rewritten to `dydx` (bit-identical b/V to an explicit dydx run); `contreatpoly(0)` silently ran as order 1 (bit-identical to `contreatpoly(1)` — old only appends terms when poly > 1); any non-`off` `safety()` token, garbage included, silently meant the safety-ON default |

### 1.17 The esf calendar-support drop ran AFTER the reference checks: silent-wrong estimates under `safety(off)` (G10B)

- **Symptom / impact:** in the `esfixedbaseperiod` branch, old ran the per-group
  reference checks BEFORE the calendar-support drop (periods with no control units).
  When the support drop deleted a calendar period containing a cohort's fixed
  reference row, that cohort survived with no reference: an opaque r(322) under
  default safety, and **silently wrong estimates at rc 0 under `safety(off)`** —
  exactly the switch the error message steers users toward.
- **Fix + proof:** `7ba1e16` (VAL-G10B, 2026-08-19): minimal support-first reorder,
  fixed-point proof, two review rounds; benign esf configs bit-identical across
  builds; HEAD-side contract pinned by `test/test_esfsupport.do`.
- **This run (coverage added in-run; deterministic 14-unit fixture, controls
  unobserved at t=3, cohorts onset 4 — reference t=3 unsupported — and onset 5):**
  a three-value triangle proves the old value WRONG, not merely different:
  - OLD, full fixture, safety(off): rc 0, att = **0.6149** — silently wrong,
    understating the DGP effect (1.5) by ~56%.
  - OLD, hand-benchmark sample (refless cohort removed by hand): att = **1.4106**
    — old's own code produces the correct value when handed the correct sample,
    isolating the sample-construction bug.
  - NEW, full fixture, safety(off): att = **1.4106**, equal to the benchmark to all
    16 printed digits.
  Mechanism visible in the dumps: old's e(sample) keeps 70 rows (the refless
  onset-4 cohort retained) vs new's 50; new's cohort drop also removes the ES pre2
  cell from b/V/wdidmainresults (1x6 → 1x5), which is why the `G10B-esf-n`
  matnames row exists. Configs `trig_esfsup_ref` (ledger rows `G10B-esf-s/m/e/n`,
  accepted blind-spot class) and `trig_esfsup_ref322` (default safety, 322 → 0,
  table 1.16 bucket A).

---

## 2. Old-build issues NOT visible as old-vs-new mismatches

These are defects or environment facts that the dual-build comparison structurally
cannot flag — both builds agree (wrongly), or the difference is environmental. They
are recorded so "0 unexplained" is not misread as "old had no other problems."

1. **`ccc_absorb()` broken on BOTH builds** (`cov_ccc_absorb`, both r198): the
   regression completes, then a bare "invalid syntax" r198 in the margins phase.
   Zero coverage existed anywhere before this suite. Standing defect ticket
   (REPORT_2026-08-18 ticket 2); the mini-assessment's T-M1 ruling (ignore, with the
   ccnoabsorb corner noted) is adjacent.
2. **`o.`-omission-marker naming fragility** (`ext_contract10_rcs_factor_fe`,
   pinned 0 → 198 with note "5o. invalid name"): new errs under Stata 19 on
   RCS+factor-fe; the colleagues' Stata 17 **old** build hit the same error (their
   `msg_old` records it); our Stata 19 old build escapes by omission-pattern luck.
   The pin keeps the run green — it does **not** certify the behavior as intended
   (REPORT_2026-08-18 ticket 3). The dupcol entry (section 1.11) notes the related
   `o.`-marker truerow-strip trigger as its own future ticket.
3. **Unit-level `i` at realistic scale exits r(907)** (ticket T-A): 2,000 counties x
   periods exceeds maxvar via the cell-dummy grid; `adm_unit_i` runs a ~280-county
   subsample instead. Correct usage at that scale is cohort-level `i`, but the
   failure is cryptic and late; a pre-flight size hint remains ticketed.
4. **Stata-17 environment-dependent rc pins** (REPORT_2026-08-18 ticket 5): the
   colleagues' bundle records different rcs on Stata 17 (ppmlhdfe separation 503 vs
   0 here; the rcs_factor_fe case above). Pins are calibrated to THIS environment
   (Stata 19 MP) with notes in `configs_ext.csv`.
5. **Corrected T1 reading** (REPORT_2026-08-18 ticket 1, corrected framing): the
   old-build missing pairwise-SG cells recovered by HEAD (section 1.13) were
   initially misread as maxdouble artifacts; they are the %21x encoding of Stata
   missing. Recorded here because the correction is part of the report of record.
6. **T-M2: continuous-SG joint pre-trend p stored as literal missing on BOTH
   builds** — the sg-continuous overjoint `r(p)` is 1x1 with a real p at [1,1]
   (instrumented probe: sg1 6.209e-06, sg2 4.109e-07) that the `[1,2]` read silently
   discards; `e(joint_Choriz_s*)` stores "." where significant p-values exist. Old
   shares the bug, so this comparison can never flag it. Author ruling: DEFER to the
   full assessment under the mandated evidence protocol
   (`review_process/mini_assessment_2026-08-23.md` section 2b); `cont_sg.csv` is the
   sole baseline pinning the ".". Adjacent open observation (solved_issues, T-M2
   final ruling): the continuous SG-equality (`allsgssame` ame) joint tests are also
   missing-valued in both builds and pinned missing — separate site family (E1 chain),
   undiagnosed, deep-phase item.
7. **Remaining validation gap running rc 0 silently on both builds**:
   `poisexpresults` without `poisson` (`ext_contract10_poisexp_without_poisson`,
   judged EXPLAINED this run via CORE10-cmdline only — i.e. both builds still accept
   it). The other two members of REPORT_2026-08-18 ticket 4
   (`contreatelasticitytype`, `safety()`) were closed by `b42b71f` and are the rc-pin
   updates of section 0.
8. **Both-error parity pins (42 configs)**: identical error rcs on both builds,
   pinned by return code only (messages out of contract) — e.g. `scr_strclus` r109
   (old bare type mismatch late, new clean at parse — same rc), `trig_pois_sep` r322
   separation parity, `scr_impsg` r198 (impvalues+subgroups refused by both), the
   `ext_val_*` validation deck. These verify old and new agree on the error contract
   where no fix intervened.

---

## 3. Intended changes (non-bug divergence classes)

### 3.1 `e(cmdline)` begins with the command name (CORE-10)

Commit `aec904e` (isolated re-bless, verified row-by-row). Standard Stata convention;
old omitted the leading `wooldid `. Ledger row `CORE10-cmdline` fired on **all 375
rc-0 configs** — a uniform 8-byte prefix, no other content change. This is the sole
reason MATCH=0 this run (section 0).

### 3.2 Stable cfx row keys (NFFAST-02)

Commit `cb1996c` batch. Old cfx matrices carried margins-tempvar row stripes
(`intensity:2.__TMP1`, `1.__TMP1:3.__TMP2`) — volatile every run and meaningless to
consumers. New uses stable `lattice_#` keys with null roweq. Ledger row
`A6-cfxstripe` fired on 6 configs (`trig_cfx`, `plot_cont`, `ext_char_cfx_nonatt`,
`ext_cont_ext_cfx`, `ext_cont_ext_cfx_subgroups`, `ext_cov10_cfx_inference`);
matnames-only, numeric cells covered by their own rows.

### 3.3 saveplots x subgroups filename scheme

Commit `06e1a5b`: the per-subgroup ES data filename now separates the subgroup
variable and value (`svpadm_ES_region1.dta` → `svpadm_ES_region_1.dta`), removing
ambiguity for multi-digit subgroup values. Ledger row `C-esfn` on `adm_svpsg`
(`files|` comparison; .gph names unchanged).

### 3.4 New-only options (old r198 unknown-option, new rc 0)

Seven rc-by-design pins exist solely because the option does not exist on v1.3:
`gt_basic` (gteffects, `b65b91b`), `perm_basic` (permute, F6 series `125bbd0`..),
`scr_antic` (anticipation, `1f521f6`), `scr_fastnose` (fastnose, `fc0cadb`),
`scr_ptp` (pretrendpower, `f9fd239`), `scr_permexact_t`, `scr_permstrata` (permute
sub-options). Old's r198 is correct old behavior, not a defect.

### 3.5 Adversarial statistic recalibration and line coordinates

Author-approved intended changes NFADV-01/02/05 (`4fcf127`) — statistic replacement,
covariance correction, `slope`/`intercept_rt0` re-expression, and unposting of
pooled-reference/degenerate-covariance p-values. Documented with the old-build
defects they ride on in section 1.8.

### 3.6 e(V)-preserve posting behavior

The `47defa0` sentinel-V design (section 1.10) is simultaneously a bugfix (old's
silent b-only) and a deliberate contract change: partially-degenerate runs now post
e(V) with sentinel diagonals plus a Caution instead of posting nothing. Recorded in
both places by design.

### 3.7 Message texts and console output

Message texts changed program-wide across the campaign (Error:/Caution: conventions,
S5/C12; mini-assessment `85da19b`). Console output and error message texts are out of
comparison contract (README, "both-error parity pins assert the return code only");
no ledger row is needed or present for them.

### 3.8 New-only e() additions (11,762 across the run — informational, never diffs)

Additions never reach the comparer's diff path; the classes: F2 event-study VCV
returns (`b_es_*`, `V_es_*`, `numpre_*`, `d4b2118`), F8 composition diagnostics
(`escomp*`, `escomp_ref*/trt*`, `d87adfd`), F11 replay state (`wdid_*` hidden
macros, `d2a398e`..), permute/gteffects/ptp returns on their configs, and NFES-02's
`es_ndegen_*` degenerate-row counts (ledger comment block documents
`ext_cont_ext_untreatedzero` gaining `e(es_ndegen_ame)` = 2).

---

## 4. Residual notes

- **Accepted blind-spot configs (5, by design):** `trig_tt_big`, `trig_allest_drop`,
  `trig_extmiss`, `trig_mintsent_off`, `cov_kitchen` carry unbounded `*` ledger rows —
  they exist to prove a documented LARGE fix fires on a total-divergence design, not
  to bound its magnitude (README, "Accepted blind spots"). Every other row is bounded
  or named-object-scoped.
- **Run-varying ULP classes (bounded, documented non-reproducible):**
  `ULP-ppw-*` (`adm_pois_pw`, ppmlhdfe+pweight IRLS; this run worst 3.21e-09),
  `ULP-upw-*` (`adm_uncond_pw`, unconditionalse 1–2 ulp class; worst 2.92e-07),
  `A12-ttp-s` (timetrends-poisson IRLS; within 1e-10 one run, 1.24e-5 another),
  `X-est-pw` (degenerate-SE pwsgs cells; observed 2.2e-7 and 1.12e-4 across runs).
  All are non-reproducible on BOTH builds run-to-run; bounds calibrated across runs.
- **Accepted blind-spot list grew to 6:** `trig_esfsup_ref` (section 1.17) joins the
  five above; README updated.
- **Coverage gaps (prereg P3; gap 1 was closed in-run by section 1.17):**
  1. The NFADV pooled-reference/degenerate-q unpost paths are ledgered only for
     `ext_*` configs (`NFADV-ext-p` removals); no core/adm config hits those paths.
  2. `cov_icbad` under-triggers the B21(b) engineered-collinearity class (levels
     ~1e6 are float-representable; a ~1e8 variant would pin it old-side).
- **Comparison-scope limits (unchanged from certification):** console text never
  diffed; `.gph` compared by name only; both-error pins compare rc only; the
  colleagues' oracle/metamorphic library is not ported (future component).
- **Suite changes landing with this report:** the three rc-pin edits
  (`configs_ext.csv`), the two G10B configs + prep (`configs.csv`, `preps.json`),
  the four `G10B-esf-*` ledger rows, and the README blind-spot-list update.

---

*Report of record for the 2026-08-23 assessment. Every claim above carries its
commit sha, solved-issues entry, ledger row id, or config name; the underlying
verification narratives (assessor/verifier waves, parent probes, battery evidence)
are in `wooldid_solved_issues.md` and the campaign ledger.*
---

# 5. New-feature statistical validation (2026-08-24 extension)

The assessment's second phase (work-file chunk "New Features") validated the surfaces
with no old-build analogue by EXTERNAL and STATISTICAL means: preregistered
Monte-Carlo size/power calibration (`test/suite/PREREG_NEWFEAT_2026-08-24.md`,
committed 4a39dc1 BEFORE any run; MC component `test/suite/mc/`, 29d3d26), external
oracles (Roth's `pretrends` package v0.5.1, honestdid v1.1.2, a reghdfe common-onset
oracle), and deterministic parity grids. Four new permanent suite scripts (f8ac51d,
110f8f5) bring the `newfeat` component to 7/7 PASS; five new battery configs pin the
previously-unpinned new-feature surfaces at rtol 0 (c5d00d2). Every MC rep records
its seed; each row reproduces standalone.

## 5.1 VALIDATED — clean bills of health

- **Randomization inference (`permute()` family): exactly sized, 7/7 preregistered
  configs PASS.** M=400 each, B=99 (or exact enumeration): shuffle baseline
  (18/400 rej at nominal .05, band [11,29]); AR(1) rho=.7 within-cluster errors
  (18) and 1..4x heteroskedastic clusters + permtstat + permtci (18) — the
  cluster-asymptotics-hostile cases where RI's finite-sample exactness is the
  selling point; strata (22); probabilities with its 3^10 assignment space (20;
  inestimable-draw exclusion behaves as the conditioning argument requires);
  permexact's 20-assignment enumeration (20/400 = .0500 on the nose); analytic
  weights (23). Secondaries all clean: the 5-95 null band covers the comparator
  356-361/400 (band [348,372]) in every config; the inverted-t permtci covers the
  true zero effect 382/400 (band [371,389]); grid-uniformity chi-square flags
  nothing. The external draw-stream re-implementation (`ri_crosscheck.do`, 5 blocks)
  passes on the same build.
- **pretrendpower power claims: PASS both arms.** Planting a linear pre-trend equal
  to the run's reported mdt_wald yields 83.8% rejection of the pre-trend Wald test
  (band [.72,.87]); planting mdt_nis yields 82.0% for the any-pre-significant (NIS)
  test. Roth's own `pretrends` package, fed wooldid's e(b_es)/e(V_es)/e(numpre),
  agrees: conventions match at 1e-12, mdt_nis vs the package's slope-for-power at
  reldif 8.2e-4..2.1e-3, package power evaluated at wooldid's mdt_nis(80) =
  .809..816. (The package reports no data-fitted slope, so wooldid's GLS
  e(roth_slope) has no package analogue to compare — structural N/A.)
- **F2 event-study b/V returns: externally reproduced at machine precision.** A
  common-onset reghdfe TWFE event study (single-treated-cohort formulation)
  reproduces e(b_es_att) at mreldif ~1e-16, e(V_es_att) BITWISE on 3 of 4 variants
  (aw, asymmetric window, poisson via ppmlhdfe included), and the att aggregation
  (treated-obs/aweight-sum weighted post average) at ~1e-16
  (`newfeat/reghdfe_oracle.do`).
- **honestdid handoff: validated on an 11-case grid** (binary base, anticipation,
  asymmetric windows, aw, per-subgroup ES, admin panel, rich datasets, dta_cont),
  delta(rm) AND delta(sd): the Original row reproduces the wooldid post0 CI at
  reldif 0.0 in 22/22 extractions; layout, mvec, and rm-monotonicity clean
  (`newfeat/honestdid_grid.do`).
- **fastnose: full parity.** 22-config grid (subgroups, poisson, contreat, weights,
  anticipation, estimands, lattice, gteffects): zero violations of the two-arm gate
  (|diff|<=1e-12 OR magnitude-rel<=1e-10); worst abs 7.1e-15, inside the documented
  fast-path envelope. All 7 incompatible-option refusals pinned at rc 198 with the
  fastnose-specific screens verified (`newfeat/fastnose_grid.do`).
- **anticipation / escomp / gteffects / replay**: shift-identity, recompute-from-data,
  aggregation-identity, and byte-replay checks all pass (existing identity_extras +
  cert suites, re-run green this chunk); new battery pins cover fastnose+gteffects,
  pretrendpower, anticipation, permtstat+permtci, permute(probabilities).

## 5.2 NF-1 (HEADLINE ISSUE): V_es-based joint tests are severely oversized under
## cohort-level clustering — adversarial, jointtests, and the roth pre-trend Wald

**Finding (preregistered FAIL class).** With `cluster(i)` at the cohort level, the
flat-null rejection rate of `adversarial` at nominal 5% is:

| design (all null true) | adversarial @.05 | joint_horiz @.05 |
|---|---|---|
| G=40, singleton cohorts, 4/4 window (A1, M=1000) | **20.5%** | 11.3% |
| G=40, planted-line null (A2) | **23.1%** | (33.8% = genuine power, its null differs) |
| altmethod, flat / planted (A3a/A3b) | **12.6% / 11.2%** | 12.3% |
| G=40, AR(1) errors (A6) | **23.2%** | 15.2% |
| G=200, singleton (D1, M=240) | **11.7%** | 8.7% |
| G=40, 8 units/cohort (D2, M=250) | **48.4%** | 36.0% |
| G=100, 8 units/cohort (D4, M=250) | **44.8%** | 32.0% |
| G=40, singleton, 2/2 window (D3) | 6.4% | 6.8% |
| G=40, 8 units/cohort, 2/2 window (D5) | **24.8%** | 21.2% |
| G=10, singleton (A5) | 2.4% (conservative) | **65.7%** |

The roth pre-trend Wald (`e(roth_prep)`) shows the same inflation (14.0% at .05 on
the ptp design's gamma=0 arms), and the reported NIS size understates the empirical
any-pre-significant rate (reported mean .1325 vs empirical .2132, n=408 — outside
the binomial band). Chi-square grid-uniformity is destroyed in every failing arm.

**Verified mechanism (direct MC measurement, V1/V2/D5 arms + INV-1 report
`review_process/inv1_v_es_ame_2026-08-24.md`).** The Mundlak specification gives
every treated cohort x relative-period cell its own coefficient; under cohort-level
clustering each cell sits inside ONE cluster, so the OLS normal equations zero the
within-cell residual moments and the plug-in CRVE score of every cell coefficient
vanishes. The damage is DIRECTIONAL: per-horizon reported variances are only mildly
low (Var_true/mean reported diag = 1.0-1.2 singleton; 1.5-2.3 multi-unit), but the
covariance's small-eigenvalue directions — exactly the cross-horizon CONTRASTS the
line-fit residual lives in — claim ~7-20x less variance than the truth (mean
reported smallest eigenvalue .0206/.0014/.0033 in V1/V2/D5 vs true adjacent-contrast
variance ~.20/.026/.024). The whitened minimum-distance statistic explodes in those
directions. Singleton designs look "better" only because their EXACTLY-zero
directions are dropped by the 1e-12 whitening tolerance (and at G=10 enough drop to
make the F(q, 9) reference conservative); multi-unit cells leave small-but-wrong
eigenvalues in play — the worst case (48-60%). Larger G barely helps (D1 ~10%,
D4 ~42%): this is not a finite-G asymptotics defect but a structurally-blind V.
The test CODE is sound: under `robust` on iid multi-unit data the same adversarial
test is correctly sized (A4 att: 5.75% at .05, PASS; joint_horiz 5.95%).

**Old-build status:** the specification, clustering, and V are identical in old
wooldid (INV-1 probe8); the oversizing is a SHARED long-standing property, newly
QUANTIFIED here, not a regression. The sthlp already cautions that these tests are
sensitive to the number of clusters; the true dependence is on cluster-saturated
cells and test dimension, and the magnitudes above are far beyond a rounding caveat.

**Disposition: DOCUMENTED, fix deferred to the bug-fixing chunk** (work-file order:
"we will work on solving identified bugs later"). Candidate remedies for that chunk,
in rising ambition: (a) sthlp rewrite stating the design conditions under which the
adversarial/joint p-values are trustworthy (clustering finer than cohort-by-period
cells, or robust/unclustered errors; narrow windows do NOT rescue multi-unit
designs); (b) detect cell-saturated clustering and print a specific caution on
every V_es-consuming p-value (adversarial, jointtests, pretrendpower, and the
honestdid handoff note); (c) a cluster-bootstrap or CR3/jackknife V for the ES
block; (d) eigenvalue-floor regularization inside wooldid_advstat (defensible but
ad hoc). Note remedy (b)'s detection predicate is cheap: every treated (i, relyr)
cell inside a single cluster of the specified cluster().

## 5.3 NF-2: rank-1 e(V_es_ame) under cohort-level clustering (defect-class);
## old build posted an INVALID adversarial ame p there — new build refuses

The extreme case of NF-1 (full report: `review_process/inv1_v_es_ame_2026-08-24.md`):
with `contreat()`, each treated cell carries its own dose slope while the pooled
baseline dose slope phi is absorbed; under cohort-level clustering the CRVE slope
block collapses EXACTLY to Var(phi_hat)*11', so e(V_es_ame) is rank 1, every ame SE
(main, pretreatment, all ES horizons) equals SE(phi_hat), and cross-horizon
contrasts carry zero reported variance while the true covariance is full rank and
~5x larger (120-rep MC). Consequences: honestdid/lincom/test consumers of
e(V_es_ame) silently ingest a meaningless covariance; the adversarial ame test has
no testable dimensions. **New wooldid correctly refuses** ("composite-null
covariance is degenerate; no p-value computed"). **Old wooldid POSTED
e(adversarial_pval_ame)=.6147 on the same degenerate V** via a generalized-inverse
test — an invalid p, i.e. the refusal is a genuine old-bug fix, invisible to the
old-vs-new comparison because no suite config combined cohort-clustered contreat
with adversarial (the ext_* configs that post ame p's are unclustered).
Unclustered/robust/finer-than-cell configs produce a full-rank, internally-sound
V_es_ame. Fix candidates fold into NF-1's; the F2 sthlp entry deserves the caveat
regardless.

**Runnable demonstration:** `review_process/demo_vesame_pathology.do` (from the repo
root: `"<stata>" -e do review_process/demo_vesame_pathology.do`; deterministic, ~15s).
On a 30-cohort x 6-unit irregular-dose null panel it shows, under `cluster(i)`: all
eight ES_ame horizons (and the main/pretreatment ame rows in the printed table) with
the IDENTICAL SE 0.059040286 while the estimates differ freely; min pairwise
correlation of e(V_es_ame) = 0.99999999999999; eigenvalues .0279 vs 6e-17 (rank 1);
the (post0 - pre2) contrast with reported variance -2.3e-17 against an estimate
difference of 0.491. The same data under `robust` gives horizon-specific SEs of
0.166-0.211 — roughly 3x the clustered figure, making the understatement visible
directly — and a full-rank V (eigenvalues .0245-.0605).

## 5.4 NF-3: ame-side tests remain oversized even under robust VCE

A4 (30 cohorts x 6 units, dose varying within cell, iid errors, `robust`): the att
side is correctly sized (5.75%), but `adversarial` on the ame event study rejects
**31.5%** at nominal 5% and the ame joint pre-trend test (`joint_Choriz`) **22.0%**
(M=504; the ame p posts every rep here — full-rank V). Root cause family: per-cell
dose slopes are estimated from few observations; their HC variance is biased down
(leverage), and the 6-8 dimensional quadratic forms amplify it. Continuous-treatment
ES inference is therefore oversized in EVERY tested configuration — degenerate
under cohort clustering (NF-2), inflated under robust. Same disposition as NF-1.

## 5.5 Documentation-grade notes and upstream findings

- **Singleton-cell V_es_att is a structurally different estimator** (oracle O16):
  with unit-level i (every treated cell a singleton), treated residuals vanish and
  V_es carries control-cluster noise only (SE ratios 0.84-1.07 vs the pooled event
  study; b and att still match at 1e-15). Candidate sthlp caveat on the
  honestdid/pretrendpower handoffs for singleton-cell designs.
- **honestdid C-LF robust confidence sets can be legitimately EMPTY** on data with
  sharply significant pre-trends (observed on binary_fine and dta_cont, pre-t
  -4.1/-5.5): a documented Rambachan-Roth outcome, not a handoff defect (Original
  row reproduces exactly on the same inputs). Encoded as an explicit allowlist in
  honestdid_grid.do (preregistered-deviation note in the script header).
- **Upstream honestdid v1.1.2 bug:** delta(sd) crashes with mata assert 3498 in
  _flciFindWorstCaseBiasGivenH on the binary_mixed_imbalanced case (parallel AND
  sequential paths); delta(rm) runs clean on the same b/V. Candidate courtesy
  report to the package maintainer.
- **Preregistration deviations, both documented in-file:** the M5 "every robust row
  finite" rule gained the legitimately-empty allowlist; A4 was fixed post-INV-1
  (before any A4 rep fired) to the sound robust-VCE class.

## 5.6 Where the artifacts live

Prereg + amendment: `test/suite/PREREG_NEWFEAT_2026-08-24.md`. MC drivers/analyzer:
`test/suite/mc/` (judged report: `mc/out/mc_report.md`; per-rep CSVs with seeds in
`mc/out/`, gitignored — regenerate with `./run_mc.sh all`). New suite scripts:
`test/suite/newfeat/{reghdfe_oracle,honestdid_grid,pretrends_crosscheck,fastnose_grid}.do`.
INV-1 mechanism report: `review_process/inv1_v_es_ame_2026-08-24.md`. Battery pins:
`test/baseline/{fn_gt,ptp_two,antic_es,perm_tstat,perm_probs}.csv`.

## 5.7 Ticket index — every documented-not-fixed item, for author consideration

All of these were deliberately left unfixed this chunk (work-file order: bug-fixing
follows in the next chunk); each cites its evidence above. Tickets NF-1a..NF-1c and
NF-2/NF-3 share one root cause (the V_es pathology) and one remedy menu (§5.2 a-d).

| ticket | surface | symptom (measured) | where |
|---|---|---|---|
| NF-1a | `adversarial` p (att) under cohort-level cluster() | 20.5-23.2% rejection at nominal 5% (G=40, incl. AR(1)); 48.4% with 8-unit cohorts; 44.8% at G=100; 24.8% even at a 2/2 window multi-unit; altmethod 11-13% | §5.2 |
| NF-1b | `jointtests` pre-trend p (`joint_horiz`; old feature, same V) | 11.3-36.0% at G=40 configs; **65.7% at G=10** (flat null) | §5.2 table |
| NF-1c | `pretrendpower` reported p/size scalars (`e(roth_prep)`, `e(roth_nissize)`) | prep rejects 14.0% at nominal 5% under the null; reported nissize .1325 vs empirical .2132 (power claims themselves PASS) | §5.2 |
| NF-2 (demo issue) | ALL displayed ame SEs + returned `e(V_es_ame)` under cohort-level cluster() | every ame SE identical (= SE of the absorbed baseline dose slope); V exactly rank 1; contrasts get zero reported variance; true variances ~5x, full rank; honestdid/lincom/test ingest it silently. Old build shares the V and additionally POSTED an invalid g-inverse adversarial ame p (.6147) where new refuses (verified in probe8_old.log) | §5.3 + demo .do |
| NF-3 | ame-side tests under ANY tested VCE (incl. robust) | adversarial-ame 31.5%, `joint_Choriz` 22.0% at nominal 5% on iid multi-unit data with robust (att side simultaneously correct at 5.75%) | §5.4 |
| (prior chunk, unchanged) | continuous-SG `joint_Choriz_s*` stored missing via the r(p)[1,2] read (T-M2) + the allsgssame-ame missing scalars | deferred 2026-08-23 by author ruling with a mandated >=3-dataset console+matrix evidence protocol; same jointtests-family neighborhood, untouched this chunk by design | §2; solved-issues ledger |

Documentation-grade notes (no wooldid defect; candidate sthlp caveats / upstream
reports) are in §5.5: the singleton-cell V_es structural property, honestdid v1.1.2's
delta(sd) mata-assert crash, legitimately-empty Rambachan-Roth C-LF sets, and the
pretrends package's absent data-fitted slope.

Explicitly NOT tickets: randomization inference, fastnose, gteffects, escomp, replay,
and anticipation came through the full validation with zero documented bugs (§5.1).
