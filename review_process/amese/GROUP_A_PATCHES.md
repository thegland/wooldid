# Group A Patches — the AMESE R1/R3/R4 implementation pool (2026-08-24/25)

The canonical handle for the family-A behavior change (author-ruled R1 + R3 + R4;
R2 skipped). Use this list to test, discuss, or revert the pool as a unit. Final
verified state: HEAD 294513a-era chain, oldnew 450/450 (prov new|294513a), cert
22/22, battery 34/34 rtol-0, newfeat 7/7.

## TESTING RULE (applies to this repo generally)

The default Stata adopath does NOT serve this repo's build: a copy of the LIVE
external wooldid (the version outside users access) sits both in the PLUS adopath
and at /Users/hegland/Documents/Research/wooldidinstall/ — BY DESIGN (it mirrors
what Stata loads by default and what external users run; it updates only when the
repo's edits go live). Therefore, when testing the development build, ALWAYS load it
explicitly: `do <repo>/wooldidinstall/wooldid.ado` before calling `wooldid`. Never
rely on adopath resolution, and treat any test result obtained without an explicit
load as suspect (this bit twice during the Group A wave: theory-agent probes and a
cert-suite path slip both silently ran the live build).

## Tier 1 — behavior-changing patches (the revert set proper)

| sha | subject | touches |
|---|---|---|
| 7054878 | Core R1/R3/R4: wooldid_amescreen subprogram (detect/vzero/sescan), containment predicate, ame-family scrubs at every capture site, F2 split (b_es_ame posts, V_es_ame withheld), ame-side test gates, R4 detection+refit, generic guards | wooldidinstall/wooldid.ado, wooldidinstall/wooldid.sthlp (2 caution paragraphs), test/test_amescreen.do (new cert), test/baseline/{cont_basic,cont_uz,cont_uz_ss,cont_lat,cont_ss,cont_sg,cont_sg_tidose,cont_sg_gap}.csv (re-bless), review_process/demo_vesame_pathology.do (header note) |
| ea69544 | R4 narrowing: refit gated on firing predicate + single-var cluster(), reordered after the predicate (fixes: binary-timetrends V resurrection, the rcs o.-crash side-cure, two-way-cluster 0->198 regression) | wooldidinstall/wooldid.ado |
| c545f03 | Verifier-wave fixes: glevelsof capture (rc-920 at scale), sescan differing-estimates signature, cluster-abbreviation parse (cl..cluster), cfx-ame matrix scrub + noci plots; ledger reason amendments + bounded R4-b rows | wooldidinstall/wooldid.ado, test/suite/ledger.csv |
| 294513a | V1 closing: refit barred under fe()/ccc_absorb (bare-factor linearization), caution texts corrected (byCohort ame-point normalization dependence; neutral vzero advice) | wooldidinstall/wooldid.ado |

## Tier 2 — suite-contract patches (revert together with Tier 1, or the suite fails)

| sha | subject | touches |
|---|---|---|
| cae6947 | 26 AMESE ledger rows scoped to the fire-list | test/suite/ledger.csv |
| c184ea0 | test_amescreen run-from-test/ cwd convention | test/test_amescreen.do |

## Tier 3 — records and evidence (KEEP on any revert; no behavior)

| sha | subject |
|---|---|
| f4aeb60 | R1 design note (pre-implementation contract) |
| a603f38 | T-M2 mandated evidence banked (pre-R3 build — deliberately interleaved in this range; NOT part of the patch set) |
| 42d12d9 | Three verifier reports preserved |
| 439d5cb | Campaign ledger + standing-issues A.6 records |
| b019ef1 | Demo Exhibit-A guard (branches on V_es_ame presence — works on both pre- and post-Group-A builds; keep on revert) |

## Revert procedure (if testing demands it)

```
git revert --no-commit 294513a c545f03 c184ea0 cae6947 ea69544 7054878
git commit -m "Revert Group A Patches (AMESE R1/R3/R4)"
```
This single commit restores: the pre-R1 ado + sthlp, the 8 pre-R1 cont battery
baselines, the pre-AMESE ledger, and removes test_amescreen.do (cert returns to
21 suites). Do NOT revert a603f38/42d12d9/439d5cb/b019ef1/f4aeb60 (records). After
reverting, verify: battery 34/34 rtol-0 (the 5 new-feature pins from c5d00d2 are
independent of Group A and must still pass), cert 21/21, and a full oldnew run
(expect the pre-Group-A 450/450 state; the AMESE ledger rows are gone with cae6947
so no stale-row fallout). The standing-issues A.6 and solved-issues entries should
then be annotated as reverted rather than deleted.

## Interactions to remember while testing

- Old-vs-new suite runs stage git HEAD, never the working tree — commit before
  expecting the suite to see a change.
- cont_uz_ss keeps its pre-existing rc pin (198->0) and is judged on rc alone.
- V_es_ame withholding is invisible to oldnew (F2 objects are new-only); it IS
  visible to the battery (the re-blessed baselines encode it).
- The T-M2 defect sites no longer execute on cohort-clustered contreat configs
  (R3 suppresses their block); the banked evidence (a603f38) predates that.

## Rider: author message revisions (2026-08-25)

The author's revised caution/error texts (all five screens) plus two display-
precedence rules (the all-zero-VCE caution and the identical-SE caution defer to the
cell-containment caution when it displayed; the refit caution became self-contained)
landed in a follow-on commit spanning both patch groups — see `git log --oneline`
for the "author message revisions" commit. Console-only plus the sthlp steer sync
(bootstrap SEs / IM p-values / randomization inference / finer clustering); no
compared surface moved (cert 13/13, battery 34/34 rtol-0). A revert of this group
should either include that commit's relevant hunks or accept the revised texts.
