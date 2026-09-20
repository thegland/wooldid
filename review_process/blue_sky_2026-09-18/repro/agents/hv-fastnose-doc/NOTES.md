# hv-fastnose-doc — verification notes

## The claim under test
recon: "2026-09-13 ruling item 11 appears unbuilt. sthlp:1325 still opens the
'covariance not returned but coefficient vector still is' list with {cmd:fastnose},
while under fastnose e(b_es_att) is NOT posted (rc 111) and e(numpre_att) is empty."

## VERDICT: REFUTED. The ruling WAS built — as option (b)+(doc trim), not (a).

### Documentary evidence
- `git log -S` over `wooldidinstall/wooldid.sthlp`: commit `6fe814d` (2026-09-14) removed
  `{cmd:impvalues}` from the parenthetical and deliberately LEFT `{cmd:fastnose}`.
- `git show 6fe814d --no-patch --format=%B` line 42: "fastnose now posts the event-study
  coefficient vector, so the help's claim becomes true as written and downstream handoffs
  have something to read; the covariance stays absent, which is the other half of that
  claim. ... impvalues is removed from the no-identified-variance list, since the code
  skips that block there."  => the author chose (b) for fastnose + the (a) half for impvalues.
- ado:1982-1990 carries the matching comment and the new local `f2bonly`
  (`= ("fastnose"!="" & "permute"=="")`), and ado:2162/2667/2930/3327 gate the F2 b-capture
  on `"marginsdof"!="nose" | f2bonly`.
- The recon's rc-111 measurement is the Stage C probe of 2026-09-13 — taken the day BEFORE
  the fix. Stale measurement vs post-fix help text.

### Measured on the pinned build (probe1_fastnose_bes.do / .log)
48 units x 16 periods, onsets 6/9/12 + never-treated, `att cluster(i) esprelength(5)
espostlength(6) esfixedbaseperiod`:

| run | e(b_es_att) | e(V_es_att) | e(numpre_att) | e(wdid_vcvstate) |
|---|---|---|---|---|
| A analytic          | rc 0   | rc 0   | 5 | analytic-blocks |
| B + fastnose        | rc 0   | rc 111 | 5 | absent |
| C + timetrends      | rc 0   | rc 111 | 5 | absent |
| D impvalues         | rc 111 | rc 111 | . | (none) |

B prints the documented caution verbatim: "Caution: no horizon of the att event study has
an identified variance, so e(b_es_att) posts without e(V_es_att)...". e(V) itself is absent
under fastnose (rc 111) and every se/t/p/ci cell in the table is missing.
So sthlp:1325 is TRUE in BOTH halves for fastnose, and TRUE for timetrends;
impvalues correctly no longer appears in the list and matches sthlp:276.

## Consequence surface (probe2_surface.do, probe3_surface2.do)
- E1: fastnose b_es == analytic b_es to 2.22e-16 max abs (last-ulp only); colnames identical
  (`pre6..pre2 post0..post5`, pre1 omitted — the honestdid signature the help describes).
- F3: e(b_es_att) cells are BYTE-IDENTICAL to the `ES_att_*` cells of e(b) in BOTH arms
  (0 differing cells), confirming the commit's "byte-identical to cells already present".
- E2: fastnose + subgroups(): e(V) absent; e(b_es_att_s0/s1) post; e(V_es_att_s0/s1) absent;
  ZERO non-missing se/t/z/p/lb/ub cells across e(wdidmainresults), both
  e(wdidsgresults_*) and e(pwsgs). BS-26's "keep the sg margins VCE for the estimability
  markers then blank it" leaks nothing. CLEAN.
- F2/F5: permute(shuffle) fast mode — with OR without a user-supplied `fastnose` —
  posts NO e(b_es_att)/e(numpre_att). Matches sthlp:283 and the `permute==""` carve-out.
- F4: honestdid after a fastnose run fails with honestdid's own
  "Last estimation vcov matrix not found; please specify matrix." r(198). Informative.
- F1: contreat + contreatcontrols(byCohort) on a dose that is 0 in every untreated period
  is refused r(459) with the documented cohort-naming message. Re-run with
  contreatcontrols(untreatedZero) in probe4.

## Harness notes
- DGP trap that cost one run: `D*(0.8+0.15*rel)` is MISSING for never-treated because
  `rel` is missing and `0*.` = `.` in Stata. wooldid then correctly reported
  "no never-treated control group" r(459). Use cond(D==1, ..., 0).
- After a failed `cap noi wooldid`, e() still holds the PREVIOUS run's results.
  `ereturn clear` before every capture-wrapped run or you will read stale matrices.

## probe5/6/7 — the NEW finding (HV-FN2)
probe5 (i = group, contreatcontrols(untreatedZero), dose constant within cohort) died r(322)
in BOTH arms — the cohort-period dummy drops before the R6 screen is reached. The geometry
that REACHES R6 is the demo_vesame_pathology shape: i = cohort, 6 units per cohort,
contreatcontrols(basic), dose a function of i only. probe6 uses it and reproduces cleanly.

probe6 VOID arm: analytic withholds e(b_es_ame) and prints every ame cell as `.` under the
R6 "NOT IDENTIFIED" caution; fastnose prints every ame cell as `0`, posts e(b_es_ame) with
8 exact zeros, and says "The estimates stand". att side identical in both arms.

probe7: permute fast mode also skips R6 (ado:5448 sets fastnose internally) but names every
affected row in its own degeneracy caution and shows p=1.000 with a [0,0] interval.
