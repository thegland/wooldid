# hv-fastnose-doc — FINDINGS

## HV-FN1 — REFUTED: the 2026-09-13 item-11 ruling WAS built (as (b), not (a))
**severity low (refutation) | confidence certain | seeded**
**surface:** `fastnose`, `e(b_es_*)`, `e(numpre_*)`, sthlp:1325
**repro:** `probe1_fastnose_bes.do`

The recon claim has three parts; only the first is true, and it is not a defect.

1. TRUE as text: `sthlp:1325` still opens the list with `{cmd:fastnose}`; `git log -S` shows
   `6fe814d` (2026-09-14) removed `{cmd:impvalues}` and left `{cmd:fastnose}` deliberately.
2. FALSE as behaviour: on the pinned build `fastnose` DOES post `e(b_es_att)` and
   `e(numpre_att)`. The recon's rc-111 measurement is the Stage C probe of **2026-09-13**,
   i.e. the day BEFORE the fix. Stale measurement compared against post-fix help text.
3. FALSE as a ruling audit: `git show 6fe814d --no-patch --format=%B` states
   "fastnose now posts the event-study coefficient vector, so the help's claim becomes true
   as written ... impvalues is removed from the no-identified-variance list, since the code
   skips that block there."  The author took option **(b) for fastnose + the (a) half for
   impvalues**, not the recommended (a). Code: `f2bonly` at `ado:1990`, gates at
   `ado:2162 / 2667 / 2930 / 3327`, posting at `ado:4148-4150`.

Measured, 48 units x 16 periods, onsets 6/9/12 + never-treated,
`att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod`:

| run | e(b_es_att) | e(V_es_att) | e(numpre_att) | e(wdid_vcvstate) |
|---|---|---|---|---|
| analytic   | rc 0   | rc 0   | 5 | analytic-blocks |
| + fastnose | rc 0   | rc 111 | 5 | absent |
| + timetrends | rc 0 | rc 111 | 5 | absent |
| impvalues  | rc 111 | rc 111 | . | (none) |

=> `sthlp:1325` is TRUE in both halves for `fastnose` and for `timetrends`, and `impvalues`
correctly no longer appears (matching `sthlp:276`).  **Known-by-design §7's last bullet
should be struck.**

Supporting measurements: fastnose `e(b_es_att)` equals the analytic twin to 2.22e-16 max
abs (last ulp) with identical colnames, and is BYTE-IDENTICAL to the `ES_att_*` cells of the
same run's `e(b)` (probe3, `F3_*_NDIFF=0`).

---

## HV-FN2 — fastnose prints and now RETURNS fabricated ame zeros that the analytic twin refuses
**severity high | confidence certain | NOT seeded — found while mapping the consequence surface**
**surface:** `fastnose` + `contreat()`; the R6/SI-A4 dose-slope void screen (`ado:11958`)
**repro:** `probe6_amevoid2.do` (self-contained; builds its own panel)

`ado:11958` gates the dose-slope void screen on `"\`fastnose'" == ""`. When every treated
cohort-period dose-slope term drops for multicollinearity — the help's own named case,
"the dose is constant within each cohort" — the analytic run sets `amevoid=1`, prints a
caution naming the non-identification, reports every ame row as missing and WITHHOLDS
`e(b_es_ame)`. Under `fastnose` that screen never runs, so `amevoid` stays 0 and the
fabricated zeros go all the way out.

Geometry: 30 cohorts x 6 units x 12 periods, `i` = cohort, `dcoh = exp(0.2*mod(i,5))`
(constant within cohort, range 1.00-2.23), true DGP `yV = y0 + 0.5*dcoh*D` — so the
structural dose coefficient is **0.5**, and the printed 0 is a collinearity artifact, not a
small estimate. Command: `att robust contreat(dcoh) contreatcontrols(basic) esprelength(4)
espostlength(4) esfixedbaseperiod`, +/- `fastnose`.

| | analytic | + fastnose |
|---|---|---|
| generic "interacted ... dropped for multicollinearity" caution (165 hits) | printed | printed (identical) |
| R6 caution ("ame NOT IDENTIFIED ... structural zeros, not estimates") | printed | **absent** |
| printed `main:ame`, `pretreatment:ame` | `.` | **`0`** |
| printed `ES_ame` pre5..post3 (8 rows) | all `.` | **all `0`** |
| `e(b_es_ame)` | rc 111 (withheld) | **rc 0, all 8 entries exactly +0X-3ff** |
| `e(numpre_ame)` | 4 | 4 |
| last word on the ame block | the R6 refusal | **"no horizon of the ame event study has an identified variance, so e(b_es_ame) posts without e(V_es_ame). The estimates stand"** |

The att side is identical between the arms (`main:att` .95929975 both), so nothing separates
the fiction from the sound half.

**What the user is and is not told.** They ARE told, in both arms, that some interacted
treatment-period coefficient was dropped for multicollinearity and that "this could imply
your results are not identified" — a generic caution that names no block and fires on many
benign designs. What disappears under `fastnose` is the caution that names THIS block and
says the ame entries are structural zeros rather than estimates; and what replaces it is a
caution asserting that the estimates stand.

**The crispest statement of the defect.** The R6 caution's own closing sentence reads: "You
can ignore this check if you specify safety(off), but be careful: with safety(off) the ame
output will contain the structural zeros." So the author has already ruled that zeros may
reach the output only behind an explicit opt-out and a warning. `fastnose` — which is not an
opt-out from safety, and whose documented promise is only that standard errors are skipped —
silently inherits the `safety(off)` outcome on this one screen, minus the warning.

**Why it matters.** `fastnose` is documented (`sthlp:40`) as the mode for specification
search. A user exploring a policy intensity defined at the cohort/state level (constant
within cohort by construction — an ordinary administrative shape, and `i` = onset cohort is
the help's own recommended choice) sees `ame = 0.00` at every horizon and a caution saying
the estimates stand. Running the same command without `fastnose` tells them the design does
not identify it. Since `6fe814d` the zeros also reach `e(b_es_ame)`, so a script that tests
`confirm matrix e(b_es_ame)` — the exact guard the help prescribes — now passes on a design
where the help says "a honestdid or pretrends handoff therefore finds no ame input at all".

**Tie-back to the sentence I was sent to verify.** `sthlp:1325` ends "the estimates remain
meaningful even where their uncertainty is unavailable." Under `fastnose` that is true of the
att family and FALSE of the ame family on this design — and it is the ame family's own
caution, quoted above, that asserts it.

**Documentation status.** `sthlp:903` carries the clause "and the check is skipped under
fastnose", so the *skip* is disclosed; but the `e(b_es_*)` entry states flatly that on such a
design "both e(b_es_ame*) and e(V_es_ame*) are withheld ... the structural zeros surviving
only in the raw e(b)". That sentence is now false under `fastnose`. The display half predates
`6fe814d` (`c_local amescrub 1` / `amescrublo 1` sit inside the same single `if` gate, so they
have never run under fastnose — inferred from the gate, not run against a pre-fix build);
what `6fe814d` added is the `e(b_es_ame)` matrix of zeros on top of it.

---

## HV-FN3 — `e(numpre_ame)` posts even when `e(b_es_ame)` is withheld by the void screen
**severity low | confidence certain | not seeded**
**surface:** `ado:4149`; `contreat()` + the R6 void
**repro:** `probe6_amevoid2.do`, ARM VOID MODE analytic (`Q_VOID_analytic_NUMPRE_ame=4`)

`ado:4148` guards the `b_es` posting with `!(amevoid & strpos(f2ky,"ame")==1)`; the very next
line, `cap ereturn scalar numpre_\`f2ky' = \`f2np'`, carries no such guard. So on a voided ame
design the analytic run withholds the matrix but still posts `e(numpre_ame) = 4`. The help
says such a run leaves "no ame input at all" for honestdid/pretrends. Harmless on its own
(honestdid still has no `b`), but it is the wrong half of the pair to keep, and it is the
scalar a user's script is most likely to branch on.

---

## Clean (probed hard, nothing found)
- `fastnose` + `subgroups()`: BS-26 keeps the subgroup-horizon `margins` VCE for the
  estimability markers and then blanks it. ZERO non-missing `se/t/z/p/lb_95ci/ub_95ci` cells
  across `e(wdidmainresults)`, both `e(wdidsgresults_*)` and `e(pwsgs)`; `e(V)` absent;
  `e(b_es_att_s0/s1)` post, `e(V_es_att_s0/s1)` absent. (probe2, `E2_LEAK_*=0`.)
- `permute(shuffle)` fast mode, with and without an explicit user `fastnose`: posts no
  `e(b_es_att)` / `e(numpre_att)` — the `permute==""` carve-out in `f2bonly` holds in both
  directions. (probe3 F2/F5.)
- honestdid after a `fastnose` run fails informatively with honestdid's own
  "Last estimation vcov matrix not found; please specify matrix." r(198). (probe3 F4.)
- Cheap audit of the other twelve 2026-09-13 rulings: item 1 `DF-1 (author ruling 2026-09-13)`
  at `ado:1258`; item 2 commit `1a52a8e`; item 4 pointer at `sthlp:~63` + the 0.66-0.75 range
  at `sthlp:255`; item 12 `e(es_ndegen_*)` named at the honestdid handoff `sthlp:1321` and the
  example carries `parallel(0)` at `sthlp:1567`; item 13 `esprelength(0)` with its reason at
  `sthlp:1494`. All present. Items 3, 5, 6, 7, 8, 9, 10 are suite/report-policy items not
  chased here.

---

## Addendum to HV-FN2 — how far the R6 skip reaches (probe7_permvoid.do)
`permute()` without `permtstat`/`permanalyticse` sets `local fastnose fastnose` internally
(`ado:5448`), so the OBSERVED permutation run also skips the R6 screen. On the same void
geometry it prints `main:ame`, `pretreatment:ame` and all four `ES_ame` horizons as
`estimate 0, p_perm 1.000, 95% CI [0, 0]`, with the R6 caution absent.

**But permute rescues itself**: it prints
"Caution: the permutation distribution is degenerate (every completed draw equal) for:
main_ame pretreatment_ame ES_ame_pre3 ES_ame_pre2 ES_ame_post0 ES_ame_post1; those
permutation p-values are uninformative." — every affected row is named, and a zero-width
interval at p = 1.000 is self-evidently degenerate. `permanalyticse` (which leaves
`fastnose` unset) restores the R6 caution and the missing ame rows.

So the ranking of exposure is: **plain `fastnose` is the worst case** — zeros printed as
point estimates, no degeneracy signal at all, and a caution that says "The estimates stand".
`permute()` fast mode prints the same zeros but names them. `bootstrap()` is CLEAN — measured, not inferred (probe8_bootvoid.do).

The display half of HV-FN2 predates `6fe814d`: `c_local amescrub 1` / `amescrublo 1`, which
blank the displayed ame cells, sit inside the same single `if "\`contreat'" != "" &
"\`fastnose'" == "" & strpos("\`bvars'","o.") != 0` gate, so they have never run under
fastnose. What `6fe814d` added is the `e(b_es_ame)` matrix of zeros on top of it.

### probe8_bootvoid.do — bootstrap() on the void geometry is CLEAN (measured)
`bootstrap(standard) reps(20) rseed(3)` on the same VOID geometry: the observed run is not a
fastnose run, so the R6 caution fires, every ame row prints `.`, `e(b_es_ame)` is withheld
(rc 111), and the engine adds its own disclosure — "Caution: the following estimates have no
identified observed value and are excluded from bootstrap inference: main_ame
pretreatment_ame ES_ame_pre3 ES_ame_pre2 ES_ame_post0 ES_ame_post1" — even though every
replication (and the rep-0 comparator) runs under fastnose and therefore returns ame = 0.
The replications' zeros never reach the table. No defect on this path.
