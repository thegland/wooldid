# Held for the author, 2026-09-13 — ALL RULED 2026-09-13, BUILT 2026-09-14

**Every item in this file was ruled.** The rulings were taken the same day in multiple-choice
batches and carried out the next: `6fe814d` ("Implement the author's 2026-09-13 rulings: degrees
of freedom, the weighted Bayesian bootstrap, and six smaller calls") opens "Fourteen rulings were
taken today in multiple-choice batches. This is the code half", and its siblings on 2026-09-14
carry the rest - `1a52a8e` for the subgroup pre-treatment drop, `81b3ce2`/`cb1796b`/`b0a783f` for
the ledger rows the degrees-of-freedom fix falsified, `64c03dd` for the two battery baselines it
moved, `af3c47c` for the suite rows the subgroup fix needed.

The file is kept as the record of what was put to the author and on what evidence. **The body
below is unchanged from the day it was written**, so its own line "nothing in it has been decided"
is true only of that morning. This header is the correction.

Where the rulings are recorded, since they are not in this file: in the code and data themselves.
`wooldid.ado:1258` carries `DF-1 (author ruling 2026-09-13)` with the full measurement for item 1;
`degeneracy_expected.csv` writes the keep-as-pins decision into its rows as
`PINNED_DEGENERATE (ruled 2026-09-13: keep as a pin ...)`, which `DEGENERACY.md` and the suite's
own degeneracy warning both echo at a total of twelve pinned rows, the nine of item 8A among them; item 13's `esprelength(0)` condition is at `sthlp:1495`.

**Verification status, stated honestly.** Items 1, 8 and 13 were confirmed directly in the code
and data files on 2026-09-16. Items 4 and 11 rest on `6fe814d`'s message plus partial code
evidence. The remainder rest on that commit message alone. A systematic pass over all thirteen
implementations has NOT been run; the project has previously found ruled items quietly
unimplemented after being declared done, so "the commit says so" is not the same as verified.

**Why this header exists.** The file was committed once and never updated to record its own
outcome, and on 2026-09-16 that stale opening line was read as evidence that thirteen decisions
were still waiting on the author. Searching later commits for the internal item codes (`df_a`,
`F058`, `T-FN1`) found nothing, because the rulings were recorded in prose - "degrees of freedom",
"fastnose" - and never by code. A status line that is only true on the day it is written should
say so.

---


Thirteen items reached during the 2026-09-12/13 assessment chunk (Stage A suite audit, the Stage B
old-vs-new proof drafts, and the Stage C benchmark feasibility probes) where the work cannot go on
without a decision: a number would move, a contract would change, or two records disagree about what
a shipped behaviour is. Each section carries the analysis and the evidence exactly as it is being put
to the author. The **ruling wanted** column is the multiple choice; nothing in it has been decided.

Each item is marked with which kind of thing it is, using the project's own three-way split:
**(i) straight bug** — doctrine says fix it without asking, listed here only because it travels with a
ruling; **(ii) number-moving** — held for the author by standing instruction; **(iii) documentation
wording** — no number moves. Items 11-13 are mostly (iii) and are listed so the author sees the
sentence that is being changed under his name.

| # | item | kind | ruling wanted (pick one) | recommendation |
|---|---|---|---|---|
| 1 | SI-20 `df_a` over-count from the absorbed ringer terms | (ii) | (a) accept the extra absorbed dof as built; (b) skip the null/lattice ringer terms when a `byX` term is already present, restoring the pristine `df_r`; (c) accept but pin the `df_r` of every lattice+basic config in the suite so the factor can never drift silently | (a), plus the (c) pins |
| 2 | SI-7 blank rows vs dropped rows: the bitwise-zero subgroup pre-treatment estimate under `safety(off)` | (ii) | (a) keep as built — report it with the `fcguard` caution; (b) drop it, as the `fullprecisionslow` path and BS-8 do, so all three routes agree on dropping; (c) blank it (row present, every cell missing) on all three routes | (b) — one rule on all three routes |
| 3 | `mc_ri` M1: the preregistration was written for `permabsp`, the shipped default is equal-tailed | (ii) | (a) amend the preregistration to the equal-tailed default (bands (8,24), R6 judged at .10 only); (b) pass `permabsp` in `mc_ri.do` so the driver measures the preregistered statistic; (c) arm both — keep the default arm at the new bands and add a `permabsp` arm; (d) change the shipped default back to `permabsp` | (a), and (c) if a second arm is affordable |
| 4 | Unit-level `i` guidance in the help contradicts the help's own clustering warning | (iii) | (a) add a pointer at `sthlp:67` to the "What level to cluster at" paragraph and scope the efficiency claim to the point estimates; (b) additionally give the measured att-row understatement range at `sthlp:1303`; (c) leave both passages as they are | (b) — it is (a) plus the number item 6 supplies |
| 5 | `lattice(all 3)` returns an all-zero VCV on realistic 5-12-unit cells, under `robust` as well | undetermined | (a) accept as the documented all-zero-VCE outcome (`sthlp:1076`) and stop there; (b) spend one diagnosis run to establish whether the 254 collinearity drops are the documented containment or a lattice-construction defect; (c) add a lattice-specific caution without diagnosing | (b) — the class cannot be assigned without it |
| 6 | The F058 82.8% cell: both standing attributions are refuted | (iii) + (ii) | (a) document it by geometry — correct the F058 record and give the measured att-row range at `sthlp:1303`; (b) run the shock-removal MC first (~4 min) to isolate the mechanism, then document; (c) correct the record only and leave the help alone | (b) then (a) |
| 7 | `test_level` L11/L12 copied-number pins | test policy | (a) keep them — a copied number is admissible when an identity sits beside it; (b) strip the four-number checks and keep the identities; (c) replace the numbers with independent-route derivations | (a) |
| 8 | The ten UNRULED degeneracy expectations (nine HIDDEN configs + `ext_b4_h55_p`) | suite curation | **8A, the nine HIDDEN rows:** (a) keep every one as a pin of its degenerate surface and let the new twins carry the live surface; (b) re-host each on a host where the surface is live and retire the degenerate config; (c) mixed: keep the three with landed twins, re-host or drop the rest. **8B, `ext_b4_h55_p` only** (this half is a code change, not suite curation): (a) accept the silent pre5 loss and close the clause as BY_CONSTRUCTION; (b) let the existing forced-reference caution through even under `safety(off)`, so the loss is announced; (c) shorten the grid corner's window so the config stops sitting on the data's edge | 8A (a); 8B (b) |
| 9 | H-9: bounded ledger rows firing at 0.47-0.58 of their bound | (ii) | (a) accept the tighter-than-policy bounds knowingly and record the ceilings; (b) widen the three 1e-6 rows to 5e-6 and the two 2e-4 rows to 1e-3, ceilings recorded; (c) widen only the three 1e-6 rows | (a) — the 18 ceiling appends ship regardless |
| 10 | NFADV-01/02/05 classification for the final report | (iii) | (a) review-verified defect for all three; (b) author-approved recalibration for all three; (c) split — NFADV-05 a defect (fabricated anchor), NFADV-01/02 a recalibration | (a) |
| 11 | T-FN1: `sthlp:1295` says `fastnose` posts `e(b_es_*)`; it does not | (iii) + (ii) | (a) correct the sentence — remove `fastnose` (and `impvalues`) from the list, leaving the all-zero/`timetrends` case that does post; (b) change the `fastnose` path to post `e(b_es_att)`/`e(numpre_att)`; (c) both | (a) |
| 12 | T-HD2: honestdid does not reliably refuse a degenerate covariance | (iii) + (ii) | (a) help sentence at the honestdid handoff naming `e(es_ndegen_*)` as the guard, plus a cert leg; (b) (a) plus withhold `e(numpre_*)` when `e(es_ndegen_*)` is nonzero, so the handoff has no inputs; (c) (a) plus a run-time caution at the handoff | (a) |
| 13 | T-BJS1: the `did_imputation` equivalence conditions omit `esprelength(0)` | (iii) | no ruling needed — straight documentation fix; ruled only if the author wants the sentence phrased differently | fix at `sthlp:1464` |

The per-item sections below are the analysis as it is being put to the author.

---

## 1 — SI-20 `df_a` over-count from the absorbed ringer terms

**Kind: (ii) number-moving.** Nothing here is a bug; the question is whether the degrees of freedom
the fixed model consumes are the right ones.

**What the ship record already said.** `review_process/si20_ship_2026-09-05.md` §4 item 1 (line 108)
lists this as the first of five "judgment calls for the author (nothing below was changed)":

> reghdfe counts the added slopes in `df_a`, so V/SE/t/p/CI on EVERY row (att rows included) scale by
> `df_r_old/df_r_new`: +1 df on the plain-basic default path (0.035% at N = 1728; V x1.115 on
> dta_tiny with ~9 residual df), +4-6 df under basic byX + lattice where the terms are redundant with
> the already-absorbed byX#lattice slopes (V x1.0028..1.0050; joint p in the 3rd digit).

**No ruling was ever made.** A recursive grep of `review_process/` (the other `df_a` hits are
unrelated), `qr_ship_2026-09-06.md`, `held_for_ruling_2026-09-09.md`, `rulings_2026-09-07.md` and the
rulings memory finds nothing (scratchpad `stageB_proofs/si20/NOTES.md` lines 48, 51). The ledger rows
document it as built.

**What the Stage B reproduction measured** (scratchpad `stageB_proofs/si20/SECTION.md` §S.1, §S.3,
§S.4; `stageB_proofs/contsg/SECTION.md` §S2, block `R4W`; both builds, separate Stata processes, dumps
bit-identical to the suite's 2026-09-11 dumps):

| config | `df_r` old -> new | the scale factor, measured |
|---|---|---|
| `si20_lat_pw` | 1268 -> 1267 | `main_att` point estimate 1.3808962103 unchanged; se ratio new/old = 1.000394555 = `sqrt(1268/1267)` to nine digits |
| `si20_lat_aw_bt` | 1414 -> 1413 | se ratio 1.000353795 = `sqrt(1414/1413)` |
| ticket fixture (`lattice(baseTreated 2)`, `robust`) | 1414 -> 1413 | se ratio 1.000353794 |
| `si20_lat_byt` (`basic byTreated` + lattice) | 1263 -> 1258 | `SI20-P-si20_lat_byt-V` bounded at 1e-2, measured 3.96e-03; `b` pinned at 7.07e-08 |
| `cont_lat` | absorbed dof 754 -> 753, recovered exactly from `rss / rmse^2` | `rmse` .0098460017 -> .0098525374 (**6.6e-4**), `r2adj` 2.2e-9, `r2_within` 2.3e-6, `ll_0` -2419.22 -> -2315.56; **no** estimate, standard error or test moves on this config |

`cont_lat` is the clean read: it is a fit-statistic-only config (`e(properties)` = `b` on both builds),
so the extra absorbed dof shows up there with nothing else riding on it. `rss / rmse^2` = 754.000000
old, 753.000000 new.

**One thing in the ship record's own option list is now stale.** Option (b) as written — "skip the
null/lattice terms when a `byX` term is present" — is annotated "restores pristine SEs on reghdfe but
re-exposes the R4 reg-refit failure, where reg dropped the byTreated ringer and re-introduced the
anchor". The R4 refit was **removed** by `b801282` on 2026-09-08 (ruling `rulings_2026-09-07.md` §J;
record `review_process/research_2026-09-08/fastnose_precision.md`). There is no refit left to
re-expose, so option (b) is cheaper today than the record says. It remains a solve-path change on a
model the author has already shipped.

**What rides on the answer.** If the author rules that the redundant absorbed slopes should *not*
count, then (a) the att-side se/t/p/ci movement recorded under every `SI20-L` ledger row, and the three
`SI20-P-*-cfxatt` rows, become a residual defect of the fixed build rather than a mechanical side
effect (`stageB_proofs/si20/SECTION.md` §S.4); (b) `cont_lat`'s new `rmse` is 6.6e-4 too large on
lattice+basic runs (a fit statistic only); and (c) the two wide bounds the harness audit
flagged — `SI20-P-si20_lat_byt-V` at 1e-2 and the three `SI20-P-*-cfxatt` rows at 1e-1 — exist only to
absorb this factor (harness_health audit §7.4, "Wide bounds"), so item 9 partly resolves with this one.

**Recommendation: (a), accept.** The non-lattice `basic` path already carried the +1 df before the fix;
the factor is a genuine consequence of moving a real control into the absorb set, not an accounting
error; and the largest measured effect on a reported number is `V x1.0050`. Add the `df_r` pins of
option (c) so the factor cannot drift unremarked.

---

## 2 — SI-7 blank rows vs dropped rows: the bitwise-zero subgroup pre-treatment estimate under `safety(off)`

**Kind: (ii) number-moving.** Not a wrong number — an `e(b)` whose width and contents change with the
ruling.

**The question, in one line.** When a subgroup's pre-treatment horizon is not estimable, wooldid today
does three different things on three paths. Which one is the rule?

| path | what happens on the non-estimable subgroup horizon |
|---|---|
| SI-7 drop (`fullprecisionslow`, and the default before route P) | the horizon is **dropped** from that subgroup's event study with a caution naming it |
| BS-8 (continuous subgroup horizons, `49bfb52`) | the twelve structural zeros are **dropped**: `cont_sg_gap` `e(b)` 45 cols -> 33 |
| the default path since route P (`4b224f0`, 2026-09-02; collapse default since `111ecbb`) | the value is **reported as a bitwise 0** under the `fcguard` caution |

**The measurement** (`stageB_proofs/contsg/SECTION.md` §S4, block `SI7` on
`ext_contract10_rcs_factor_fe`: `wooldid y_rcs state_group t ttre [pw=survey_weight], att
subgroups(group) fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state)
safety(off)` on `test/suite/data/survey_rcs.dta`, both builds, plus a new-only `fullprecisionslow` leg):

```
SI7_b_pretreatment_SG1_att  OLD ABSENT (ncols 11)     NEW 0 (firstpos 6 of 12)
SI7_sg1_shape               OLD 3x7 [att pre2 post0]  NEW 4x7 [att att pre2 post0]
SI7_b_ES_att_SG1_pre2       OLD 8.8817842e-16         NEW 8.8817842e-16   (retained on both)
SI7_wdid_fcguard            OLD ABSENT                NEW [caution]
SI7FPS_sg1_shape  NEW 3x7 [att att post0];  SI7FPS_b_ES_att_SG1_pre2 ABSENT;  SI7FPS_wdid_fcguard ABSENT
```

The `fullprecisionslow` leg shows exactly the behaviour the four `SI7-*` ledger reasons describe; the
default path does not. Why: under route P the cell collapse runs with `subgroups()` behind a per-read
degeneracy screen, and under `safety(off)` "there is no fallback — the run proceeds on the collapsed
rows and prints a caution"; on the collapsed rows the SG1 pre2 margin is **not omission-marked**, so
the machine zero survives as an estimate.

**The author was asked once and there is no ruling.** `review_process/qr_ship_2026-09-06.md` line 107:
"rows versus dropped rows is a display choice; SI-7 consistency argued for dropping." Nothing in
`rulings_2026-09-07.md` answers it.

**Cost of each option.** (b) *drop everywhere* is the consistent rule and matches the two paths that
already drop, but it needs a detection mechanism the collapsed path does not have — the margins
omission marker the SI-7 reader keys on is absent there, so the screen would have to test the value or
the cell count instead. **evidence: thin —** the code site was not instrumented in this pass; §S4
records only that the pre-treatment-row match at ado 3500-3514 is byte-identical to old's 2815-2829,
so the flip is a numerical knife edge on the margins colname, and the exact colname old saw when it
lost the SG1 row was never bisected. (c) *blank everywhere* keeps `e(b)`'s width stable across
`fullprecisionslow` and would require re-pinning `SI7-eb-*` and `SIA4`/`BS-8` rows in the other
direction. (a) *keep as built* leaves three rules in one package and the four `SI7-*` ledger reasons
permanently describing a path the default run does not take.

**Whichever way this goes**, the four `SI7-*` reason texts are rewritten (they describe the horizon
drop, 11 -> 10 names; the run shows the `pretreatment_SG1_att` addition, 11 -> 12) or the config is
pinned with `fullprecisionslow` so the rows describe what fires.

**Recommendation: (b).** One rule on all three routes, and it is the one two of the three already
follow.

---

## 3 — `mc_ri` M1: the preregistration measures a convention the package no longer ships by default

**Kind: (ii) number-moving.** This is a blocker for Stage C: R6 cannot pass at any seed.

**The facts.** HEAD's default permutation p is the equal-tailed add-one rule (ado L7601-7606;
`sthlp:426-431`). `permabsp` is HEAD's name for the older `|s| >= |obs|` rule (ado L15374-15378).
`test/suite/PREREG_NEWFEAT_2026-08-24.md:35-37` (§M1, "exactness being tested") pins the old formula:
*"tail count is inclusive (|s| >= |obs|); sampled p = (1+nex)/(nuse+1); exact p = nex/nuse"*.
The convention changed after the 2026-08-24 run (`1561386`/`dc3f345`; `permexact` BS-23 `7e30e85`,
`4484493`).

**Why R6 cannot pass.** R6 is `permexact` on a complete 20-assignment enumeration, and the observed row
is in the exact reference set (`wooldid_permmemb`, ado L15412-15417), so the equal-tailed p is
`2 min(R, 21-R)/20 >= .10` on every draw. Rejections at .05 are therefore identically 0 against the
preregistered band [11, 29]: **R6 FAILs on every run, at every seed**. Re-measured at seed 11010001:
p .60 equal-tailed vs .55 under `permabsp` (critique probe `probe1.log:17228`, `:17268`).

**And the other six configs are mis-banded.** At B = 99 the equal-tailed add-one attains
`P(p <= .05) = 4/100 = .04` exactly (R in {1, 2, 99, 100}), not .05. Under Binomial(400, .04) the
preregistered band (11, 29) is missed with probability .0742 per config (.0733 below, .0009 above), so
across R1-R5 and R7 the chance of at least one spurious size FAIL per production run is
`1 - (1 - .0742)^6 = .37`, on top of R6's certain FAIL. Route (a)'s replacement band (8, 24) gives
`P(<8) = .0089`, `P(>24) = .0200`.

**One thing the audit got wrong, corrected by the critique.** The audit called the `permabsp` edit
"inert under `permtstat permtci`" on R3. It is not: the audit's patched copy never changed R3's line.
Measured at seed 11010001, R3 with `permabsp` moves `e(perm_pconv)` from `equaltailed` to `absolute`
and p from .68 to .63, with `lb`/`ub`/`tlb`/`tub`/comparator bit-identical (`%21x`). So route (b)'s R3
edit is load-bearing, and R3's `TCI_COV` secondary is convention-invariant either way.

**Option (d) — change the shipped default back — would overturn a standing ruling** (the equal-tailed
default was built at `1561386`/`dc3f345` and the help documents it). It is listed because it is the
only option that makes the preregistration true as written without editing it.

**Note on the lettering:** the options above are lettered independently of `mc_drivers/REPORT.md` F1,
which calls this document's (a) its *route (b)* and this document's (b) its *route (a)*. The audit's
own recommendation ("(b) plus keeping R3 as is") is the same substantive choice as (a) here.

**Recommendation: (a).** The suite should certify what users actually get. Amend the preregistration
with a documented deviation note (its own L25 rule requires one), set `RI_BAND_05 = (8, 24)`, and judge
R6 at .10 only. Add (c)'s `permabsp` arm if a second arm is affordable, since after (a) the `permabsp`
path has no size guard at all.

---

## 4 — Unit-level `i` guidance in the help contradicts the help's own clustering warning

**Kind: (iii) documentation wording.** No number moves. The audit's "D-1 SE understatement" was filed
as a blocker and **refuted by the critique** as a rediscovery of the documented saturated-cell CRVE at
unit-level `i`; what survives is that the help steers users into the configuration it separately warns
about.

**The two passages, verbatim.**

`wooldidinstall/wooldid.sthlp:67` (the tail of the `i` entry, lines 65-67):

> Simulations suggest the estimator is most statistically efficient when {it:i} is set at the level at
> which you cluster.

`wooldidinstall/wooldid.sthlp:251` (the "What level to cluster at" paragraph opening the Standard
Errors section):

> The intended practice is to cluster {ul:below} the level of {it:i} - at a level whose units are
> nested inside cohorts, so that a treated cohort-by-period cell contains several clusters. At or above
> the level of {it:i} the cluster-robust covariance stops seeing the treated cells' own variation, and
> two things follow that the table does not announce: the pointwise standard errors on the {it:att}
> rows are understated (measured at roughly 0.66 to 0.75 of the reference on the diagnostic fixtures),
> and the joint pre-trend tests are oversized [...]

Two further passages carry the same warning from the other side: `sthlp:1075` ("binary-treatment (att)
standard errors under cohort-level clustering are also somewhat understated on such saturated designs
(a milder, documented phenomenon); wooldid does not suppress them") and `sthlp:1303`, condition **Two**
of the `e(V_es_*)` design conditions ("with a unit-level `i`, each treated cohort-by-period cell is a
single observation, the treated residuals vanish, and the covariance carries control-cluster noise
only [...] the event-study standard errors sit at 0.84 to 1.07 of the pooled event study's [...] the
numbers look entirely ordinary and the object is not measuring what it appears to").

**These are not literally contradictory** — `:67` is an efficiency claim about the point estimates,
`:251` is a validity claim about the covariance. Both can be true at once. The defect is that `:67`
tells the user to set `i` at the cluster level and says nothing about what `:251`, `:1075` and `:1303`
then say about his standard errors, and neither `:67` nor `:251` mentions `robust` or cluster-above-i.

**What was measured** (critique_datasets D-1, reproduced to the last printed digit; `crit1.log:17080`,
`:17134`, `:17209`, `:17296`): on `test/data/dta_binary.dta`, `att cluster(i) esprelength(2)
espostlength(2)` with `i` = the unit id gives post0 se .073601, against .103944 for the direct
cohort-cell `lincom` on the same points. On the P1 fixture at seed 1001, .062324 vs .097415. The
control that isolates the mechanism: pass a **cohort-level** id and cluster on the unit, and the ratio
collapses to a uniform 1.017 across att, pre1 and post0 (.105715 / .097868 / .071001 against
.103944 / .096228 / .069812) — an ordinary small-sample factor from the different absorbed set, points
identical. So there is nothing wrong between reghdfe's `e(V)` and wooldid's `e(V_es_att)`; the axis is
the level of `i`. v1.3 gives the same .0623 to 7 digits, so this is pre-existing and invisible to
old-vs-new parity. The 40-seed MC behind "sd/se 1.14-1.50, coverage ~86%" was not re-run by the critic
and is taken as reported.

**What the fix is.** A qualifier and a cross-reference at `:67`, not a deletion. Option (b) adds the
number item 6 supplies, which is the only att-row magnitude the help would then carry for the
singleton-cell geometry.

**Recommendation: (b).** This is the cheapest change in the batch and it closes the one place where the
help's guidance and the help's warning pull in opposite directions.

---

## 5 — `lattice(all 3)` returns an all-zero VCV on realistic multi-unit cells, under `robust` too

**Kind: undetermined.** It is either the documented all-zero-VCE outcome or a lattice-construction
defect, and the evidence does not separate them.

**The measurement** (datasets audit item 12; probes C1-C3 on `p1_probe.dta`, the realistic
12-state / 5-12-counties-per-cell fixture, unit-level lognormal dose with 10% time noise;
`probe_p1p2` run5, 17 s):

| run | result |
|---|---|
| `y_dose st year ttre, att contreat(dose_it) contreatcontrols(basic) lattice(all 3) cluster(unit) esprelength(2) espostlength(2)` | "all-zero variance-covariance matrix" caution, **12 of 12 SEs missing**; ame -2.16, pretreatment ame -10.3 (wild points) |
| the same with `robust` (C2) | all-zero again |
| the identical spec **without** lattice, `cluster(unit)` (C1) | identified: ame .394 (se .39) against truth .40 |
| the identical spec **without** lattice, `robust` (C3) | identified: ame .394 (se .145) |

**Why the obvious explanation does not fit.** The mechanism the SI-20 ruling gave for
`ext_san_input_missing` — clusters fewer than parameters — cannot be what zeroes V under `robust`.

**What the critique added** (critique_datasets D-6, verbose leg): the run is rc 0, N 2061, `df_r` 1442,
and the fitted model is the full lower-order `##` expansion
`reghdfe y_dose i.cell#i.lat##i.lat##c.dose_it ..., absorb(i.st i.year i.lat##c.(base) ...) vce(robust)`
over st x year x bin cells of 1-4 counties, with **"the multicollinearity caution above was triggered
254 times"** before the all-zero caution. `sthlp:1076` already documents the all-zero state as a known
outcome ("On some of these designs the estimation command returns an ALL-ZERO covariance outright [...]
wooldid says so and continues"). `all` is a valid `latticetype` (ado L999 comment, L1024-1044), so the
syntax is not at fault.

**The open question is exactly one thing:** whether those 254 drops are the containment `sthlp:1076`
describes, or the lattice expansion generating columns it should not. Reproduction is two 17-second
commands on `p1_probe.dta`.

**Effort-cap doctrine applies but does not settle it.** This is not a hopeless dataset — 12 states with
5-12 counties per cell is the realistic case the P1 fixture was built to be — so the "don't spend heavy
effort on pathological data" cap does not cover it. The failure is observable (a caution fires and
every SE is missing), which is the doctrine's other test, and that argues for (a).

**Recommendation: (b).** One diagnosis run assigns the class. Until then the item cannot honestly be
called either a bug or a documented property, and the suite has routed `real_cont` around it
(`configs.csv`: `real_cont` runs `contreat` on P1 with **no** lattice, citing this finding).

---

## 6 — The F058 82.8% cell: both standing attributions are refuted

**Kind: (iii) record correction, plus (ii) if the help gains a number.**

**The cell.** `test/suite/mc/out2/bootstrap/tables.md:31`, `:46`, `:63`, `:111` and
`review_process/research_2026-09-06/F058_bootstrap.md:206`, `:221`, `:238`, `:286`, `:316`: design B,
plain analytic `att`, mean SE / SD 0.702 (0.104), coverage@95 **82.8% (24 of 29)**, -3.02 MC standard
errors. It is the only cell in the F058 study past two MC standard errors on its own n, and it is the
analytic path, not a bootstrap one.

**Two attributions were on record, and the discriminating MC refutes both.**
- The **F058 record itself** (`F058_bootstrap.md:318`) attributes it to the "unequal-cluster design";
  the design is described at `:322` ("six treated clusters of sizes 1, 3, 5, 7, 9, 11 with the effect
  proportional to cluster size").
- The **Stage A audit** (KF-4) attributed it to SI-A1 condition One and proposed retiring the cell as a
  ruled, disclosed property, mapping 0.702 onto the help's 0.66-0.75 band at `sthlp:251`.

The critique ran the audit's own proposed discriminator plus three controls (`kf4/kf4_disc.do`, plain
analytic `att cluster(state) safety(off)` exactly as the BA/plain arm, 100 reps per design, 500 fits
all rc 0, 530 s in one Stata session; design-B draws reuse the F058 seeds 7600001+m and draw 1
reproduces the F058 BA/plain row bit for bit: b 1.47511213895, se .22466688461):

| design | n | truth | mean b | SD(b) | mean SE | SE/SD (MCSE) | coverage@95 (MCSE) |
|---|---|---|---|---|---|---|---|
| B0: G=12, unequal sizes 1..12 (= F058 design B) | 100 | 0.79444 | 0.7744 | 0.2840 | 0.2019 | 0.711 (0.051) | 85.0 (3.6) |
| B1: G=12, **equal** sizes | 100 | 0.60000 | 0.5987 | 0.2698 | 0.2010 | 0.745 (0.053) | 87.0 (3.4) |
| B2: G=48, unequal | 100 | 0.79444 | 0.7877 | 0.1299 | 0.1055 | 0.813 (0.058) | 91.0 (2.9) |
| B3: G=48, equal | 100 | 0.60000 | 0.5957 | 0.1336 | 0.1069 | 0.800 (0.057) | 89.0 (3.1) |
| A0: F058 design A — the **same** i = unit / cluster-above-i / singleton geometry | 100 | 2.00000 | 2.0116 | 0.3585 | 0.3585 | 1.000 (0.071) | 92.0 (2.7) |

Reading: (1) cluster-size **inequality is not the cause** — equalising the sizes moves nothing (B1 vs
B0 and B3 vs B2 agree within one MCSE), so the F058 record's attribution is wrong. (2) **Few clusters
are part of it** — quadrupling G lifts the ratio from 0.71 to 0.81 and coverage from 85 to 91. (3) A
residual understatement of about **0.80** (-3.2 MCSE from 1) persists at G = 48 on both size patterns,
so it is structural to design B — while design A, which shares the geometry the audit blamed, is
calibrated at 1.000. So the audit's retirement basis is refuted by its own control.

**What is left is missing coverage, not a retired item.** The candidate mechanism is the help's
singleton-cell condition **Two** (`sthlp:1303`): with `i` = unit the treated unit-by-relative-period
dummies saturate every treated cell (wooldid absorbs `i.i` and `i.t`, ado:1395), the treated states'
post-period state-by-period shocks (sd 0.5 in this DGP) are absorbed into those dummies, and the
cluster-robust covariance cannot see them. That mechanism was **not isolated** by this MC. And the help
gives no att-row magnitude for condition Two at all — its 0.84-1.07 figure there is for **event-study**
standard errors, and the 0.66-0.75 figure at `:251` was measured on the AMESE grid (i = cohort,
6 units per treated cell), a different condition.

**The one cheap test that would settle it:** re-run design B with the cluster-by-period shock removed.
If the ratio returns to ~1, the absorbed-shock reading is confirmed. ~4 minutes.

**What the answer feeds.** Item 4's option (b) wants a measured att-row range at `:1303`; this item is
where that number comes from (0.71-0.81 of the truth on this DGP, if the author accepts documenting it
without the mechanism isolated).

**Recommendation: (b) then (a).** Run the 4-minute shock-removal MC, then correct
`F058_bootstrap.md` §1 and give `sthlp:1303` the att-row range. The sentence the audit proposed for the
F058 record must **not** be inserted: its anchor ("covers 83%.") is unique, but the claim is false.

---

## 7 — `test_level` L11/L12: copied-number pins in legs born inside the review window

**Kind: test policy.** Neither a bug nor a moving number — a question about what the cert suite is
allowed to assert.

**What the legs are.** `test/test_level.do` L11 checks four copied constants at 1e-7
(`.61519121`, `.87494362`, `.63671413`, `.85342071`; file lines 210-211) and L12 four more
(`.07961397`, `.20181734`, `.69718718`, `.22199476`; line 236). Each sits **beside** an identity: L11
also checks `hw/se == invttail(49, .025)` and `invttail(49, .05)`, and that `b`..`p` are unchanged
across `level()`; L12 also checks that `nissize` and `roth_power` are invariant to `level()` while the
slope CI still moves.

**Why it is being raised.** The audit's headline question was whether any re-pin in the 2026-09-05..11
window substituted a copied number for a derivation; the answer is no. But the critique established
that these two legs were **created inside the window**: `test_level.do` was born 2026-09-06
(`a969126`); L11's four numbers arrived in `1acef0d` (BS-5, 2026-09-06) and L12's in `5f5dc6c`
(2026-09-09). They are new legs born as copied-number pins, not re-pins that replaced derivations.

**The critique states the fork exactly** (critique_newfeat_vacuity A6): *"If the coordinator's rule is
'no copied number without a derivation beside it', they comply; if it is 'no copied numbers at all in
legs born after 2026-09-05', they do not."* Both readings are defensible and the project has never
written the rule down.

Two other legs are in the same class and move with the answer: `test_advstat` T2 (3.6886e-06 @ 1e-3,
beside an independent-route reconstruction) and T3 (`cb1996c` bitwise parent pins) — both **predate**
the window (verified: `25727a6` changed only `advt_primary` and T7).

**Recommendation: (a), keep.** A copied number with an identity beside it is a regression pin doing
what a regression pin is for; removing it loses the ability to detect a silent renumber that the
identity would still satisfy.

---

## 8 — The ten UNRULED degeneracy expectations

**Kind: suite curation.** No package number moves under any option; what changes is which config pins
which surface.

`test/suite/degeneracy_expected.csv` was seeded 2026-09-13 from the Stage A degeneracy scan (165 rows).
An `UNRULED` citation means: a degenerate surface that is neither the dataset's design nor a ruled
behaviour, recorded rather than silently accepted.

**Count the UNRULED citations, not the class column.** The file carries **ten** of them — its own seed
comment at line 11 says so ("165 rows, 10 UNRULED citations pending rulings"), and a `grep` for
`UNRULED` returns ten data rows (lines 18, 26, 27, 29, 30, 127, 128, 129, 130, 174) beside two comment
lines. Nine of the ten are class `HIDDEN`; the tenth, `ext_b4_h55_p` at line 127, is class
`PART_RULED`, so anyone who counts `HIDDEN` rows is off by one. That is **8B** below, and the choice
for it is a different choice: its `ADV_POOLED_NOP` stem is already ruled by ledger row `NFADV-ext-p`,
so the only open question is whether the silent `WINDOW_SHORT` loss is acceptable — not whether to keep
or re-host the config. 8B is also the one half of item 8 that is a **code** change rather than suite
curation.

### 8A — the nine HIDDEN rows

| config(s) | stems | what is degenerate | twin that carries the live surface |
|---|---|---|---|
| `simple_esf` | `WINDOW_SHORT` | subgroup(ttre) level 3 (onset t=3) has relyr -2,-1 only, so under `esfixedbaseperiod` the requested pre4/pre3 horizons have no cells; pooled and SG7 are full. **0 cautions** | **none proposed** |
| `ext_b4_h55_f1`, `_f2`, `_f4` (the `_p` sibling is **8B**) | `WINDOW_SHORT` | the onset-6 cohort has 5 pre periods, the ado wants pre6..pre2 under `esfixedbaseperiod`; SG1 is one horizon short. **0 cautions**. The h55 grid corner sits exactly at the data's edge | **none proposed** |
| `trig_tt_year` | `MAIN_INFDEG VES_ABSENT V_ABSENT` | `timetrends` at unit-level `i` under `cluster(i)`: no identified variance on either build | `trig_tt_year_rob` (**robust**, landed at `configs.csv:167`) — not a cluster-below-i twin. The genuine cluster-below-i `timetrends` host on stochastic data is `real_tt` (`i` = st, `cluster(unit)`), which does post a variance |
| `trig_tt_pois` | same | `poisson` x `timetrends`, same geometry | `trig_tt_pois_rob` (**robust**, `configs.csv:168`); the config note records that `cluster(person_id)` is dead too, so robust is the only live host |
| `trig_bigcodes` | `WINDOW_SHORT` | requests `esprelength(2) espostlength(2)` on a YYYYMMDD time grid and realizes **0 pre / 1 post** — `numpre_att` 0, no `pretreatment_*` row, **0 cautions**; the config pins the no-window estimate | `trig_bigids_month` (`configs.csv:169`: same extreme IDs, month-index time, `cluster(i)`, window realized) |
| `trig_sgadv_deeppre` | `PWSGS_DEAD SGDEGEN` | `altmethod_adversarial` x `subgroups`: the pairwise table is blanked by the rank-deficiency screen. The only config with both. Blanked under `robust` too, so it is the design, not the VCE | `ext_b4_h22_f2_alt` — **landed** at `configs_ext.csv:363` |
| `ext_cov10_rcs_twoway_cluster_c2` | `MAIN_INFDEG PWSGS_DEAD SG_DEAD VES_ABSENT V_ABSENT` | the only multi-way clustering config in the suite; b-only on **both** builds, `N_clust` 8, 7 cautions. `sthlp:253` documents multi-way clustering on the default path | `twoway_binary` (`dta_binary`, `cluster(i t)`, probe Q5 `analytic-blocks`) was **proposed and has NOT landed** — `grep -c twoway_binary` = 0 in both config files |

**The package-side observation behind the two silent ones.** `trig_bigcodes` and the `WINDOW_SHORT`
family fire at **zero cautions**: on a non-consecutive integer time grid the requested window is
silently realized as 0/1 horizons with no console message, and horizons with no cells are dropped
without one (the SI-7 caution at ado ~12102 covers margins-omitted horizons, not cell-less ones). That
silence is itself an author's call, separate from this ruling.

**What each option costs.** (a) *keep as pins* costs nothing — every one of the nine agrees on both
builds, so old-vs-new parity is unaffected, and the three landed twins already carry the live surfaces.
(b) *re-host or drop* loses six regression pins, including the only all-zero-VCV caution pin
(`trig_tt_*`) and the only both-builds b-only multi-way pin. (c) *mixed* leaves `simple_esf`, the three
`ext_b4_h55_f*` and `ext_cov10_rcs_twoway_cluster_c2` needing twins that do not exist yet; the
`twoway_binary` row is written and ready to append.

**Recommendation for 8A: (a), keep every one as a pin.** If the author also wants multi-way clustering
shown live, append `twoway_binary`; it is one line and its probe is measured.

### 8B — `ext_b4_h55_p`: the tenth citation, and a different question

The row, verbatim (`test/suite/degeneracy_expected.csv:127`):

```
ext_b4_h55_p,PART_RULED,ADV_POOLED_NOP WINDOW_SHORT,"NFADV-ext-p;UNRULED: WINDOW_SHORT att_s1 pre4/5 - SG1 loses pre5 to the forced reference, silent under the grid's safety(off) (audit F8, critique A4/A5)"
```

**Half of it is already ruled.** The `ADV_POOLED_NOP` stem — no adversarial p on a pooled-reference run
— is cited to ledger row `NFADV-ext-p`, the same unpost that item 10 classifies. Nothing is open there.

**The open half is the `WINDOW_SHORT` stem, and its mechanism is NOT the one the 8A siblings have.**
`ext_b4_h55_p` is the **pooled** arm — its options are `att esprelength(5) espostlength(5)
subgroups(sg_ttre) jointtests adversarial makeplots safety(off)`, with **no** `esfixedbaseperiod`
(scan output `degeneracy_table.csv`, row `ext_b4_h55_p`). Its realized window is
`att=5p5; att_s1=4p5; att_s2=5p5; att_s3=5p5`: the pooled event study realizes pre5..pre1 in full and
**subgroup 1 loses `pre5`**, because with every pre period requested as a horizon the ado **forces the
first observed period back into the reference** for the units that have no panel length left. The
`_f1/_f2/_f4` siblings in 8A lose a horizon for a different reason — theirs is cell-less under
`esfixedbaseperiod` (they realize pre6..pre2) — and `simple_esf` is cell-less too.

**The decisive fact, and it changes the menu: the caution already exists and `safety(off)` silences
it.** `wooldidinstall/wooldid.ado:11594-11598` prints *"Caution: reference period for some units forced
(for want of panel length) to include time period for which treatment effects were requested"* — but
the branch at `:11587-11589`, whose condition is `safety` == off **or** `esprelength` == 0, takes the
silent replace instead, and **every** row of the ext grid carries `safety(off)` (critique_degeneracy_scan A5).
So this is not a missing message; it is a message the grid's own option suppresses. The config's one
recorded caution (`ncaution` = 1, scan output `cautions.csv`) is the **adversarial** pooled-reference
caution, not a window one, so the window loss is indeed unannounced on this run.

**Why the 8A menu does not apply.** Re-hosting is not on the table: the grid corner *is* the test
subject, and the three `f1/f2/f4` siblings in 8A pin the same edge on the non-pooled arms. The only
question is whether a requested horizon disappearing into a forced reference should stay silent when
the user typed `safety(off)`.

**evidence: thin —** two things were not established for this row specifically. (1) Whether v1.3
realizes the same 4-of-5 window: the "both builds agree" reading is recorded for `simple_esf` (the scan
report's old-dump STALE line), not for `ext_b4_h55_p`; the config is `rc0:bV` on both sides with no
ledger row on `mat b`, which is consistent with agreement but is not a check of the window. (2) Whether
`safety(off)` silencing this particular caution was ever a deliberate choice — the ado comment at
`:11586` records only the BS-2 gating fix, not the `safety(off)` disjunct's rationale.

**Recommendation for 8B: (b).** Let the forced-reference caution through even under `safety(off)`: it
reports a fact about the realized window, not a refusal, and `safety(off)` is documented as a bypass for
screens that stop a run, not for disclosures about what was estimated. It costs one branch condition and
moves no number. If the author rules (a), the clause should be rewritten from `UNRULED:` to
`BY_CONSTRUCTION` so the scanner stops counting it and the next reader is not sent looking for a ruling
that was already made.

---

## 9 — H-9: bounded ledger rows firing at half their bound

**Kind: split.** The 18 ceiling appends are straight hygiene and need no ruling. The widenings are
(ii): a harness tolerance is a number.

**The policy.** The ledger's own rule is that a bounded row's tolerance should sit at roughly **10x**
the measured ceiling, so that ordinary run-to-run noise cannot push a row over and a real regression
cannot hide under it.

**Rows firing at more than half their bound** on the 2026-09-11 dual run (301b706 vs 14c30d3;
reproduced independently by the critique, which confirmed H-9 with nothing to amend):

| row | worst fired | bound | ratio |
|---|---|---|---|
| `A8-sgit-s` | 5.8e-7 | 1e-6 | 0.58 |
| `A8-B21-w` | 5.2e-7 | 1e-6 | 0.52 |
| `A8-sgit-w` | 5.1e-7 | 1e-6 | 0.51 |
| `X-est-w` | 9.5e-5 | 2e-4 | 0.47 |
| `X-est-b` | 9.5e-5 | 2e-4 | 0.47 |

**The hard constraint on any widening.** The comparer's own code comment
(`test/suite/compare_oldnew.py:218-227`) explains why a STALE bounded row means the opposite of what it
says — *"Never repair a STALE line by widening the bound until it fires - that is exactly how a genuine
numeric regression gets blessed by a row whose reason column documents a different, smaller effect."* **Never widen a row that has gone STALE** —
widening turns a row that has stopped detecting anything into a row that will never detect anything.

**Two other bounds are wide for a different reason, and item 1 owns them.**
`SI20-P-ext_cont_ext_cfx-cfxatt`, `...cfx_subgroups-cfxatt` and `...cfx_inference-cfxatt` sit at
**1e-1** on every `cfxatt_*` cell (the `b` column unchanged, se/t moving 3.9e-4, `p` tail-amplified to
2.4e-2), and `SI20-P-si20_lat_byt-V` sits at **1e-2**. A 10% bound on the `b` column of a cfx matrix is
not a noise-vs-wrong-answer separator; the honest fix is column-aware patterns, and the underlying
movement is item 1's df factor.

**Caveat on the evidence:** all five ratios are from **one** run. A second dual run could move any of
them either way; that is precisely the argument option (b) rests on.

**Shipping regardless of the ruling:** the 18 bounded rows carrying no measured ceiling get their
ceilings appended (audit §7.3, byte-exact appends listed there), and the 21 rows whose reason is the
bare word "same" get their class named. Neither changes a tolerance.

**Recommendation: (a), accept the tight bounds knowingly** and record the ceilings, so that a future
maintainer reading a 0.5 ratio knows it was looked at. Widening on one run's evidence trades a real
detector for a comfort margin.

---

## 10 — NFADV-01/02/05 classification for the final report

**Kind: (iii) documentation wording.** The classification decides which section of the old-bug report
these three land in; no code moves.

**The disagreement, stated plainly.** Three sources call this an author-approved recalibration:
the ledger rows `NFADV-core-p` and `NFADV-ext-p` (whose reason text says "author-approved"), the
2026-08-23 report (§1.8, §3.5), and v1.3's own help, which documented the planted-zero method at
`sthlp:323`. Two sources introduce it as a review-verified defect: `wooldid_solved_issues.md:416`
records the ruling verb as **"FIX BOTH adversarial pre-trend statistics (primary measured
over-rejecting ~1.4x nominal under realistic heteroskedastic V — not reliably conservative; altmethod
anti-conservative, corrected covariance moves the pinned adv_altb cell from p=1.0e-22 to ~3.7e-6)"**,
and `4fcf127`'s commit message records **"null Monte Carlo 0.04 vs the old statistic's 0.90"**
(`solved_issues:402`).

**What the reproduction isolated** (`stageB_proofs/adversarial/SECTION.md` §S2, §S3; old vs mid
= `7e30e85` = `25727a6^`, which isolates NFADV-01/02/05 from BS-12):
- `esf_adv`: p 3.53039879927e-10 (old) -> 8.29342036543e-11 (mid), the line re-published from
  `[t _cons]` to `[slope intercept_rt0]` with the slope bitwise-invariant (0.16798562157366959 ->
  0.16798562139139653).
- NFADV-05, `cov_advpool`: old posted `adversarial_pval_att` 2.1774111081e-22 on a pooled-reference
  design; new posts none, with the caution. 145 removal cells across the 36 pooled ext configs. The
  argument for defect here is the strongest of the three: **no single omitted reference period exists**
  on a pooled-reference run, so the anchoring point the old statistic fitted through was fabricated.

**One thing that complicates (b).** NFADV-01's replacement primary statistic was itself replaced three
weeks later by BS-12 (`25727a6`, author ruling 2026-09-08: the one-parameter line through the
reference, q = rank(V) - 1). So "author-approved recalibration" describes a statistic that no longer
ships either; the classification is about the 2026-08-19 transition only.

**Both-builds caveat that belongs in the report whichever way this goes.** On the continuous (ame) ext
value cells, ADVBOUND found the recalibrated test oversized at small cells (31.55% at 6 obs/cell under
`robust`; 5.95% only at 50) and v1.3's frozen-line statistic is worse calibrated still (T7 0.90). Those
cells therefore move from one poorly calibrated p to a better-but-gated one, and neither side's ame p
is a validated number.

**Recommendation: (a), review-verified defect for all three.** The ruling verb was "fix both", the
miscalibration was measured before the ruling, and NFADV-05's old output rested on a fabricated anchor.
"Author-approved" describes how the fix was authorised, not what was wrong.

---

## 11 — T-FN1: `sthlp:1295` says `fastnose` posts the coefficient vector; it does not

**Kind: (iii) for the sentence; (ii) only if the author picks (b).**

**The sentence** (`wooldidinstall/wooldid.sthlp:1295`, the `e(b_es_type)` Stored Results entry):

> When NO horizon has an identified variance ({cmd:impvalues}, {cmd:fastnose}, typically
> {cmd:timetrends}, and any run whose estimation returns an all-zero covariance), the covariance is not
> returned but the coefficient vector still is, under a Caution naming what has no input on that run:
> the estimates remain meaningful even where their uncertainty is unavailable.

**The measurement** (Stage C `staggered_refs` probe, §8): `wooldid y i t ttre, att cluster(i)
esprelength(5) espostlength(6) esfixedbaseperiod fastnose` then `confirm matrix e(b_es_att)` returns
**rc 111**, `e(numpre_att)` is empty, and `e(wdid_vcvstate)` reads `absent`. The identical run
**without** `fastnose` returns rc 0 with `numpre_att` 5.

**The help already contradicts itself on two of the four named cases.** `sthlp:276` (impvalues): *"Every
IM `margins` call runs with `nose`, so the event-study b/V capture never fires: `e(b_es_)`, `e(V_es_)`
and `e(numpre_)` are not returned and the honestdid and pretrends handoffs have no inputs."*
`sthlp:283` (permute fast mode): *"`e(b_es_)`, `e(V_es_)`, and `e(numpre_)` are not returned"*. Both
agree with the code and disagree with `:1295`.

**The case `:1295` correctly describes is the fourth one.** The Stage B reproduction of the R4
withdrawal records the new build's own caution on an all-zero-VCV run: *"no horizon of the att event
study has an identified variance, so `e(b_es_att)` posts without `e(V_es_att)`"*
(`stageB_proofs/contsg/SECTION.md` §S2). So the sentence is right for the all-zero/`timetrends` case
and wrong for `fastnose` and `impvalues`.

**The practical consequence, and why it surfaced.** Any Monte Carlo over `fastnose` runs must read the
`ES_att_*` columns of `e(b)`, which **do** post, rather than `e(b_es_att)`, which does not.

**Option (b) is the number-moving one:** making the `fastnose` path post `e(b_es_att)` and
`e(numpre_att)` changes what a `fastnose` run returns, and would then also have to be decided for
`impvalues` and permute fast mode, whose own help paragraphs currently say the opposite.

**Recommendation: (a).** Trim the parenthetical at `:1295` to the all-zero / `timetrends` case; the
other three paragraphs already state the truth. This is a straight documentation fix and needs no
ruling if the author agrees the sentence, not the code, is wrong.

---

## 12 — T-HD2: honestdid does not reliably refuse a degenerate covariance

**Kind: (iii) for the help sentence and the cert leg; (ii) only for option (b).**

**Why it matters.** `honestdid` and `pretrendpower()` consume `e(V_es_*)` directly and inherit its
identification conditions (`sthlp:1303-1305`). The natural assumption — that a degenerate covariance
will make honestdid complain — is false, and the direction of the failure is the quiet one.

**The measurement** (honestdid 1.1.2, Caceres; Stage C `staggered_refs` §1. Probes zeroed one ES
row/column of V and set its diagonal to the `-1e100` sentinel, to 0, or to missing. Every cell below is
a real honestdid call in the session log: probe loop `refs.log:17441-17460`, console output
`refs.log:17470-17560`, return codes and CI rows `results.txt:190-229`):

| V diagonal | row | `delta(rm)` | `delta(sd)` |
|---|---|---|---|
| `-1e100` | a pre row | **rc 0**, "(*) CI open at both endpoints", CI = the grid bounds +/-2.7496 (garbage, but flagged) | rc 3598 |
| `-1e100` | a post row | rc 198 (ECOS) | **rc 0, CI identical to the clean run** |
| 0 | a pre row | **rc 0, no warning at all**, silently different CI | rc 3498 |
| 0 | a post row | **rc 0, CI identical to the clean run** | **rc 0, identical** |
| missing | either | rc 198 | rc 3598 / 3498 |

The decisive pair: `hdrc|sent_zero_row1_rm|0` with CI rows `.6688125360899168648` /
`1.362395906849830585` at M = 0.5 (`results.txt:198-203`) against the clean run's
`.6467940163832528455` / `1.362395906849830585` (`results.txt:147`) — same upper endpoint, lower
endpoint moved 0.022, **no warning, rc 0**.

**The route asymmetry, which inverts the intuition.** `e(V)` (the PS-62 splice) keeps the `-1e100`
sentinel, because the splice skips non-positive-finite diagonals; `e(V_es_att)` carries 0 or missing on
the same row. So in the structural-zero case the **documented** route 2 (`e(V_es_att)`) is the quieter
failure and the spliced `e(V)` is the louder one. A missing diagonal errors on both routes and every
delta; a zeroed **post** row is silently accepted by both.

**So the guard has to be wooldid's own `e(es_ndegen_*)`** (absent == clean), not honestdid's return
code. That scalar is already documented at `sthlp:1295` as "the count of degenerate rows, posted only
when such rows exist — its absence marks a clean covariance"; what is missing is a sentence saying so
**beside the honestdid handoff**, and a cert leg that pins honestdid's behaviour so a version bump is
noticed.

**Related, and already fixed, so no ruling:** all three suite honestdid call sites omitted
`parallel(0)` while the `parallel` package is installed, so each `delta(rm)` call silently spawned four
extra Stata processes — a breach of the 8-instance cap from inside the suite. Fixed in the working tree
on 2026-09-13 (`honestdid_grid.do:126`, `honestdid_crosscheck.do:73`, `identity_extras`).

**Also validated, so no ruling:** the PS-62 spliced `e(V)` route and the documented `e(V_es_att)` route
give **bit-identical** honestdid CI tables under `delta(rm)` (mreldif 0; `delta(sd)` 2.7e-9, which is
honestdid's own FLCI grid search).

**Option (b) is the number-moving one:** withholding `e(numpre_*)` when `e(es_ndegen_*)` is nonzero
would leave the handoff with no inputs at all rather than with inputs the user is warned about.

**Recommendation: (a).** Document the guard and pin honestdid's behaviour in a cert leg; do not change
what wooldid returns. The failure is upstream and the disclosure wooldid already computes is sufficient
once it is pointed at.

---

## 13 — T-BJS1: the `did_imputation` equivalence conditions omit `esprelength(0)`

**Kind: (iii) documentation wording. This is a straight documentation fix and needs no ruling** — it is
listed so the author sees the sentence being changed and can phrase it himself if he wants to.

**The sentence** (`wooldidinstall/wooldid.sthlp:1461-1465`; the conditions are at line **1464**):

> ... given that one uses the {it:att} estimand-type, uses a pooled reference period (i.e., not using
> {it:esfixedbaseperiod}), does not use unconditional standard errors, and does not use controls.

**Why the list is incomplete.** A pooled reference is necessary but not sufficient, because the
reference pool is *the horizons the run asks for*. Asking for pre horizons takes those periods out of
the pool and shifts every post coefficient. Equivalence needs **`esprelength(0)`**.

**The measurement** (`espre0.do`, 2026-09-13, one self-contained run on the `cohom` panel: common onset
t=7, 40 units, T=12, tau=0.5; `did_imputation y i t ttre, horizons(0/5) minn(0) cluster(i)` with and
without `pretrends(4)` against `wooldid y i t ttre, att cluster(i) espostlength(6)` at
`esprelength(0)` and `esprelength(4)`; results at `espre0.txt:36-48`):

| comparison | worst \|absdiff\| over h = 0..5 | verdict line |
|---|---|---|
| wooldid `esprelength(0)` vs BJS | **4.0738e-08** | `EQUIVALENCE_HOLDS` |
| wooldid `esprelength(4)` vs BJS | **4.7619e-02** (the same shift at every horizon) | `EQUIVALENCE_BROKEN` |
| BJS `pretrends(4)` vs BJS without it | **0.0000e+00** | `POST_TAUS_UNCHANGED` |

The third row is what makes this a trap rather than a difference of convention: **asking the reference
tool for pre-trends does not move its post estimates at all, and asking wooldid for pre horizons moves
every one of them.** A user running the two side by side as an event study — the natural thing to do —
sees a uniform 0.0476 gap and has nothing in the help to explain it.

Across the four fixtures: `cohom` 4.76e-2 (constant), `cohet` 4.76e-2, `costag` 1.17e-2 (**not**
constant), `dta_binary` 7.05e-2 (post0 shift 3.71e-2, not constant). At `esprelength(0)` agreement runs
5.1e-8 to 9.8e-8, and over 100 MC draws the worst per-horizon disagreement is 1.16e-07 (`stag`/`post7`)
— the two solvers' FE convergence tolerance, not a conceptual gap.

**The fix:** add `esprelength(0)` (a pooled reference **with no pre horizons requested**) to the
condition list at `sthlp:1464`.

**One adjacent finding that is a record, not a defect** (`T-BJS2`): BJS `pretrends(k)` against
wooldid's pooled-reference pre coefficients is a *conditional* identity — exact to 1e-13 under common
onset, but 0.012-0.029 absolute on `costag` and 0.006-0.049 on `dta_binary` under staggered adoption,
because BJS fits a different regression on the untreated observations. Do not write a Stage C equality
test on pre coefficients outside common onset.
