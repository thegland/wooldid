# hv-fastnose-doc-refute — SKEPTIC VERDICTS

Every number below is from my own data and my own do-files in this directory.
I read the prober's probes but ran none of them.

## Shared design
40 states x 16 years x 20 individuals = 12,800 rows. `i` = state, `t` = year,
`ttre` = adoption (2010 / 2013 / 2016; 12 never-treated states). `dose` = state
policy generosity, constant within state (measured max within-state SD 1.9e-14),
range 0.5-3.0; structural dose effect 0.40. This is the help's own d3 case
(sthlp:902, "a dose constant within each cohort"), on an administrative shape a
real applied economist would have.

---

## HV-FN1 — CONFIRMED (severity low)
The prober's finding is itself a refutation of a recon lead; my run confirms the
refutation. Repro: `p1_fn1_binary.do`.

| | analytic | + fastnose |
|---|---|---|
| `e(b_es_att)` | rc 0 | **rc 0** |
| `e(V_es_att)` | rc 0 | rc 111 |
| `e(V)` | rc 0 | rc 111 |
| `e(numpre_att)` | 3 | **3** |
| `e(wdid_vcvstate)` | analytic-blocks | absent |
| colnames | pre4 pre3 pre2 post0..post3 | identical |

Max abs difference between the two coefficient vectors: **4.441e-16**. The fastnose
arm prints "Caution: no horizon of the att event study has an identified variance,
so e(b_es_att) posts without e(V_es_att)."

So sthlp:1325 is TRUE of `fastnose` on the pinned build and
`known_bydesign.md §7`'s last bullet ("Item 11 looks unfinished ... check whether
fastnose actually posts e(b_es_att)/e(numpre_att); if it does not, this is a live
doc-vs-code divergence") is REFUTED BY MEASUREMENT. It does.

Paper trail I verified independently: the last column of the ruling table in
`review_process/held_for_ruling_2026-09-13.md:49` is headed **"recommendation"**,
not "ruling", so the (a) on item 11 is the recommendation, not the author's choice;
`6fe814d`'s message says "fastnose now posts the event-study coefficient vector ...
the covariance stays absent ... impvalues is removed from the no-identified-variance
list". The file header added in `40d46ac` says item 11 rested on "the commit message
plus partial code evidence" — this run supplies the code evidence.

---

## HV-FN2 — PARTIAL (severity revised high -> low)

**The behaviour reproduces exactly on my own data.** `p2_fn2_d3.do`, arm A2
(`att robust contreat(dose) contreatcontrols(basic) ... fastnose`):
the R6 caution is ABSENT; `main:ame` = 0, `pretreatment:ame` = 0, all seven
`ES_ame` rows print `0` in the table and in `e(wdidmainresults)`; `e(b_es_ame)`
posts with seven exact zeros (`+0.0000000000000X-3ff`); `e(numpre_ame)` = 3;
the only ame message is "no horizon of the ame event study has an identified
variance ... The estimates stand". Arm A1 (analytic `robust`): the full R6 caution
prints, every ame cell is `.`, `e(b_es_ame)` rc 111. The att side is the same
number in both arms (.70128565). Mechanism as stated: ado:11958 gates the R6 screen
on `"`fastnose'" == ""`.

**What is wrong with the finding is its framing, and it is load-bearing.**

1. **The skip is documented — in the option entry for the very option that creates
   the design.** `sthlp:902`, `{opt contreat(varname)}`: "...under default safety
   wooldid names this situation, reports the ame results as missing rather than as
   the regression's structural zeros ... Only the raw e(b) coefficient vector
   retains the structural zeros, since Stata cannot post missing coefficients,
   **and the check is skipped under {it:fastnose}.**" The prober cited sthlp:40 (the
   one-line overview) and sthlp:1325, and missed this. So this is not a
   doc-vs-code divergence; the code does what the option's own paragraph says.

2. **The stated expectation is false.** The finding says these zeros "may reach the
   output only behind an explicit opt-out AND a warning". I measured arm A3,
   `robust safety(off)`: ame rows print **0**, `e(b_es_ame)` posts, and **no caution
   names them** — safety(off) is equally silent, because the sentence that mentions
   safety(off) lives *inside* the R6 caution that safety(off) skips. `sthlp:1209`
   states this in advance: "the ame rows print the regression's structural zeros
   instead of missing values, and without the caution that names them. This is the
   same output an explicit safety(off) produces." fastnose is therefore not
   "inheriting the safety(off) outcome minus the warning" — it is inheriting the
   safety(off) outcome exactly.

3. **The convention is deliberate and repeated.** I probed the sibling screen the
   prober did not: the BS-3 `interactivecontrols()` twin at ado:12005 carries the
   identical gate, and behaves identically (`p3_permute_and_ic.do`: the caution
   "every treated cohort-period slope term for the interactive control urate was
   dropped ... this run is identical to controls(i.t#c.urate)" prints in the analytic
   arm B1, is absent in the fastnose arm B2). And `sthlp:1068` documents THAT one
   too: "(that check does not run under {it:fastnose} or {it:safety(off)}, where the
   reduction is silent)." Two screens, two explicit documentation sentences.

4. **Bounded blast radius.** The fastnose paragraph itself instructs "re-run the
   chosen specification without fastnose before reporting results", and that re-run
   refuses the ame outright. Both arms print the generic multicollinearity caution
   ("triggered 420 times total").

**What survives.** A real UX hazard, at low severity: a spec-search user on a d3
design reads `ame = 0` where the analytic run reads `.`, under a caution that says
"The estimates stand". A point estimate moving from `.` to `0` is also far outside
the tolerance the same help paragraph advertises for fastnose ("absolute differences
on the order of 1e-15"). And `sthlp:1324`'s "A honestdid or pretrends handoff
therefore finds no ame input at all on such a design" reads as false under fastnose
and under safety(off) — though it is scoped to "Where the dose-slope screen voids
the ame block", which under both of those is precisely what does not happen.
Worth a sentence, not a code change.

---

## HV-FN3 — CONFIRMED (severity low, as claimed)
`p2_fn2_d3.do` arm A1: `e(b_es_ame)` rc 111 (withheld under amevoid) while
`e(numpre_ame)` = **3** posts. ado:4148 guards the matrix, ado:4149 does not guard
the scalar. Exactly as claimed. Consequence is a noisy failure (a script that
branches on `e(numpre_ame) < .` then hits a missing matrix), not a wrong number.

---

## Negative results from this pass (nothing found)
* **Per-subgroup event-study covariance under fastnose** (`p5_fn_subgroups.do`).
  ado:1981 keeps a *variance* on the subgroup-by-horizon margins call under
  fastnose, so I checked whether a per-subgroup `e(V_es_att_s*)` leaks into `e()`
  against the help's flat "e(V_es_) is not returned". It does not: analytic arm
  posts b/V/numpre for att, s1 and s2 (rc 0 across the board); fastnose arm posts
  `b_es_att`, `b_es_att_s1`, `b_es_att_s2` and all three `numpre_` (= 3) with
  **every** `V_es_*` rc 111 and `e(wdid_vcvstate)` = absent. The subgroup variance
  is used for the estimability screen and discarded, exactly as documented.
* **The interactivecontrols screen twin** (above): documented at sthlp:1068, so not
  filed as a finding despite reproducing.

---

## The permute() arm (p4_permute_d3.do) — checked because it is the one arm that could
## have escalated FN2, and it does not

`permute(shuffle) permreps(20)` on the same d3 design borrows `fastnose` internally
(ado:5448), so R6 is skipped there too and the ame rows DO reach a reportable table:
`main:ame` = 0.000000, `p_perm` = 1.000000, interval [0.000000, 0.000000], same for
`pretreatment:ame` and all seven `ES_ame` rows. BUT the engine names every one of them:

> Caution: the permutation distribution is degenerate (every completed draw equal) for:
> main_ame pretreatment_ame ES_ame_pre4 ES_ame_pre3 ES_ame_pre2 ES_ame_post0
> ES_ame_post1 ES_ame_post2 ES_ame_post3; those permutation p-values are uninformative.

`e(b_es_ame)` rc 111 and `e(numpre_ame)` empty there, per sthlp:283 (permute fast mode
returns none of the three). So the reportable path is protected by a row-naming caution;
FN2's reach claim is confirmed and its severity is NOT escalated. Plain `fastnose`
remains the worst case, and it is the exploratory mode the help tells you to re-run.

## Dating the documentation (seals "documented, not divergence")
`git log -S'skipped under {it:fastnose}' -- wooldidinstall/wooldid.sthlp` -> `81c284b`
(the coherency slate documentation half, 2026-09-10). The identical query for the
interactivecontrols sentence 'does not run under {it:fastnose}' -> the same commit.
Both sentences predate the pinned build by a week; neither is a post-hoc patch.

## The one thing the author should actually fix (brief category 5, intra-help)
The `fastnose` option paragraph (sthlp:286-297) is the paragraph a fastnose user reads,
and it promises "only point estimates are produced" and that "Point estimates can differ
from a full run in roughly the last digit of precision (absolute differences on the order
of 1e-15)". On a d3 design the difference is `.` -> `0`. The paragraph never mentions the
collinearity screens; the sentence that does lives in the `contreat()` entry at sthlp:902.
That is an internal inconsistency in the help, not a doc-vs-code divergence, and one
clause in the fastnose paragraph closes it.
